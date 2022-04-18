import json
import math
import affine
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from osgeo import gdal,ogr,osr
import csv
import os
import rasterio
from sklearn import preprocessing

import earthpy

# Print ENVI header metadata; Open ENVI image with GDAL 

os.environ['PROJ_LIB'] =(r'C:\Users\narayan\anaconda3\envs\rsgis\Library\share\proj')

with open(r'D:\ang20180415t061649_rfl_v2s1\ang20180415t061649_corr_v2s1_img.hdr')  as f:
    print("Metadata keys:\n"+", ".join(
        [ln.strip().split(" = ")[0] for ln in f.readlines() if " = " in ln]))

img = gdal.Open(r'D:\ang20180415t061649_rfl_v2s1\ang20180415t061649_corr_v2s1_img') 

nbands = img.RasterCount
nrows = img.RasterYSize
ncols = img.RasterXSize

print("\n".join(["Bands:\t"+str(nbands),"Rows:\t"+str(nrows),"Cols:\t"+str(ncols)]))


# Get band information 

band_dictionary = {
    "visible-violet": {'lower': 375, 'upper': 450, 'color': 'violet'},
    "visible-blue": {'lower': 450, 'upper': 485, 'color': 'blue'},
    "visible-cyan": {'lower': 485, 'upper': 500, 'color': 'cyan'},
    "visible-green": {'lower': 500, 'upper': 565, 'color': 'green'},
    "visible-yellow": {'lower': 565, 'upper': 590, 'color': 'yellow'},
    "visible-orange": {'lower': 590, 'upper': 625, 'color': 'orange'},
    "visible-red": {'lower': 625, 'upper': 740, 'color': 'red'},
    "near-infrared": {'lower': 740, 'upper': 1100, 'color': 'gray'},
    "shortwave-infrared": {'lower': 1100, 'upper': 2600, 'color': 'white'}
    }

# Band Informatin

between = lambda wavelength, region: region['lower'] < wavelength <= region['upper']
def classifier(band):
    for region, limits in band_dictionary.items():
        if between(band, limits):
            return(region)


band_numbers = [int(b.split("_")[1]) for b in img.GetMetadata().keys() if b != "wavelength_units"]
band_centers = [float(b.split(" ")[0]) for b in img.GetMetadata().values() if b != "Nanometers"]
em_regions = [classifier(b) for b in band_centers]

bands = pd.DataFrame({ 
    "Band number": band_numbers, 
    "Band center (nm)": band_centers, 
    "EM region": em_regions }, index = band_numbers).sort_index()

bands.head(10)

# bands.head(425).to_csv(r'D:\ang20180415t061649_rfl_v2s1\wavelenght.csv', index = False)


 #  Extract image data at specific locations 


sites = (r"D:\Sample file/Export_Output__project.shp")
os.environ['PROJ_LIB'] =(r'C:\Users\narayan\anaconda3\envs\rsgis\Library\share\proj')

driver = ogr.GetDriverByName("ESRI Shapefile")
permafrost_sites = driver.Open(sites, 0)

site1 = json.loads(permafrost_sites[0].GetFeature(0).ExportToJson())

print(site1)


print("ENVI image WKT: \n"+str(img.GetProjectionRef()))
print("\nShapefile WKT: \n"+str(permafrost_sites[0].GetSpatialRef()))


                           

lyr = permafrost_sites.GetLayer() 
feat = lyr.GetFeature(1)         
geom = feat.GetGeometryRef()      


from_srs = lyr.GetSpatialRef()                                         
to_srs = osr.SpatialReference()                                     
to_srs.ImportFromEPSG(4326)                                           
xytransform = osr.CoordinateTransformation(from_srs,to_srs)            

utm_coordinate_pairs = {}
ll_coordinate_pairs = {}
for feature in lyr:
    geom = feature.GetGeometryRef()                                  
    utm_coordinate_pairs[feature['CLASS_NAME']] = (geom.GetX(), geom.GetY()) 
    geom.Transform(xytransform)                                       
    ll_coordinate_pairs[feature['CLASS_NAME']] = (geom.GetX(), geom.GetY())  


x, y = utm_coordinate_pairs['ROI _1' ]
affine_transform = affine.Affine.from_gdal(*img.GetGeoTransform())     
inverse_transform = ~affine_transform                                  
px, py = inverse_transform * (x, y)                               
px, py = int(px + 0.5), int(py + 0.5)                                 



print( "\n".join(["Site 1 UTM coordinates (x,y): "+"\t"*4+str((x,y)),
                  " are equal to geographic coordinates (lng,lat): \t"+str(ll_coordinate_pairs['ROI _1']),
" and fall within image coordinates (pixel,line):\t"+str((px,py))]) )


image_coordinate_pairs = {site: inverse_transform * pair for site,pair in utm_coordinate_pairs.items()}
print("site id: (col, row)")
for site,coord in image_coordinate_pairs.items(): 
        print(site + ": (" + str(int(coord[0] + 0.5)) + ", " + str(int(coord[1] + 0.5)) +")")     
      
        
      
 
    # bands.head(425).to_csv(r'D:\ang20180415t061649_rfl_v2s1\ROI_excle.csv', index = False)    
      
      
  # Plot the spectral curves for the first two permafrost monitoring sites   
  
# Repeat for site 1:
      
site1name = list(image_coordinate_pairs.keys())[0]
site1xy = list(image_coordinate_pairs.values())[0] 
px, py = int(site1xy[0] + 0.5), int(site1xy[1] + 0.5)

band1_array = img.GetRasterBand(1).ReadAsArray()
print("Band 1 reflectance at site 1: "+str(band1_array[py,px]))

get_pixel = lambda img,band,y,x: img.GetRasterBand(band).ReadAsArray()[y,x]

_bands = bands
_bands[site1name+" reflectance"] = [get_pixel(img,b,py,px) for b in range(1,nbands+1)]
_bands.head(1)
 

# _bands.head(425).to_csv(r'D:\ang20180415t061649_rfl_v2s1\ROI_excle_ref.csv', index = False)   



# Repeat for site 2:
    
    
site2name = list(image_coordinate_pairs.keys())[1]    
site2xy = list(image_coordinate_pairs.values())[1]    
px, py = int(site2xy[0] + 0.5), int(site2xy[1] + 0.5) 

_bands[site2name+" reflectance"] = [get_pixel(img,b,py,px) for b in range(1,nbands+1)]
_bands.head(1)




# Plot the spectral curves with matplotlib:

#matplotlib inline

titlefont = {'fontsize':16,'fontweight':2,
             'verticalalignment':'baseline','horizontalalignment':'center'}
plt.rcParams['figure.figsize'] = [16, 8]


ax1 = plt.subplot(211)
ax2 = plt.subplot(212, sharex=ax1)

_bands.plot(x='Band center (nm)', y=site1name+" reflectance", 
            ax=ax1, c='black', label='_nolegend_', legend=False)
_bands.plot(x='Band center (nm)', y=site2name+" reflectance", 
            ax=ax2, c='black', label='_nolegend_', legend=False)

for i, ax in enumerate([ax1,ax2]): 
    for region,limits in band_dictionary.items():
        ax.axvspan(limits['lower'], limits['upper'], alpha=0.2, 
                   color=limits['color'], label=region)
        

    ax.axvspan(1340, 1445, alpha=0.1, color='green', label='water vapor regions')
    ax.axvspan(1790, 1955, alpha=0.1, color='green')
    

    ax.set_ylim(0,0.6)
    ax.set_xlim(min(band_centers),max(band_centers))
    ax.set_ylabel("reflectance", fontsize=16)
    ax.set_xlabel("wavelength (nm)", fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid('on', alpha=0.25)
    ax.set_title("Jhariya sample site "+str(i+1), fontdict = titlefont, pad = 10)

ax1.legend(prop={'size': 12}, loc = 'upper right', 
           bbox_to_anchor=(1, 1.2), ncol = 2, framealpha = 1)




# Plot individual bands and multiband composites 



get_band_number = lambda w: bands.iloc[(bands["Band center (nm)"]-w).abs().argsort()[1]]

Ri, Gi, Bi = get_band_number(642.32), get_band_number(552.16), get_band_number(472.02)


print(str("\n"+"---"*20+"\n").join([str(Ri),str(Gi),str(Bi)]))

get_band = lambda b: img.GetRasterBand(int(b["Band number"])).ReadAsArray()


Ra, Ga, Ba = get_band(Ri), get_band(Gi), get_band(Bi)



Ra[Ra == -9999.], Ga[Ga == -9999.], Ba[Ba == -9999.] = 0, 0, 0

print(Ra,Ga, Ba)


scale8bit = lambda a: ((a - a.min()) * (1/(a.max() - a.min()) * 255)).astype('uint8')
Ra8, Ga8, Ba8 = scale8bit(Ra), scale8bit(Ga), scale8bit(Ba)


Ra8[Ra == 0], Ga8[Ga == 0], Ba8[Ba == 0] = 0, 0, 0

Ra8 = np.ma.masked_equal(Ra8, 0)
Ga8 = np.ma.masked_equal(Ga8, 0)
Ba8 = np.ma.masked_equal(Ba8, 0)


rgb_stack = np.zeros((nrows,ncols,3),'uint8')
rgb_stack[...,0], rgb_stack[...,1], rgb_stack[...,2] = Ra8, Ga8, Ba8

print("---"*20+"\nRGB stack (rows,cols,bands): " + str(rgb_stack.shape))

xmin, xres, xrot, ymax, yrot, yres = img.GetGeoTransform()
xarr = np.array([int(xmin+i*xres) for i in range(0,ncols)])
yarr = np.array([int(ymax+i*yres) for i in range(0,nrows)])

print("the first 10 x coordinates:")
xarr[:10]



# Plot the R,G,B bands and the RGB stack:
    
 
plt.rcParams['figure.figsize'] = [20, 20]
gs = gridspec.GridSpec(1, 4)

plotdict = { 'Red': { 'subplot': 0, 'array': Ra8, 'colormap': 'Reds_r' },
             'Green': { 'subplot': 1, 'array': Ga8, 'colormap': 'Greens_r' },
             'Blue': { 'subplot': 2, 'array': Ba8, 'colormap': 'Blues_r' },
             'RGB': { 'subplot': 3, 'array': rgb_stack, 'colormap': None } }

fig1 = plt.figure()
for band,data in plotdict.items():
    clim = None if band == "RGB" else (0,255)
    ax = fig1.add_subplot(gs[ 0, data['subplot'] ])
    p = ax.imshow(data['array'], cmap=data['colormap'], clim=clim,
                  extent=[xmin,xmin+ncols*xres,ymax,ymax+nrows*yres])
    ax.set_title(band, pad = 20, fontdict = titlefont)  
    
plt.imshow(rgb_stack, extent=[xmin,xmin+ncols*xres,ymax,ymax+nrows*yres])



# Equalize histograms ¶

for i in range(rgb_stack.shape[2]):
 b = rgb_stack[:,:,i]

 b_histogram, bins = np.histogram(b.flatten(), 256)
b_cumdistfunc = b_histogram.cumsum()
b_cumdistfunc = 255 * b_cumdistfunc / b_cumdistfunc[-1]
b_equalized = np.interp(b.flatten(), bins[:-1], b_cumdistfunc)

rgb_stack[:,:,i] = b_equalized.reshape(b.shape)


plt.rcParams['figure.figsize'] = [20, 20]
gs = gridspec.GridSpec(1, 4)

plotdict = { 'Red': { 'subplot': 0, 'array': rgb_stack[:,:,0], 'colormap': 'Reds_r' },
             'Green': { 'subplot': 1, 'array': rgb_stack[:,:,1], 'colormap': 'Greens_r' },
             'Blue': { 'subplot': 2, 'array': rgb_stack[:,:,2], 'colormap': 'Blues_r' },
             'RGB': { 'subplot': 3, 'array': rgb_stack, 'colormap': None } }



fig1 = plt.figure()
for band,data in plotdict.items():
    clim = None if band == "RGB" else (0,255)
    ax = fig1.add_subplot(gs[ 0, data['subplot'] ])
    p = ax.imshow(data['array'], cmap=data['colormap'], clim=clim,
                  extent=[xmin,xmin+ncols*xres,ymax,ymax+nrows*yres])
    ax.set_title(band, pad = 20, fontdict = titlefont)
plt.imshow(rgb_stack)


# NDVI

ndvi_red, ndvi_nir = get_band_number(685), get_band_number(900)

print(str("\n"+"---"*30+"\n").join([str(ndvi_red),str(ndvi_nir)]))


R685, R900 = get_band(ndvi_red), get_band(ndvi_nir)

R685[R685 == -9999.], R900[R900 == -9999.] = np.nan, np.nan


ndvi_array =(R900-R685)/(R900+R685)

ndvi_array[ndvi_array>=1], ndvi_array[ndvi_array<=-1] = np.nan, np.nan

print(str("---"*20+"\n")+", ".join([
    "ndvi_array stats --- mean: "+str(np.nanmean(ndvi_array)), 
    "std: "+str(np.nanstd(ndvi_array)), 
    "min: "+str(np.nanmin(ndvi_array)), 
    "max: "+str(np.nanmax(ndvi_array)),
    "median: "+ str(np.median(ndvi_array))
]))
        
NDVI_normalized= ((ndvi_array-np.nanmin(ndvi_array)))/((np.nanmax(ndvi_array)-np.nanmin(ndvi_array))) 



plt.savefig('filename.tif', dpi=500)
plt.imshow(NDVI_normalized,extent=[-1,1,-1,1])
plt.imshow(NDVI_normalized, cmap='RdYlGn'), 
plt.colorbar(label="Index Values ( 0 to 1)", orientation="vertical")
plt.title('ARI_1 Index {}'.format(0))





plt.xlabel('Flosh ')
plt.ylabel('Flosh')




#sd = np.std(ndvi_array)
#mean = np.mean(ndvi_array)
#outliers = (ndvi_array > (mean + 2 * sd)) | (ndvi_array < (mean - 2 * sd))
#ndvi_array[outliers] = np.nan""

#plt.imshow(ndvi_array1)
# plt.colorbar(ndvi_array)



# SR

R750, R705 = get_band(get_band_number(750)), get_band(get_band_number(705))

R750[R750 == -9999.], R705[R705 == -9999.] = np.nan, np.nan

sr_array = R750/R705

sr_array[sr_array>=30], sr_array[sr_array<=0] = np.nan, np.nan

SR_normalized= ((sr_array-np.nanmin(sr_array)))/((np.nanmax(sr_array)-np.nanmin(sr_array)))



print(str("---"*20+"\n")+", ".join([
    "SR_normalized  stats --- mean: "+str(np.nanmean(SR_normalized)), 
    "std: "+str(np.nanstd(SR_normalized)), 
    "min: "+str(np.nanmin(SR_normalized)), 
    "max: "+str(np.nanmax(SR_normalized)) 
]))







# SAVI

savi_array = 0.5*(R900-R685)/(R900+R685+0.5)

savi_array[savi_array == -9999.] =np.nan

print(str("---"*20+"\n")+", ".join([
    "Ssavi_array stats --- mean: "+str(np.nanmean(savi_array)), 
    "std: "+str(np.nanstd(savi_array)), 
    "min: "+str(np.nanmin(savi_array)), 
    "max: "+str(np.nanmax(savi_array)) 
]))

SAVI_normalized= ((savi_array-np.nanmin(savi_array)))/((np.nanmax(savi_array)-np.nanmin(savi_array)))

print(str("---"*20+"\n")+", ".join([
    "SAVI normalized stats --- mean: "+str(np.nanmean(SAVI_normalized)), 
    "std: "+str(np.nanstd(SAVI_normalized)), 
    "min: "+str(np.nanmin(SAVI_normalized)), 
    "max: "+str(np.nanmax(SAVI_normalized)) 
]))

# VARI

VARI_red, VARI_green,VARI_blue = get_band_number(665), get_band_number(560),get_band_number(490),

print(str("\n"+"---"*30+"\n").join([str(VARI_red),str(VARI_green),str(VARI_blue)]))

R665, R560 ,R490 = get_band(VARI_red), get_band(VARI_green),get_band(VARI_blue)

VARI_array=(R560-R665)/(R560+R665-R490)

#VARI_array[VARI_array == -9999.] =np.nan

VARI_array[VARI_array== -00000.] = np.nan

# VARI_array[VARI_array>=1], VARI_array[VARI_array<=-1] = np.nan, np.nan


print(str("---"*20+"\n")+", ".join([
    "VARI_array stats --- mean: "+str(np.nanmean(VARI_array)), 
    "std: "+str(np.nanstd(VARI_array)), 
    "min: "+str(np.nanmin(VARI_array)), 
    "max: "+str(np.nanmax(VARI_array)) 
]))


VARI_normalized= ((VARI_array-np.nanmin(VARI_array)))/((np.nanmax(VARI_array)-np.nanmin(VARI_array)))





# RVSI

R714, R752,R733 = get_band_number(714), get_band_number(752),get_band_number(733),
print(str("\n"+"---"*30+"\n").join([str(R714),str(R752),str(R733)]))

R714, R752 ,R733 = get_band(R714), get_band(R752),get_band(R733)

RVSI_array=(R714 + R752/2) -R733

RVSI_array[RVSI_array== -4999.5] = np.nan

print(str("---"*20+"\n")+", ".join([
    "RVSI_array stats --- mean: "+str(np.nanmean(RVSI_array)), 
    "std: "+str(np.nanstd(RVSI_array)), 
    "min: "+str(np.nanmin(RVSI_array)), 
    "max: "+str(np.nanmax(RVSI_array)) 
]))


# RVSI_array[RVSI_array>=1], RVSI_array[RVSI_array<=-1] = np.nan, np.nan


RVSI_normalized= ((RVSI_array-np.nanmin(RVSI_array)))/((np.nanmax(RVSI_array)-np.nanmin(RVSI_array)))




#  CI-red_edge

R705,R842 = get_band_number(705), get_band_number(842)

print(str("\n"+"---"*30+"\n").join([str(R705),str(R842)]))

R705, R842 = get_band(R705), get_band(R842)


CI_array= R842/R705-1

CI_array[CI_array== 0.] = np.nan     

print(str("---"*20+"\n")+", ".join([
    "CI_array stats --- mean: "+str(np.nanmean(CI_array)), 
    "std: "+str(np.nanstd(CI_array)), 
    "min: "+str(np.nanmin(CI_array)), 
    "max: "+str(np.nanmax(CI_array)) 
]))


CI_array[CI_array>=2], CI_array[CI_array<=-2] = np.nan, np.nan


CI_normalized= ((CI_array-np.nanmin(CI_array)))/((np.nanmax(CI_array)-np.nanmin(CI_array)))


# NDLI

R1754,R1680 = get_band_number(1754), get_band_number(1680)

print(str("\n"+"---"*30+"\n").join([str(R1754),str(R1680)]))

R1754, R1680 = get_band(R1754), get_band(R1680)

R1754[R1754 == -9999] = np.nan

R1680[R1680 == -9999] =np.nan

NDLI_array = (np.log1p(R1754)-np.log1p(R1680))/(np.log1p(R1754)+np.log1p(R1680))

print(str("---"*20+"\n")+", ".join([
    "NDLI_array stats --- mean: "+str(np.nanmean(NDLI_array)), 
    "std: "+str(np.nanstd(NDLI_array)), 
    "min: "+str(np.nanmin(NDLI_array)), 
    "max: "+str(np.nanmax(NDLI_array)) 
]))


NDLI_array[NDLI_array>=0.10], NDLI_array[NDLI_array<=-0.20] = np.nan, np.nan


NDLI_normalized= ((NDLI_array-np.nanmin(NDLI_array)))/((np.nanmax(NDLI_array)-np.nanmin(NDLI_array)))


# WBI

R970,R902 = get_band_number(970), get_band_number(902)

print(str("\n"+"---"*30+"\n").join([str(970),str(902)]))

R970, R902 = get_band(R970), get_band(R902)

WBI_array = R970/R902

WBI_array[WBI_array == 1] = np.nan


print(str("---"*20+"\n")+", ".join([
    "WBI_array stats --- mean: "+str(np.nanmean(WBI_array)), 
    "std: "+str(np.nanstd(WBI_array)), 
    "min: "+str(np.nanmin(WBI_array)), 
    "max: "+str(np.nanmax(WBI_array)) 
]))

WBI_array[WBI_array>=2], WBI_array[WBI_array<=-0.03] = np.nan, np.nan

WBI_normalized= ((WBI_array-np.nanmin(WBI_array)))/((np.nanmax(WBI_array)-np.nanmin(WBI_array)))

# CARI

R700,R670,R550 = get_band_number(700), get_band_number(670), get_band_number(550)

print(str("\n"+"---"*30+"\n").join([str(700),str(670),str(550)]))

R700,R670,R550 = get_band(R700), get_band(R670), get_band(R550)

CARI_array =(R700-R670)-0.2*(R700-R550)

CARI_array[CARI_array == 0.0] = np.nan


print(str("---"*20+"\n")+", ".join([
    "CARI_array stats --- mean: "+str(np.nanmean(CARI_array)), 
    "std: "+str(np.nanstd(CARI_array)), 
    "min: "+str(np.nanmin(CARI_array)), 
    "max: "+str(np.nanmax(CARI_array)) 
]))

CARI_array[CARI_array>=0.26], CARI_array[CARI_array<=-0.20] = np.nan, np.nan

CARI_normalized= ((CARI_array-np.nanmin(CARI_array)))/((np.nanmax(CARI_array)-np.nanmin(CARI_array)))

# PSRI

R680,R500,R750 = get_band_number(680), get_band_number(500), get_band_number(750)

print(str("\n"+"---"*30+"\n").join([str(680),str(500),str(750)]))

R680,R500,R750 = get_band(R680), get_band(R500), get_band(R750)

PSRI_array =(R680-R500)/R750

PSRI_array[PSRI_array ==-0] = np.nan

print(str("---"*20+"\n")+", ".join([
    "PSRI_array stats --- mean: "+str(np.nanmean(PSRI_array)), 
    "std: "+str(np.nanstd(PSRI_array)), 
    "min: "+str(np.nanmin(PSRI_array)), 
    "max: "+str(np.nanmax(PSRI_array)) 
]))


PSRI_array[PSRI_array>=1], PSRI_array[PSRI_array<=-0.10] = np.nan, np.nan

PSRI_normalized= ((PSRI_array-np.nanmin(PSRI_array)))/((np.nanmax(PSRI_array)-np.nanmin(PSRI_array)))

# CAI

R2000,R2200,R2100= get_band_number(2000), get_band_number(2200), get_band_number(2100)

print(str("\n"+"---"*30+"\n").join([str(2000),str(2200),str(2100)]))

R2000,R2200,R2100 = get_band(R2000), get_band(R2200), get_band(R2100)


CAI_array =0.5*(R2000+R2200)-R2100

CAI_array[CAI_array ==0] = np.nan


print(str("---"*20+"\n")+", ".join([
    "CAI_array stats --- mean: "+str(np.nanmean(CAI_array)), 
    "std: "+str(np.nanstd(CAI_array)), 
    "min: "+str(np.nanmin(CAI_array)), 
    "max: "+str(np.nanmax(CAI_array)) 
]))

CAI_array[CAI_array>=1], CAI_array[CAI_array<=-0.20] = np.nan, np.nan

CAI_normalized= ((CAI_array-np.nanmin(CAI_array)))/((np.nanmax(CAI_array)-np.nanmin(CAI_array)))


# PSND


R800,R675,R650= get_band_number(800), get_band_number(675), get_band_number(650)

print(str("\n"+"---"*30+"\n").join([str(800),str(675),str(650)]))

R800,R675,R650 = get_band(R800), get_band(R675), get_band(R650)



PSND_array= R800-R675/R800-R675; R800-R650/R800+R650


PSND_array[PSND_array ==-1] = np.nan

print(str("---"*20+"\n")+", ".join([
    "PSND_array stats --- mean: "+str(np.nanmean(PSND_array)), 
    "std: "+str(np.nanstd(PSND_array)), 
    "min: "+str(np.nanmin(PSND_array)), 
    "max: "+str(np.nanmax(PSND_array)) 
]))

PSND_array[PSND_array>=1], PSND_array[PSND_array<=-2] = np.nan, np.nan

PSND_normalized= ((PSND_array-np.nanmin(PSND_array)))/((np.nanmax(PSND_array)-np.nanmin(PSND_array)))


# PSRI


R680,R500,R750= get_band_number(680), get_band_number(500), get_band_number(750)


print(str("\n"+"---"*30+"\n").join([str(680),str(500),str(750)]))

R680,R500,R750 = get_band(R680), get_band(R500), get_band(R750)

PSRI_array= R680-R500/R750

PSRI_array[PSRI_array ==-10000] = np.nan

print(str("---"*20+"\n")+", ".join([
    "PSRI_array stats --- mean: "+str(np.nanmean(PSRI_array)), 
    "std: "+str(np.nanstd(PSRI_array)), 
    "min: "+str(np.nanmin(PSRI_array)), 
    "max: "+str(np.nanmax(PSRI_array)) 
]))


PSRI_array[PSRI_array>=1], PSRI_array[PSRI_array<=-2] = np.nan, np.nan

PSRI_normalized= ((PSRI_array-np.nanmin(PSRI_array)))/((np.nanmax(PSRI_array)-np.nanmin(PSRI_array)))



# NDWI


R857,R1241= get_band_number(857), get_band_number(1241)

print(str("\n"+"---"*30+"\n").join([str(857),str(1241)]))

R857,R1241 = get_band(R857), get_band(R1241)

NDWI_array= R857-R1241/R857+R1241

NDWI_array[NDWI_array ==-19999] = np.nan 

print(str("---"*20+"\n")+", ".join([
    "NDWI_array stats --- mean: "+str(np.nanmean(NDWI_array)), 
    "std: "+str(np.nanstd(NDWI_array)), 
    "min: "+str(np.nanmin(NDWI_array)), 
    "max: "+str(np.nanmax(NDWI_array)) 
]))

NDWI_array[NDWI_array>=1], NDWI_array[NDWI_array<=-3] = np.nan, np.nan


NDWI_normalized= ((NDWI_array-np.nanmin(NDWI_array)))/((np.nanmax(NDWI_array)-np.nanmin(NDWI_array)))

# RGRI


R665,R560= get_band_number(665), get_band_number(560)

print(str("\n"+"---"*30+"\n").join([str(665),str(560)]))

R665,R560 = get_band(R665), get_band(R560)

RGRI_array= R665/R560

RGRI_array[RGRI_array ==1] = np.nan

print(str("---"*20+"\n")+", ".join([
    "RGRI_array stats --- mean: "+str(np.nanmean(RGRI_array)), 
    "std: "+str(np.nanstd(RGRI_array)), 
    "min: "+str(np.nanmin(RGRI_array)), 
    "max: "+str(np.nanmax(RGRI_array)) 
]))

RGRI_array[RGRI_array>=2], RGRI_array[RGRI_array<=0] = np.nan, np.nan


RGRI_normalized= ((RGRI_array-np.nanmin(RGRI_array)))/((np.nanmax(RGRI_array)-np.nanmin(RGRI_array)))

# PSSR

R800,R675,R650= get_band_number(800), get_band_number(675),get_band_number(650)

print(str("\n"+"---"*30+"\n").join([str(800),str(675),str(650)]))

R800,R675,R650 = get_band(R800), get_band(R675),get_band(R650)

PSSR_array= (R800/R675);(R800/R650)

PSSR_array[PSSR_array ==1] = np.nan

print(str("---"*20+"\n")+", ".join([
    "PSSR_array stats --- mean: "+str(np.nanmean(PSSR_array)), 
    "std: "+str(np.nanstd(PSSR_array)), 
    "min: "+str(np.nanmin(PSSR_array)), 
    "max: "+str(np.nanmax(PSSR_array)) 
]))

PSSR_array[PSSR_array>=9], PSSR_array[PSSR_array<=0] = np.nan, np.nan


PSSR_normalized= ((PSSR_array-np.nanmin(PSSR_array)))/((np.nanmax(PSSR_array)-np.nanmin(PSSR_array)))

# EVI

R900,R685,R490= get_band_number(900), get_band_number(685),get_band_number(490)

print(str("\n"+"---"*30+"\n").join([str(900),str(685),str(490)]))

R900,R685,R490 = get_band(R900), get_band(R685),get_band(R490)

EVI_array= 2.5*(R900-R685)/(R900+6* R685-7.5*R490+1)

EVI_array[EVI_array ==0] = np.nan


print(str("---"*20+"\n")+", ".join([
    "EVI_array stats --- mean: "+str(np.nanmean(EVI_array)), 
    "std: "+str(np.nanstd(EVI_array)), 
    "min: "+str(np.nanmin(EVI_array)), 
    "max: "+str(np.nanmax(EVI_array)) 
]))

EVI_array[EVI_array>=1], EVI_array[EVI_array<=0] = np.nan, np.nan

EVI_normalized= ((EVI_array-np.nanmin(EVI_array)))/((np.nanmax(EVI_array)-np.nanmin(EVI_array)))


 # GCI

R842,R560= get_band_number(842), get_band_number(560)

print(str("\n"+"---"*30+"\n").join([str(842),str(560)]))

R842, R560=  get_band(R842), get_band(R560)

GCI_array= (R842/R560)-1

GCI_array[GCI_array == 0] = np.nan


print(str("---"*20+"\n")+", ".join([
    "GCI_array stats --- mean: "+str(np.nanmean(GCI_array)), 
    "std: "+str(np.nanstd(GCI_array)), 
    "min: "+str(np.nanmin(GCI_array)), 
    "max: "+str(np.nanmax(GCI_array)) 
]))


GCI_array[GCI_array>=6.45], GCI_array[GCI_array<=0] = np.nan, np.nan

GCI_normalized= ((GCI_array-np.nanmin(GCI_array)))/((np.nanmax(GCI_array)-np.nanmin(GCI_array)))


# LAI

LAI_array=(3.618*EVI_array-0.118)

print(str("---"*20+"\n")+", ".join([
    "LAI_array stats --- mean: "+str(np.nanmean(GCI_array)), 
    "std: "+str(np.nanstd(LAI_array)), 
    "min: "+str(np.nanmin(LAI_array)), 
    "max: "+str(np.nanmax(LAI_array)) 
]))

LAI_array[LAI_array>=3], LAI_array[LAI_array<=0] = np.nan, np.nan

LAI_normalized= ((LAI_array-np.nanmin(LAI_array)))/((np.nanmax(LAI_array)-np.nanmin(LAI_array)))





# GLI

R560,R665,R490= get_band_number(560), get_band_number(665),get_band_number(490)

print(str("\n"+"---"*30+"\n").join([str(560),str(665),str(490)]))

R560,R665,R490=  get_band(R560), get_band(R665),get_band(R490)

GLI_array= (R560-R665)+(R560-R490)/(2*R560)+R665+R490


GLI_array[GLI_array == -19998] = np.nan


print(str("---"*20+"\n")+", ".join([
    "GLI_array stats --- mean: "+str(np.nanmean(GLI_array)), 
    "std: "+str(np.nanstd(GLI_array)), 
    "min: "+str(np.nanmin(GLI_array)), 
    "max: "+str(np.nanmax(GLI_array)) 
]))


GLI_array[GLI_array>=1], GLI_array[GLI_array<=0] = np.nan, np.nan

GLI_normalized= ((GLI_array-np.nanmin(GLI_array)))/((np.nanmax(GLI_array)-np.nanmin(GLI_array)))


#mARI


R530,R570, R690, R710,  R1400= get_band_number(530), get_band_number(570),get_band_number(690), get_band_number(710),get_band_number(1400)

print(str("\n"+"---"*30+"\n").join([str(530),str(570),str(690),str(710),str(1400)]))

R530,R570, R690, R710, R1400 =  get_band(R530), get_band(R570),get_band(R690), get_band(R710),get_band(R1400)

mARI_array=(1/R530-R570-1/R690-R710)* R1400

mARI_array[mARI_array == -1.9996e+08] = np.nan


print(str("---"*20+"\n")+", ".join([
    "mARI_array stats --- mean: "+str(np.nanmean(mARI_array)), 
    "std: "+str(np.nanstd(mARI_array)), 
    "min: "+str(np.nanmin(mARI_array)), 
    "max: "+str(np.nanmax(mARI_array)) 
]))



# mARI_array[mARI_array>=1], mARI_array[mARI_array<=-3] = np.nan, np.nan

mARI_normalized= ((mARI_array-np.nanmin(mARI_array)))/((np.nanmax(mARI_array)-np.nanmin(mARI_array)))


# MSI

R1599,R819= get_band_number(1599), get_band_number(819)

R1599,R819=  get_band(R1599), get_band(R819)

MSI_array= (R1599/R819)


MSI_array[MSI_array == 1] = np.nan

MSI_normalized= ((MSI_array-np.nanmin(MSI_array)))/((np.nanmax(MSI_array)-np.nanmin(MSI_array)))


#IRECI

R832,R664,R704,R782= get_band_number(832), get_band_number(664),get_band_number(704),get_band_number(782)

R832,R664,R704,R782=  get_band(R832), get_band(R664),get_band(R704),get_band(R782)

IRECI_array= (R832-R664)/(R704/R782)

IRECI_array[IRECI_array ==0] = np.nan

IRECI_normalized= ((IRECI_array-np.nanmin(IRECI_array)))/((np.nanmax(IRECI_array)-np.nanmin(IRECI_array)))

















# plot. all of this is matplotlib


gs = gridspec.GridSpec(1, 3)
plotdict1 = { 'NDVI\n(R900-R685)/(R900+R685)': { 'subplot': 0, 'array': NDVI_normalized },
             'SR\nR750/R705': { 'subplot': 1, 'array': SR_normalized },
             'SAVI\nL*(R900-R685)/(R900+R685+L)': { 'subplot': 2, 'array': SAVI_normalized } }

print(plotdict1)

fig2 = plt.figure()
for band,data in plotdict1.items():
    ax = fig2.add_subplot(gs[ 0, data['subplot'] ])
    p = ax.imshow(data['array'], cmap=plt.get_cmap("RdYlGn"),
                  extent=[xmin,xmin+ncols*xres,ymax,ymax+nrows*yres])
    ax.set_title(band, pad = 20, fontdict = titlefont)
    plt.colorbar(p)



# Save outputs; Common raster operations

outraster1 = gdal.GetDriverByName('GTiff').Create(
   (r'D:\SAC project work\Vegetation index\WBI_normalized_11111.tif'),                                    
    img.RasterXSize,                                
    img.RasterYSize,                               
    1 ,                                            
    gdal.GDT_Float32)                               

geo = img.GetGeoTransform()                       
outraster1.SetGeoTransform(geo)                    
wkt = img.GetProjection()                           
outraster1.SetProjection(wkt)                      

outraster1b = outraster1.GetRasterBand(1) 


WBI_normalized[np.isnan(WBI_normalized)] = -9999.           
outraster1b.WriteArray(WBI_normalized)                 
outraster1b.SetNoDataValue(-9999)                   

outraster1.FlushCache()                             
outraster1 =None                            




Study_shape= (r"D:\ang20180415t061649_rfl_v2s1\Study area shape\Shape_0_roation.shp")
os.environ['PROJ_LIB'] =(r'C:\Users\narayan\anaconda3\envs\rsgis\Library\share\proj')

shp = ogr.Open(Study_shape)                              
lyr = shp.GetLayer()   




outraster2 = gdal.Warp(
    "savi_clipped.tif",                                    
    "SAVI_normalized.tif",                                            
    format = 'GTiff',                                    
    cutlineDSName = Study_shape,                           
    cutlineLayer = lyr.GetName(),                         
    cropToCutline = True,                                
    dstNodata = -9999)                                   













