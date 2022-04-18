

import numpy as np
import h5py as h5py
import gdal,osr,os
import matplotlib.pyplot as plt
import warnings
f = h5py.File(r'C:\Users\narayan\Downloads\data\NEON_hyperspectral_tutorial_example_subset.h5','r')
print(f)

def list_dataset(name,node):
    if isinstance(node, h5py.Dataset):
        print(name)

f.visititems(list_dataset)

SJER_refl = f['SJER']['Reflectance']
print(SJER_refl)

SJER_reflArray = SJER_refl['Reflectance_Data']
print(SJER_reflArray)

refl_shape = SJER_reflArray.shape
print('SJER Reflectance Data Dimensions:',refl_shape)

wavelengths = SJER_refl['Metadata']['Spectral_Data']['Wavelength']
print('wavelengths:',wavelengths)

print('min wavelength:', np.amin(wavelengths),'nm')
print('max wavelength:', np.amax(wavelengths),'nm')

SJER_mapInfo = SJER_refl['Metadata']['Coordinate_System']['Map_Info']
print('SJER Map Info:',SJER_mapInfo.value)

mapInfo_string = str(SJER_mapInfo.value)
mapInfo_split = mapInfo_string.split(",") 
print(mapInfo_split)

res = float(mapInfo_split[5]),float(mapInfo_split[6])
print('Resolution:',res)

scaleFactor = SJER_reflArray.attrs['Scale_Factor']
noDataValue = SJER_reflArray.attrs['Data_Ignore_Value']
print('Scale Factor:',scaleFactor)
print('Data Ignore Value:',noDataValue)

b56 = SJER_reflArray[:,:,55].astype(float)
print('b56 type:',type(b56))
print('b56 shape:',b56.shape)
print('Band 56 Reflectance:\n',b56)











b100 = SJER_reflArray[:,:,55].astype(float)
print('b100 type:',type(b56))
print('b100 shape:',b56.shape)
print('Band 100 Reflectance:\n',b100)





b56[b56==int(noDataValue)]=np.nan
b56 = b56/scaleFactor
print('Cleaned Band 56 Reflectance:\n',b56)


SJER_plot = plt.imshow(b56,cmap='Greys') 


SJER_plot = plt.imshow(b56 ,cmap='Greys',clim=(0,0.4)) 
plt.title('SSJER Band 56 Reflectance');


plt.hist(b56[~np.isnan(b56)],10000);
plt.title('Histogram of SJER Band 56 Reflectance')
plt.xlabel('Reflectance'); plt.ylabel('Frequency')



SJER_fig = plt.figure(figsize=(20,10))
ax1 = SJER_fig.add_subplot(1,2,1)

SJER_plot = ax1.imshow(b56,cmap='jet') 
cbar = plt.colorbar(SJER_plot,aspect=50); cbar.set_label('Reflectance')
plt.title('SJER Band 56 Reflectance');
ax1.ticklabel_format(useOffset=False, style='plain') 
rotatexlabels = plt.setp(ax1.get_xticklabels(),rotation=270) 


ax2 = SJER_fig.add_subplot(2,2,2)
ax2.hist(b56[~np.isnan(b56)],50); 
plt.title('Histogram of SERC Reflectance')
plt.xlabel('Reflectance'); plt.ylabel('Frequency')


ax3 = SJER_fig.add_subplot(2,2,4)
ax3.hist(b56[~np.isnan(b56)],50); 
plt.title('Histogram of SJER Reflectance, 0-0.5')
plt.xlabel('Reflectance'); plt.ylabel('Frequency')
ax3.set_xlim([0,0.5])


SJER_fig2 = plt.figure(figsize=(25,25))
ax1 = SJER_fig2.add_subplot(1,3,1)
SJER_plot = ax1.imshow(b56,cmap='gray',clim=(0,0.3)) 
SJER = plt.colorbar(SJER_plot,aspect=50); cbar.set_label('Reflectance')
plt.title('clim = 0-0.3'); 
ax1.ticklabel_format(useOffset=False, style='plain') 
rotatexlabels = plt.setp(ax1.get_xticklabels(),rotation=270) 


ax2 = SJER_fig2.add_subplot(1,3,2)
SJER_plot = ax2.imshow(b56,cmap='gray',clim=(0,0.2)) 
SJER = plt.colorbar(SJER_plot,aspect=50); cbar.set_label('Reflectance')
plt.title('clim = 0-0.2'); 
ax1.ticklabel_format(useOffset=False, style='plain') 
rotatexlabels = plt.setp(ax2.get_xticklabels(),rotation=270) 

x3 = SJER_fig2.add_subplot(1,3,3)
SJER_plot = ax3.imshow(b56,cmap='gray',clim=(0,0.1)) 
cbar = plt.colorbar(SJER_plot,aspect=50); cbar.set_label('Reflectance')
plt.title('clim = 0-0.1'); 
ax1.ticklabel_format(useOffset=False, style='plain') 
rotatexlabels = plt.setp(ax3.get_xticklabels(),rotation=270) 











