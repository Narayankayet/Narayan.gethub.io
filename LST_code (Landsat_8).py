
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

img = gdal.Open(r'H:\Corba\Landsat data (1977 to 2021)\2021\2021_subset_landsat') 

nbands = img.RasterCount
print(nbands)
nrows = img.RasterYSize
print(nrows)
ncols = img.RasterXSize
print(ncols)

print ("band count: " + str(img.RasterCount))


b4 = img.GetRasterBand(4).ReadAsArray().astype(np.float32)

b5 = img.GetRasterBand(5).ReadAsArray().astype(np.float32)

b10 = img.GetRasterBand(10).ReadAsArray().astype(np.float32)

# Step-1: Calculation of TOA (Top of Atmospheric) spectral radiance

# TOA (L) = ML * Qcal + AL
# ML = Band-specific multiplicative rescaling factor from the metadata (RADIANCE_MULT_BAND_x, where x is the band number).
#Qcal = corresponds to band 10.
# AL = Band-specific additive rescaling factor from the metadata (RADIANCE_ADD_BAND_x, where x is the band number).

# RADIANCE_MULT_BAND_10 = 3.3420E-04
# RADIANCE_ADD_BAND_10 = 0.10000

TOA_array = 0.0003342 * b10 + 0.10000

TOA_array[TOA_array== 0.1] = np.nan

L=TOA_array


# Step-2  TOA to Brightness Temperature conversion

# BT = (K2 / (ln (K1 / L) + 1)) − 273.15

#K1 = Band-specific thermal conversion constant from the metadata (K1_CONSTANT_BAND_x, where x is the thermal band number).
# K2 = Band-specific thermal conversion constant from the metadata (K2_CONSTANT_BAND_x, where x is the thermal band number).

# K1_CONSTANT_BAND_10 = 774.8853
#  K2_CONSTANT_BAND_10 = 1321.0789
# L = TOA
# Therefore, to obtain the results in Celsius, the radiant temperature is adjusted by adding the absolute zero (approx. -273.15°C).


BT_array= 1321.0789/(np.log1p(774.8853/L + 1)) -273.15
                        


#Calculate the NDVI

# NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)

ndvi_array=(b5-b4)/(b5+b4)

print(str("---"*20+"\n")+", ".join([
    "ndvi_array stats --- mean: "+str(np.nanmean(ndvi_array)), 
    "std: "+str(np.nanstd(ndvi_array)), 
    "min: "+str(np.nanmin(ndvi_array)), 
    "max: "+str(np.nanmax(ndvi_array)),
    "median: "+ str(np.median(ndvi_array))
]))

# Calculate Prop VEG using NDVI's min max, Use propveg to calculate the LSE:
    

LSE_array = np.square (ndvi_array) +(np.nanmin(ndvi_array)) /(np.nanmax(ndvi_array))- (np.nanmin(ndvi_array))



# NDVI_normalized= ((ndvi_array-np.nanmin(ndvi_array)))/((np.nanmax(ndvi_array)-np.nanmin(ndvi_array))) 


# Calculate Emissivity ε
# E_array= 0.004 *LSE_array + 0.986

E_array = 0.004 *LSE_array +0.986



# Calculate the Land Surface Temperature


LST_array = (BT_array /1 + 0.00115 *BT_array/1.4388 *np.log1p( E_array))


print(str("---"*20+"\n")+", ".join([
    "LST_array stats --- mean: "+str(np.nanmean(LST_array)), 
    "std: "+str(np.nanstd(LST_array)), 
    "min: "+str(np.nanmin(LST_array)), 
    "max: "+str(np.nanmax(LST_array)),
    "median: "+ str(np.median(LST_array))
]))


outraster1 = gdal.GetDriverByName('GTiff').Create(
   (r'D:\SAC project work\Vegetation index\LST_1.tif'),                                    
    img.RasterXSize,                                
    img.RasterYSize,                               
    1 ,                                            
    gdal.GDT_Float32)                               

geo = img.GetGeoTransform()                       
outraster1.SetGeoTransform(geo)                    
wkt = img.GetProjection()                           
outraster1.SetProjection(wkt)                      

outraster1b = outraster1.GetRasterBand(1) 

LST_array[np.isnan(LST_array)] = -9999.           
outraster1b.WriteArray(LST_array)                 
outraster1b.SetNoDataValue(-9999)               


