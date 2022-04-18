

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


os.environ['PROJ_LIB'] =(r'C:\Users\narayan\anaconda3\envs\rsgis\Library\share\proj')

with open(r'D:\ang20180415t061649_rfl_v2s1\ang20180415t061649_corr_v2s1_img.hdr')  as f:
    print("Metadata keys:\n"+", ".join(
        [ln.strip().split(" = ")[0] for ln in f.readlines() if " = " in ln]))

img = gdal.Open(r'D:\Vegetation indeces\Data_Mask') 

nbands = img.RasterCount
nrows = img.RasterYSize
ncols = img.RasterXSize

print("\n".join(["Bands:\t"+str(nbands),"Rows:\t"+str(nrows),"Cols:\t"+str(ncols)]))

