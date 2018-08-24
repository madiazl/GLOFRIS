
#function to open the files according the different format type

import numpy as np
import netCDF4 as nc
from scipy.interpolate import interp1d
import pdb
# Keep track of computation time
from datetime import datetime 
tstart = datetime.now()
print('asasasasasasssss')
#from osgeo import gdal
# Get all the needed modules
import numpy as np
import os
import csv
import pandas as pd
import subprocess
#from pcraster import *


#gdal.AllRegister()
#driver = gdal.GetDriverByName('Geo')
#driver.Register()

def fun_read_rasters(options_format,options_file_path,options_name_var):
    if (options_format=='.nc'):
        print(options_format)
        print(options_file_path)
        file = nc.Dataset(options_file_path, 'r', format='NETCDF4');
        variable=file.variables[options_name_var]
        print "value", variable._FillValue
        value=np.int(variable._FillValue)
        variable.set_auto_maskandscale(False)
        if (len(variable.shape)==3):                                     #3 dimensional arrays coming from edwin
            variable=variable[0,:,:]
        print variable.shape        
        variable=np.array(variable)
        variable[variable[:][:]==value]=0
        if (file.variables['lat'][0]<0):                            #in case the variable is upside down
            variable=variable[::-1]  # Upside down rasters
        file.close()
    if (options_format=='.map'):
        Dem = readmap(options_file_path)
        variable=pcr2numpy(Dem, -9999) #second value correspon to the missing values and numpy2pcr(dataType, array, mv) to other way around
    if (options_format=='.h5'):
        print('hasta el h5')
        file = nc.Dataset(options_file_path, 'r', format='h5');
        variable=file.variables[options_name_var]
        variable.set_auto_maskandscale(False)
        #if is close the file before extract the array it generates problems
        variable=variable[:,:]
        file.close()
    # if (options_format=='.tif'):
    #     file = gdal.Open(options_file_path)
    #     #print file.GetMetadata()  #retreive metadat in the file
    #     #print file.RasterCount    #to get the number of bands in th file
    #     band_file = file.GetRasterBand(1)
    #     #print band_file.GetStatistics( True, True )
    #     variable = (band_file.ReadAsArray(0, 0, file.RasterXSize, file.RasterYSize).astype(np.float64))
    #     variable[variable==band_file.GetNoDataValue()]=0  #if the non data value is not weird as a float number it works, other wise need to be implemented a trick
          #print(band_file.GetNoDataValue())
          #print(np.nansum(variable))
    return variable;

