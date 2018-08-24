# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 08:11:19 2015

@author: winsemi
"""

### EXAMPLE FOR GFP COURSE ON NETCDF
import netCDF4 as nc
import numpy as np

def prepare_nc(nc_file,x, y, time_list,RP, metadata, units='Days since 2012-01-01 00:00:00', calendar='gregorian',Format="NETCDF4",zlib=True, clim=False):
    """
    This function prepares a NetCDF file with given metadata, for a certain year, daily basis data
    The function assumes a gregorian calendar and a time unit 'Days since 2012-01-01 00:00:00'
    inputs:
        trg_file:     path to new netcdf file
        time_list:    list with times (in datetime format) which will give the time axis
        x:            xaxis
        y:            yaxis
        metadata:     dictionary with global attributes
        units:        time units to use in time axis
    """
    
    print('Setting up "' + nc_file + '"')
    #startDayNr = nc.date2num(time_list[0], units=units, calendar=calendar)
    #endDayNr   = nc.date2num(time_list[-1], units=units, calendar=calendar)
    #time       = np.arange(startDayNr, endDayNr+1)
    nc_trg     = nc.Dataset(nc_file, 'w', format=Format, zlib=True)
    #Return_Period       = np.arange(1,RP+1) # in here is defined an array equal to 9 return periods

    print('Setting up dimensions and attributes. lat: ' + str(len(y))+ " lon: " + str(len(x)))
    #nc_trg.createDimension('Return_Period', RP) #NrOfDays*8
    nc_trg.createDimension('lat', len(y))
    nc_trg.createDimension('lon', len(x))
    #DateHour = nc_trg.createVariable('Return_Period','f8',('Return_Period',))
    #DateHour.units = units
    #DateHour.calendar = calendar
    #DateHour.standard_name = 'Return_Period'
    #DateHour.long_name = 'Return_Period'
    #DateHour.axis = 'T'
    #DateHour[:] = Return_Period
    y_var   = nc_trg.createVariable('lat','f4',('lat',))
    y_var.standard_name = 'latitude'
    y_var.long_name = 'latitude'
    y_var.units = 'degrees_north'
    y_var.axis = 'Y'
    x_var = nc_trg.createVariable('lon','f4',('lon',))
    x_var.standard_name = 'longitude'
    x_var.long_name = 'longitude'
    x_var.units = 'degrees_east'
    x_var.axis = 'X'
    y_var[:] = y
    x_var[:] = x
    projection= nc_trg.createVariable('projection','c')
    projection.long_name = 'wgs84'
    projection.EPSG_code = 'EPSG:4326'
    projection.proj4_params = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    projection.grid_mapping_name = 'latitude_longitude'

    # now add all attributes from user-defined metadata
    for attr in metadata:
        nc_trg.setncattr(attr, metadata[attr])
    nc_trg.sync()
    nc_trg.close()

#def append_nc(nc_file, var_name, dtype='f4', chunksizes=(1, 128, 128), fill_value=-9999, metadata={}):
def append_nc(nc_file, var_name, dtype='f4', chunksizes=(128, 128), fill_value=-9999, metadata={}):

    """
    Write a new (empty) variable to target NetCDF file. 
    input:
    ::

        nc_file:         NetCDF object, referring to target file
        var_name:       String, variable name of source NetCDF
        metadata:       dictionary of attributes, belonging to the variable. 

    """
    
    # add the variable
    nc_obj = nc.Dataset(nc_file, 'a')
    #variab = nc_obj.createVariable(var_name, dtype,('Return_Period', 'lat', 'lon',),chunksizes=chunksizes,fill_value=fill_value,zlib=True)
    variab = nc_obj.createVariable(var_name, dtype,('lat', 'lon',),chunksizes=chunksizes,fill_value=fill_value,zlib=True)

    # add some general attributes usually used in lat lon data
    variab.coordinates   = 'lat lon'
    # if a attributes dictionary exists, then append attributes from this dictionary
    if metadata:
        for attribute in metadata:
            variab.setncattr(attribute, metadata[attribute])
    nc_obj.sync()
    nc_obj.close()
