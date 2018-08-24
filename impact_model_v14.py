import os
os.chdir(r'/home/madiazl/IWC')
#loading all the packages required to run the model
from func_aggregation import func_aggregation
from func_indicators import func_indicators
from func_write_nc import func_write_nc
from func_write_aggregation import *
from func_check_newnc import func_check_newnc
from fun_read_rasters import fun_read_rasters
from scipy.io import netcdf
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import scipy.io
import itertools
from scipy.interpolate import interp1d
from scipy import interpolate
from xlrd import open_workbook
from netcdf_funcs import *
import numpy as np
import time
import sys
import shutil
from optparse import OptionParser
import datetime
import numpy as np
#import pyproj
import impact_model_v12_lib as cl
#import netcdf_funcs as nc_func
import netCDF4 as nc
import pdb

#print('here2')

# here start the main function of the GLOFRIS model

def main():
    #print('here2')
    ### Read input arguments: from all defined flags in the configuration and batch file
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-q', '--quiet',
                      dest='verbose', default=True, action='store_false',
                      help='do not print status messages to stdout')
    parser.add_option('-e', '--exposure',
                      nargs=1, dest='exposure_file',
                      default='exposure_file',
                      help='exposure file from by PBL')
    parser.add_option('-f', '--flood',
                      nargs=1, dest='flood_file',
                      default='flood_file',
                      help='flood file from Deltares')
    parser.add_option('-l', '--land_uses',
                      nargs=1, dest='land_uses_file',
                      default='land_uses_file',
                      help='land uses file from PBL')
    parser.add_option('-c', '--con', dest='inifile',
                      default='impact_model_configuration.ini', nargs=1,
                      help='ini configuration file')
    parser.add_option('-o', '--output',
                      dest='output', default='output_file',
                      help='output file')
    parser.add_option('-g', '--subsidence',
                      dest='subsidence_map', default='',
                      help='Subsidence map (NetCDF)')
    parser.add_option('-t', '--time',
                      dest='time', default='2010-01-01 00:00:00',
                      help='Time stamp of flood condition')
    parser.add_option('-r', '--return_period',
                      dest='RP', default='2,5,10,25,50,100,250,500,1000',
                      help='Time stamp of flood condition')
    parser.add_option('-m', '--maximum_damage',
                      dest='max_dam_file', default='maximum_damage',
                      help='is the maximum damage per cell in indicator Urban Damage is providen by maps. Indicator Building Damage work with a table of values per country ')
    parser.add_option('-v', '--name_var_flood',
                      nargs=1, dest='name_var_flood',
                      default='name_var_flood',
                      help='flood file from Deltares')


    # pass arguments from batch file to options
    (options, args) = parser.parse_args()
    #print(options)
    # required input
    if not options.output:   # if destination is not given
        parser.error('destination output not given') #Does that line make any sense? Because in options.output will always be at least output_file? And what is parser.error, where do we pick it up? #questionJE
    if not os.path.exists(options.inifile):
        print 'path to ini file cannot be found'
        sys.exit(1)

    # file names and directory bookkeeping
    options.output = os.path.abspath(options.output)
    options.dest_path = os.path.split(options.output)[0]
    datelog = time.ctime()
    datelog=datelog.replace(':','_')
    logfilename = os.path.join(options.dest_path, 'GLOFRIS_logfile_'+str(datelog.replace(' ','_'))+'.log')
    #print('here3')

    # create dir if not exist
    if not os.path.isdir(options.dest_path): #not sure if we need this, isn't the directory the same on as the one we specify os.chdir at the begining of the script? #questionJE
        os.makedirs(options.dest_path)
    # delete old destination and log files
    else:
        if os.path.isfile(options.output):
            os.unlink(options.output)
        if os.path.isfile(logfilename):
            os.unlink(logfilename)


    # set up the logger
    logger, ch = cl.setlogger(logfilename, 'IMPACT_MODEL', options.verbose)

    # write to log file
    logger.info('$Id: $') #What is this? #questionJE
    logger.info('Exposure map: {:s}'.format(options.exposure_file))
    logger.info('Flood map: {:s}'.format(options.flood_file))

    logger.info('Land Uses map: {:s}'.format(options.land_uses_file))

    #logger.info('Impact Indicator: {:s}'.format(options.impact_indicator)) #I would put this entire log file block after the config-file has been read so we can include indicators, geounits etc. and have the begining of the log file writing in one block and not in pieces
    logger.info('Subsidence map: {:s}'.format(options.subsidence_map))
    logger.info('Output file: {:s}'.format(options.output))
    logger.info('Time of flood conditions: {:s}'.format(options.time))

    ### Read config file
    # open config-file
    config = cl.open_conf(options.inifile)

    # read settings
    options.cell_size_file = cl.configget(config,'maps','cell_size_file',True)
    if (options.land_uses_file=='land_uses_file'):
        options.land_uses_file = cl.configget(config,'maps','land_uses_file',True)    #land uses variable with the flag -l
    #print(options.max_dam_file)
    if (options.max_dam_file=='maximum_damage'):
        options.max_dam_file = cl.configget(config,'maps','max_dam_file',True)
    options.geogunit = cl.configget(config,'conf','geogunit',True)
    options.geog_path = cl.configget(config,'conf','geog_path',True)
    options.name_var_geogu = cl.configget(config,'conf','name_var_geogu',True)
    options.name_var_landu = cl.configget(config,'conf','name_var_landu',True)
    if (options.name_var_flood=='name_var_flood'):
        options.name_var_flood = cl.configget(config,'conf','name_var_flood',True)
    options.name_var_pop = cl.configget(config,'conf','name_var_pop',True)
    options.name_var_gdp = cl.configget(config,'conf','name_var_gdp',True)
    options.dam_curves_path = cl.configget(config,'conf','dam_curves_path',True)
    options.flood_units = cl.configget(config,'conf','flood_units',True)

    options.flood_input_format = cl.configget(config,'formats','flood_input_format',True)
    options.pop_input_format = cl.configget(config,'formats','pop_input_format',True)
    options.gdp_input_format = cl.configget(config,'formats','gdp_input_format',True)
    options.landuses_input_format = cl.configget(config,'formats','landuses_input_format',True)
    options.maxdam_input_format = cl.configget(config,'formats','maxdam_input_format',True)
    options.geogunits_input_format = cl.configget(config,'formats','geogunits_input_format',True)

    #print(options.geogunits_input_format)
    #print(type(options.geogunits_input_format))
    formatarray=[]
    formatarray.append([options.flood_input_format,options.pop_input_format,options.gdp_input_format,options.landuses_input_format,options.maxdam_input_format,options.geogunits_input_format])
    formatarray=formatarray[0]

    options.ind = cl.configget(config,'indicators','ind',True)
    options.sub_ind = cl.configget(config,'indicators','sub_ind',True)
    options.aggregation = cl.configget(config,'aggregation','aggregation_by_geogunits',True)

    # terminate due to missing data
    if not os.path.exists(options.cell_size_file):
        logger.error('path to cell size file {:s} cannot be found'.format(options.cell_size_file))
        sys.exit(1)
    if not os.path.exists(options.land_uses_file):
        logger.error('path to land uses file {:s} cannot be found'.format(options.land_uses_file))
        sys.exit(1)
    if not os.path.exists(options.max_dam_file):
        logger.error('path to max. damage file {:s} cannot be found'.format(options.max_dam_file))
        sys.exit(1)
    #print(options.dam_curves_path)
    if not os.path.exists(options.dam_curves_path):
        logger.error('path to damage curves {:s} cannot be found'.format(options.dam_curves_path))
        sys.exit(1)

    # write to log file
    logger.info('Cell Size file: {:s}'.format(options.cell_size_file))
    logger.info('Land Uses file: {:s}'.format(options.land_uses_file))
    logger.info('Max. Damage file: {:s}'.format(options.max_dam_file))
    logger.info('Damage Curves file: {:s}'.format(options.dam_curves_path))
    logger.info('Flood hazard map units: {:s}'.format(options.flood_units))

    # create metadata
    metadata_global = {}
    # add metadata from the section [metadata]
    meta_keys = config.options('metadata')
    for key in meta_keys:
        metadata_global[key] = config.get('metadata', key)
    # add a number of metadata variables that are mandatory
    metadata_global['config_file'] = os.path.abspath(options.inifile)
    metadata_global['history']=time.ctime()

    # split up user choises from the config input and assign them into lists
    geougunit = map(int,options.geogunit.split(','))
    indicators = map(str,options.ind.split(','))
    sub_indicators = map(str,options.sub_ind.split(','))
    aggregation = map(str,options.aggregation.split(','))

    flood_units = map(str,options.flood_units.split(','))
    flood_input_format = map(str,options.flood_input_format.split(','))
    pop_input_format = map(str,options.pop_input_format.split(','))
    gdp_input_format = map(str,options.gdp_input_format.split(','))
    landuses_input_format = map(str,options.landuses_input_format.split(','))
    maxdam_input_format =  map(str,options.maxdam_input_format.split(','))
    geogunits_input_format = map(str,options.geogunits_input_format.split(','))
    #print (geogunits_input_format)
    #print(type(geogunits_input_format))

    RP = map(int,options.RP.split(','))


    ### GLOFRIS Impact model
    time1 = time.time()
    #print(options.flood_file)

    # create a function call open file according to if nc do if h5 do

    #inun_file = nc.Dataset(options.flood_file, 'r', format='NETCDF4');
    #variable=inun_file.variables[options.name_var_flood]
    #variable.set_auto_maskandscale(False)
    #if (len(variable.shape)==3):                                     #3 dimensional arrays coming from edwin
    #    variable=variable[0,:,:]
    #if (inun_file.variables['lat'][0]<0):
    #    variable=variable[::-1]  # Upside down rasters
    print('flood map')
    print(options.flood_input_format,options.name_var_flood)
    variable=fun_read_rasters(options.flood_input_format,options.flood_file,options.name_var_flood)
    print "pass 12"
    variable[variable[:][:]>20000]=0
    inun_index = np.where(variable[:][:]>0)
    condition = variable[:][:]>0
    inun_list = np.extract(condition,variable[:][:])
    inun_list = np.array(inun_list,dtype=np.float64)

    # if units cm then inun_list/100
    if (options.flood_units=='cm'):
        print(options.flood_units)
        print('divided by 100')
        inun_list = inun_list/100     #constantmaskedvalue ##  is unmasked the raster the values are multiply by 10
    # here we should put a if in which acording the geounit will agregate or not also should be work for geoug 10

    #if not geogunit = 10:
        #for geogunits
        #   do aggretation
    file_cell_size_km2 = nc.Dataset(options.cell_size_file, 'r', format='h5');
    cell_size_km2=file_cell_size_km2.variables['data']
    cell_size_km2.set_auto_maskandscale(False)
    # prepare a nice x and y axis
    x = np.linspace(-180+1./240, 180-1./240, 43200)
    y = np.linspace(-90+1./240, 90-1./240, 21600)
    time_list = [datetime.datetime(2016, 1, 1, 0, 0)]

    # create empty lists for store indicator calculations output. the output is a vector with the flooded cells results of the indicator analisis.
    out_list=list(range(0,len(indicators),1)) # list of output from indicator calculation
    var_w=[0] * len(indicators)  # created to later check for which indicators nc-files have been produced
    # loop for impact indicator calculation per cell
    for t in range(0,len(indicators)):
        #print('length')
        #print(len(indicators))
        metadata_list=list(range(0,t+1,1))
        # impact indicator
        output = func_indicators(formatarray,indicators[t],condition,options.exposure_file,inun_list,options.land_uses_file,options.max_dam_file,cell_size_km2,options.dam_curves_path,sub_indicators,options.geog_path,options.name_var_geogu,options.name_var_gdp,options.name_var_pop,options.name_var_landu)
        logger.info('Calculated indicator: {:s}' .format(indicators[t]))
        print "########################"#testJE
        print "indicator: " + indicators[t]
        #print 'np.nansum(output):', np.nansum(output)  # why do we need this? (what were the results for this in the matlab version?) Is this something we want to keep here in the end? #questionJE
        out_list[t]=output
        write_nc = func_write_nc(options.output,out_list,indicators[t],x, y,time_list,RP,metadata_global,inun_index,var_w,t)
        # check if a nc file for the indicator have been created
#        created_nc = func_check_newfile(indicators[t], options.output)
#        logger.info(created_nc)

    logger.info('Finished the calculation of impact indicators')
    file_cell_size_km2.close()

    # aggregate impacts per geounit
    if (aggregation[0]=='yes'):
        logger.info('Start geounit aggregation')
        # loop for all geogunits
        for i_geogunits in range(0,len(geougunit)):
            # prepare variables for current geougunit
            #geogunits=nc.Dataset(options.geog_path+'\geogunit_'+str(geougunit[i_geogunits]) + '.h5', 'r', format='h5');
            #data=geogunits.variables[options.name_var_geogu]
            #print(options.geog_path+'\geogunit_'+str(geougunit[i_geogunits]) + '.nc')
            geogunits=nc.Dataset(options.geog_path+'/geogunit_'+str(geougunit[i_geogunits]) + '.nc', 'r', format='NETCDF4');
            data=geogunits.variables[options.name_var_geogu]
            data.set_auto_maskandscale(False)
            #print('max geog')
            #print(np.nanmax(data))
            if (geogunits.variables['lat'][0]<0):
                data=data[::-1]  # Upside down rasters
            e=1
            if variable.shape == data.shape:
                print "same dimensions"
            if (data.shape[0]-data.shape[1]) == (variable.shape[1]-variable.shape[0]):
                condition = condition.T
                e=2
                print "invert dimensions", e
            geogunit_list = np.extract(condition,data[:][:])
            #checking if the dimensions of the array are not inverted
            if (e == 2):
                condition = condition.T
            # close the connection
            geogunits.close()
            # read the geougunit from the excel list
            book = open_workbook(options.geog_path+'/geogunit_'+str(geougunit[i_geogunits])+ '_list.xlsx', on_demand=True)
            book_sheet = book.sheet_by_index(0)
            geogunitxls = [[book_sheet.cell_value(r,col)
            for col in range(book_sheet.ncols)]
                for r in range(book_sheet.nrows)]
            geogunitxls= np.array(geogunitxls)
            geogunit_idnums=np.array(geogunitxls[:,0],dtype=np.float64)
            geogunit_idnames=np.array(geogunitxls[:,1],dtype=np.str)
            logger.info('Prepared geougunit {0} for aggregation' .format(geougunit[i_geogunits]))

            #print('Preparing some random data') #Does this have any meaning or is this an old comment to a test you previously did? #questionJE
            # aggregate
            write_aggregation = func_write_aggregation(aggregation,options.output,geougunit,i_geogunits,geogunit_list,geogunit_idnums,out_list,indicators);
            #print(write_aggregation)
            logger.info('Aggregated geougunit: {0}' .format(geougunit[i_geogunits]))

    #free memory space
    #inun_file.close()
    del(geogunits)
    #del(inun_file)



    # close logger file
    logger, ch = cl.closeLogger(logger, ch)
    del logger, ch
    sys.exit(0)


    # LEGENDS of input variables
    # Geographical units (geounits)
    # 1: FPU
    # 2: Countries: GADM2
    # 3: World regions: file of Matti Kummu
    # 4: OLD: Basins file of WDI (new = geogunit 7)
    # 5: Nigeria, States (from World Bank)
    # 6: River basins: from Hessel...
    # 7: WRI Basins
    # 8: WRI Administrative Units
    # 9: ECA shapefile f

if __name__ == "__main__":
    main()

