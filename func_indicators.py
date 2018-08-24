#import os
#os.chdir(r'/home/adiaz/WorldBank')

import numpy as np
from xlrd import open_workbook
import netCDF4 as nc
from scipy.interpolate import interp1d
import pdb
from fun_read_rasters import fun_read_rasters




def func_indicators(formatarray,indicator,condition,options_exposure_file,inun_list,options_land_uses_file,options_max_dam_file,cell_size_km2,options_dam_curves_path,sub_indicators,options_geog_path,options_name_var_geogu,options_name_var_gdp,options_name_var_pop,options_name_var_landu):
    """ Passes arguments on to the respective impact indicator function """
    #print('urb damage1')
    #print(indicator)
    if (indicator=='Urban_Damage'):
        print('urb damage')
        if (len(inun_list)==0):
            output=np.zeros(condition.shape)
        else:
            output=func_indicator_Urban_Damage(formatarray,inun_list,condition,options_land_uses_file,options_max_dam_file,cell_size_km2,options_dam_curves_path,options_name_var_landu)
    if (indicator=='GDPexp'):
        if (len(inun_list)==0):
            output=np.zeros(condition.shape)
        else:
            output=func_indicator_GDPexp(formatarray[2],options_exposure_file,condition,options_name_var_gdp)
    if (indicator=='POPexp'):
        if (len(inun_list)==0):
            output=np.zeros(condition.shape)
        else:
            output=func_indicator_POPexp(formatarray[1],options_exposure_file,condition,options_name_var_pop)
    if (indicator=='Building_Damage'):
        if (len(inun_list)==0):
            output=np.zeros(condition.shape)
        else:
            output=func_indicator_Building_Damage(formatarray,inun_list,condition,options_land_uses_file,cell_size_km2,options_dam_curves_path,sub_indicators,options_geog_path,options_name_var_geogu,options_name_var_landu)
    return output;

def func_indicator_GDPexp (formatarray,options_exposure_file,condition,options_name_var_gdp):

    """ Calculation of GDP exposure indicator grid cell """
    #name_file = options_exposure_file.replace('.nc','gdp.h5') #change to nc when available
    name_file = options_exposure_file.replace('.nc','tgdp.nc') #change to nc when available
    # Get population data...
    #gdp_file = nc.Dataset(name_file, 'r', format='h5'); #change to NETCDF4
    #gdp_file = nc.Dataset(name_file, 'r', format='NETCDF4');
    #gdp_file_variable=gdp_file.variables[options_name_var_gdp] # change o the name of the variable
    #gdp_file_variable.set_auto_maskandscale(False)
    #if (gdp_file.variables['lat'][0]<0):
    #    gdp_file_variable=gdp_file_variable[::-1]  # Upside down rasters
    # que se lea desde el archivo bat el .nc
    gdp_file_variable=fun_read_rasters(formatarray,name_file,options_name_var_gdp)
    gdp_file_variable_list = np.extract(condition,gdp_file_variable[:][:])
    #gdp_file.close()
    gdp_file_variable_list[gdp_file_variable_list==-9999]=0
    del(gdp_file_variable)
    # Carry out analysis...
    return gdp_file_variable_list;


def func_indicator_POPexp (formatarray,options_exposure_file,condition,options_name_var_pop):
    """ Calculation of Population exposure indicator for each grid cell """
    #to change the issue of the name
    #name_file = options_exposure_file.replace('.tif','.tif') #change to nc when available
    #name_file = options_exposure_file.replace('.nc','pop.h5') #change to nc when available
    name_file = options_exposure_file.replace('.nc','tpop.nc') #change to nc when available
    # Get population data...
    #pop_file = nc.Dataset(name_file, 'r', format='NETCDF4');
    #pop_file = nc.Dataset(name_file, 'r', format='h5');
    #pop_file_variable=pop_file.variables[options_name_var_pop]
    #pop_file_variable.set_auto_maskandscale(False)
    #if (pop_file.variables['lat'][0]<0):
    #    pop_file_variable=pop_file_variable[::-1]  # Upside down rasters
    print (formatarray,name_file,options_name_var_pop)
    pop_file_variable=fun_read_rasters(formatarray,name_file,options_name_var_pop)
    #print(pop_file_variable)
    #print(np.sum(pop_file_variable))
    pop_file_variable_list = np.extract(condition,pop_file_variable[:][:])
    pop_file_variable_list[pop_file_variable_list==-9999]=0
    #pop_file.close()
    del(pop_file_variable)
    # Carry out analysis...
    return pop_file_variable_list;


def func_indicator_Urban_Damage(formatarray,inun_list,condition,options_land_uses_file,options_max_dam_file,cell_size_km2,options_dam_curves_path,options_name_var_landu):
    """ Calculation of Urban Damage indicator grid cell """
    #print('aquiii 22')
    #Import cell-size data and make list...
    print(condition)
    print(cell_size_km2.shape)
    cell_size_km2_list = np.extract(condition,cell_size_km2[:][:])
    print'uno despues', cell_size_km2_list
    # Get curves data...   land uses will be variable not fix
    book = open_workbook(options_dam_curves_path.replace('curves_2010.xlsx','CURVES_2010.xlsx'), on_demand=True)
    #book = open_workbook(options_dam_curves_path, on_demand=True)
    book_sheet = book.sheet_by_index(0)
    curves = [[book_sheet.cell_value(r,col)
    for col in range(book_sheet.ncols)]
        for r in range(book_sheet.nrows)]
    curves=np.array(curves)
    curves_infinite=np.array([[0],[0]],dtype=np.float64)
    curves_infinite[0][0]=999999999999
    curves_infinite[1][0]=curves[1][len(curves[1][:])-1]
    curves=np.concatenate((curves,curves_infinite),axis=1)
    n_landuses = len(curves)-1
    n_depths = len(curves[0][:])
    # Prepare output file...
    for i in range(1,(n_landuses+1)):
        # Get land use data...
        # que lea los formatos .tif
        #print ('llego a leer los usos del suelo')
        print(options_name_var_landu)
        variable_landuse=fun_read_rasters(formatarray[3],options_land_uses_file,options_name_var_landu)
        #print(np.nanmax(variable_landuse))
        #file_landuse = nc.Dataset(options_land_uses_file, 'r', format='h5');
        #file_landuse = nc.Dataset(options_land_uses_file, 'r', format='NETCDF4');
        #landuse = file_landuse.variables[options_name_var_landu];
        #landuse = file_landuse.variables['Urban Land Use'];
        #landuse.set_auto_maskandscale(False);
        #if (file_landuse.variables['lat'][0]<0):
        #    landuse=landuse[::-1]  # Upside down rasters
        landuse = np.extract(condition,variable_landuse[:][:]);
        #file_landuse.close();
        # Get maximum damage data...
        #file_maxdam = nc.Dataset(options_max_dam_file, 'r', format='h5');
        #maxdam = file_maxdam.variables['data'];
        #maxdam.set_auto_maskandscale(False);
        maxdam=fun_read_rasters(formatarray[4],options_max_dam_file,'data')
        maxdam_list = np.extract(condition,maxdam[:][:]);
        maxdam_list = maxdam_list*1000*1000;   # converts maxdam from USD/m2 to USD/km2...
        #file_maxdam.close()
        # Calculate damage...
        interpolation =  interp1d(curves[0][:],curves[i][:], kind='linear')
        alpha =  interpolation(inun_list);
        damage = np.multiply(landuse, alpha);
        damage = np.multiply(damage, maxdam_list);
        damage = np.multiply(damage,cell_size_km2_list);
        wheredamage=np.where(np.isnan(damage))
        #out_list_Damage =  np.array([0] * len(geogunit_list));    
        out_list_Damage=np.array(damage);
        out_list_Damage[wheredamage]=0;          # This just gives total damage: later may want to make output per land use class?         #  Carry out analysis...
        #del (damage)
        #del (wheredamage)
    return out_list_Damage;


def func_indicator_Building_Damage(formatarray,inun_list,condition,options_land_uses_file,cell_size_km2,options_dam_curves_path,sub_indicators,options_geog_path,options_name_var_geogu,options_name_var_landu):
    """ Calculation of Urban Damage indicator grid cell """
    #print(type(cell_size_km2))
    #Import cell-size data and make list...
    cell_size_km2_list = np.extract(condition,cell_size_km2[:][:])
    # Get curves data...   land uses will be variable not fix
    book = open_workbook(options_dam_curves_path.replace('curves_2010.xlsx','GLOFRIS_GlobalDamageFunctions.xlsx'), on_demand=True)
    #book = open_workbook('C:\Glofris\WRI\GLOFRIS_GlobalDamageFunctions.xlsx', on_demand=True)
    book_sheet = book.sheet_by_index(9)
    curves = [[book_sheet.cell_value(r,col)
    for col in range(book_sheet.ncols)]
        for r in range(book_sheet.nrows)]
    curves=np.array(curves)
    #print(curves)
    # curves_infinite=np.array([[0],[0]],dtype=np.float64)
    # curves_infinite[0][0]=999999999999
    # curves_infinite[1][0]=curves[1][len(curves[1][:])-1]
    # curves=np.concatenate((curves,curves_infinite),axis=1)
    # n_landuses = len(curves)-1
    # n_depths = len(curves[0][:])
    # Prepare output file...

    #geogunit_2=nc.Dataset(options_geog_path+'\geogunit_'+str('2') + '.h5', 'r', format='h5');
    #data=geogunit_2.variables[options_name_var_geogu]
    #data.set_auto_maskandscale(False);

#formatarray.append([options.flood_input_format,options.pop_input_format,options.gdp_input_format,options.landuses_input_format,options.maxdam_input_format,options.geogunits_input_format])

    name=options_geog_path+'\geogunit_'+str('2') +formatarray[5]

    #geogunit_2=nc.Dataset(options_geog_path+'\geogunit_'+str('2') + '.nc', 'r', format='NETCDF4');
    #data=geogunit_2.variables[options_name_var_geogu]
    #data.set_auto_maskandscale(False)

    data=fun_read_rasters(formatarray[5],name,options_name_var_geogu)

    #if (geogunit_2.variables['lat'][0]<0):
    #    data=data[::-1]  # Upside down rasters
    #geogunit_2=nc.Dataset(options_geog_path+'\geogunit_'+str('2') + '.nc', 'r', format='NETCDF4');
    #data=geogunit_2.variables[options.name_var_geogu]
    #data.set_auto_maskandscale(False);
    e=1
    if condition.shape == data.shape:
        print "same dimensions"
    if (data.shape[0]-data.shape[1]) == (condition.shape[1]-condition.shape[0]):
        condition = condition.T
        e=2
        print "invert dimensions", e
    #checking if the dimensions of the array are not inverted
    #print('condition cells',condition.shape,condition[1][1],'geoug cell', data.shape)
    data2 = np.extract(condition,data[:][:]);
    if (e == 2):
        condition = condition.T
    # close the connection
    #print(np.where(data2[:]==157))
    #print(np.nansum(data2)) 
    #in row 999 in the excel file is defined th nan
    data2[np.isnan(data2)]=999;
    #print('the number of geoug cell are:', data.shape)
    #geogunit_2.close();
    #book = open_workbook('C:\Glofris\WRI\GLOFRIS_GlobalMaxDamageValues.xlsx', on_demand=True)
    book = open_workbook(options_dam_curves_path.replace('curves_2010.xlsx','GLOFRIS_GlobalMaxDamageValues.xlsx'), on_demand=True)
    book_sheet = book.sheet_by_index(13)
    maxdam_excel = [[book_sheet.cell_value(r,col)
    for col in range(book_sheet.ncols)]
        for r in range(3,book_sheet.nrows)]
    maxdam_excel = np.array(maxdam_excel)
    #print(maxdam_excel.shape) #testJE
    #maxdam_excel[np.isnan(maxdam_excel)]=0
    # Get land use data...
    landuse=fun_read_rasters('.tif',options_land_uses_file,options_name_var_landu)

    #file_landuse = nc.Dataset(options_land_uses_file, 'r', format='NETCDF4');
    #file_landuse = nc.Dataset(options_land_uses_file, 'r', format='h5');
    #landuse = file_landuse.variables['Urban Land Use'];
    #landuse = file_landuse.variables[options_name_var_landu];
    #landuse = file_landuse.variables['data'];
    #landuse.set_auto_maskandscale(False);
    #if (file_landuse.variables['lat'][0]<0):
    #    landuse=landuse[::-1]  # Upside down rasters
    landuse = np.extract(condition,landuse[:][:]);
    #file_landuse.close();
    out_list=np.array(landuse)
    out_list[:]=0
    for i in range(0,len(sub_indicators)):
        # Get maximum damage data...
        # file_maxdam = nc.Dataset(options_max_dam_file, 'r', format='h5');
        # maxdam = file_maxdam.variables['data'];
        # maxdam.set_auto_maskandscale(False);
        # maxdam_list = np.extract(condition,maxdam[:][:]);
        # maxdam_list = maxdam_list*1000*1000;   # converts maxdam from USD/m2 to USD/km2...
        # file_maxdam.close()
        # Calculate damage...
        position_curve=np.where(curves[:,0]==sub_indicators[i])
        position_curve=np.array(position_curve)       
        interpolation =  interp1d(np.array(curves[0][2:len(curves[0][:])],dtype=np.float64),np.array(curves[position_curve[0][0]][2:len(curves[0][:])],dtype=np.float64), kind='linear')
        alpha =  interpolation(inun_list);
        if (sub_indicators[i]=='Residential'):
            column = 3
        if (sub_indicators[i]=='Commercial'):
            column = 4
        if (sub_indicators[i]=='Industrial'):
            column = 5
        print(maxdam_excel)
        print(maxdam_excel[999][:])

        #print(data.tolist())
        #print()
        maxdam_vector = maxdam_excel[data2.tolist(),column]
        #print(maxdam_vector)
        #print (maxdam_vector[0],maxdam_vector[1],maxdam_vector[2])
        #print(type(landuse))
        #print(maxdam_vector)
        maxdam_vector=np.array(maxdam_vector,dtype=np.float64)
        maxdam_vector = maxdam_vector*1000*1000; #from euros/m to euros/km #why did you comment that out earlier? #questionJE 
        damage = np.multiply(landuse, alpha);
        damage = np.multiply(damage, maxdam_vector);
        damage = np.multiply(damage,cell_size_km2_list);
        wheredamage=np.where(np.isnan(damage))
        #out_list_Damage =  np.array([0] * len(geogunit_list));
        out_list_Damage=np.array(damage);
        del(damage)
        out_list_Damage[wheredamage]=0;          # This just gives total damage: later may want to make output per land use class?         #  Carry out analysis...
        out_list=out_list+out_list_Damage
        #print "**********************"#testJE
        #print(sub_indicators[i]) #testJE
        #print 'np.nansum(out_list_Damage):', np.nansum(out_list_Damage) #testJE
        #print 'np.nansum(out_list):', np.nansum(out_list) #testJE

        #del (damage)
        #del (wheredamage)
    return out_list;

