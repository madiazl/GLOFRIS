import numpy as np
from netcdf_funcs import *
import time
from func_write import func_write
from func_aggregation import func_aggregation

def func_write_nc(options_output,out_list,indicators,x, y,time_list,RP,metadata_global,inun_index,var_w,t):
    """ Prepare indicator impact data to be written to nc file """
    assert isinstance(t, object)
    print('en el loop write nc')
    #print(t)
    #print(indicators)
    #print(var_w)
    #print(np.nanmax(out_list[t]))
    if (indicators=='Building_Damage'):
        #t=0;
        out_list_Damage=np.array(out_list[t]);
        print(out_list_Damage)
        metadata_var = {'units': 'USD','standard_name': 'Risk_Results','long_name': 'Urban Damage $','comment': 'Potential economical loses in the {:s}'} #name invented where to put all the information (money at 2005 year)
        filename = options_output.replace('.nc','_Building_Damage.nc')
        write_nc = func_write(filename, x, y, time_list,RP, metadata_global,metadata_var,inun_index,out_list_Damage)
        #out_list[t]=0
        var_w[t]="Building_Damage_nc"
    if (indicators=='Urban_Damage'):
        #t=0;
        out_list_Damage=np.array(out_list[t]);
        metadata_var = {'units': 'USD','standard_name': 'Risk_Results','long_name': 'Urban Damage $','comment': 'Potential economical loses in the {:s}'} #name invented where to put all the information (money at 2005 year)
        filename = options_output.replace('.nc','_Urban_Damage.nc')
        write_nc = func_write(filename, x, y, time_list,RP, metadata_global,metadata_var,inun_index,out_list_Damage)
        #out_list[t]=0
        var_w[t]="Urban_Damage_nc"
    if (indicators=='GDPexp'):
        #t=1;
        out_list_GDP = np.array(out_list[t]);
        metadata_var = {'units': 'GDP exposed','standard_name': 'Risk_Results','long_name': 'GDP exposed'} #name invented
        filename = options_output.replace('.nc','_GDPexp.nc')
        write_nc = func_write (filename, x, y, time_list,RP, metadata_global,metadata_var,inun_index,out_list_GDP)
        #out_list[t]=0
        var_w[t]="GDPexp_nc"
    if (indicators=='POPexp'):
        #t=2;
        out_list_POP = np.array(out_list[t]);
        print(out_list_POP.shape)
        metadata_var = {'units': 'POP exposed','standard_name': 'Risk_Results','long_name': 'Population exposed'} #name invented
        filename = options_output.replace('.nc','_POPexp.nc')
        write_nc = func_write (filename, x, y, time_list,RP, metadata_global,metadata_var,inun_index,out_list_POP)
        #out_list[t]=0
        var_w[t]="POPexp_nc"  
    return (var_w);

