
import numpy as np

def func_aggregation (out_list,geogunit_list,geogunit_idnums):
    #agregation dependig on th condition of the agregation
    impact_list_index=np.where(out_list>0); # hacerlo como condicion revisar los tamanyos de los arreglos
    impact_list_index=np.array(impact_list_index);
    impact_list=out_list[impact_list_index];
    impact_list_geogunit = geogunit_list[impact_list_index];
    out=np.zeros((len(geogunit_idnums),2));
    out[:,0]=np.array(geogunit_idnums.astype(int));
    for m in range (0,len(geogunit_idnums)):
        auxiliar=np.where(impact_list_geogunit==geogunit_idnums[m]);
        out[m,1]=np.nansum(impact_list[auxiliar]);
    # Carry out analysis...   
    return out;
