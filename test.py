# -*- coding: utf-8 -*-
import sys,glob,os,cPickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pylab as plt

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

os.chdir('/Users/peterpfleiderer/Documents/Projects/wlcalculator/app/')
import wacalc.CmipData as CmipData; reload(CmipData)
import wacalc.hadcrut_warming as hadcrut_warming; reload(hadcrut_warming)
os.chdir('/Users/peterpfleiderer/Documents/Projects/wlcalculator/app/')

os.chdir('/Users/peterpfleiderer/Documents/Projects/')

COU=country_analysis('SEN','country_analysis/data/SEN/',seasons={'year':range(1,13),'Apr-Jul':[4,5,6,7]})
#COU.load_data(filename_filter='tas')
COU.load_data(quiet=True,load_mask=True,load_raw=True,load_area_averages=False,load_region_polygons=True)

COU.get_warming_slices(wlcalculator_path='/Users/peterpfleiderer/Documents/Projects/wlcalculator/app/',model_real_names={'IPSL':'ipsl-cm5a-lr','HADGEM2':'hadgem2-es','ECEARTH':'ec-earth','MPIESM':'mpi-esm-lr'})


# ens_mean=COU.selection(['tas','CORDEX_BC','ensemble_mean'])[0]
# COU.period_statistics(periods=COU._warming_slices,selection=COU.selection(['tas','CORDEX_BC']),ref_name='ref')
# COU.period_model_agreement(ref_name='ref')
# ens_mean.display_map(period='diff_2.5-ref')

# plot a transient
# tas=COU.selection(['tas'])
# COU.unit_conversions()

# COU.merge_adm_regions(['Mono','Kouffo'])
# COU.create_mask_admin(tas[0].raw_file,'tas',regions=['Kouffo+Mono'])
# COU.area_average('lat_weighted',overwrite=False,regions=['Alibori'],selection=tas)
# COU.area_average('lat_weighted',overwrite=False,regions=['Kouffo+Mono'],selection=tas)

# tas[0].plot_transients()


#tas[0].plot_map(to_plot='empty',show_region_names=True,color_bar=False)

#COU.period_statistics(periods={'badf':[1980,2000],'dasdas':[2000,2020]})

# COU.get_region_area('Atakora')

# COU.get_warming_slices(warming_lvls=[1.5,2,2.5,3],ref_period=[1986,2006],model_real_names={'IPSL':'ipsl-cm5a-lr','HADGEM2':'hadgem2-es','ECEARTH':'ec-earth','MPIESM':'mpi-esm-lr'},wlcalculator_path='/Users/peterpfleiderer/Documents/Projects/wlcalculator/app/')

# print COU._warming_slices




