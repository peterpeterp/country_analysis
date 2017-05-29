import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pylab as plt

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

os.chdir('/Users/peterpfleiderer/Documents/Projects/')

BEN=country_analysis('BEN','country_analysis/BEN/',seasons={'year':range(1,13),'Apr-Jul':[4,5,6,7]})
BEN.load_data()

#BEN.get_warming_slices(GMT_path='wlcalculator/data/cmip5_ver002/',model_real_names={'IPSL':'ipsl-cm5a-lr','HADGEM2':'hadgem2-es','ECEARTH':'hadgem2-es','MPIESM':'mpi-esm-lr'})

# plot a transient
pr=BEN.selection(['pr'])
BEN.area_average(selection=pr,mask_style='lat_weighted')
pr[7].plot_transients(running_mean_years=20)