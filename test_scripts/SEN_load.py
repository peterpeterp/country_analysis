import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/allgemeine_scripte/')
from plot_functions import *
sys.path.append('/Users/peterpfleiderer/Documents/')

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')






SEN=country_analysis('SEN','Projects/country_analysis/')
SEN.load_from_tar('Projects/country_analysis/SEN.tar.gz')




# plot
fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(4,4))
for name in SEN._regions:
    SEN.plot_transient(meta_data=['SPEI_1m','NCEP'],mask_style='lat_weighted',region=name.encode('utf-8'),ax=axes,running_mean=240,show=False,ylabel='SPEI',label='asdas')

plt.show()

