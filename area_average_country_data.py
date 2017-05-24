import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

os.chdir('/Users/peterpfleiderer/Documents/')

country_iso='BEN'

COU=country_analysis(country_iso,'Projects/country_analysis/'+country_iso+'/')
COU.load_data()

# for data in COU.selection(['RX5']): data.year_max('year_RX5')
# for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
# COU.ensemble_mean()

COU.summary()

COU.area_average('lat_weighted',overwrite=True,selection=COU.selection(['tas','CORDEX_BC']))

COU.area_average('lat_weighted',overwrite=False)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)



