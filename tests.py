import sys,glob,os
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/country_analysis/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

GHA_SPEI=country_analysis('GHA','SPEI_1m','Scripts/country_analysis/')


GHA_SPEI.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',0.0)


in_files=glob.glob('data/raw/mon_rx5/CMIP5/*/mon_rx5_*_1950-2099.nc4')
for in_file in in_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA_SPEI.country_zoom(in_file,'mon_rx5',meta_data=['CMIP5',rcp,model])

GHA_SPEI.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['OBS','NCEP'])


