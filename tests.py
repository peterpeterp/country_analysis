import sys,glob,os
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

GHA=country_analysis('GHA','SPEI_1m','Scripts/country_analysis/')

GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',-180.0,mask_style='pop1990_weigthed',pop_mask_file='masks/population/population_1990_incrLat.nc')

GHA.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',0.0)

GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',-180.0)


in_files=glob.glob('data/raw/mon_rx5/CMIP5/*/mon_rx5_*_1950-2099.nc4')
for in_file in in_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'mon_rx5',meta_data=['rx5','CMIP5',rcp,model])

GHA.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['SPEI_1m','OBS','NCEP'])






tmp=GHA._data.copy()

meta=[]
for key1 in tmp.keys():
	for key2 in tmp[key1].keys():
		if 'data' in tmp[key1][key2].keys():
			meta.append([key1,key2])
		else:
			for key3 in tmp[key1][key2].keys():
				if 'data' in tmp[key1][key2][key3].keys():
					meta.append([key1,key2,key3])
				else:
					for key4 in tmp[key1][key2][key3].keys():
						if 'data' in tmp[key1][key2][key3][key4].keys():
							meta.append([key1,key2,key3,key4])
