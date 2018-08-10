# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,num2date
import pandas as pd

iso = 'BEN'
print(iso)
os.chdir('/Users/peterpfleiderer/Projects/country_analysis/')
import country_analysis as country_analysis; reload(country_analysis)

# data will be stored in working_directory
COU=country_analysis.country_analysis(iso,working_directory='data/'+iso+'/',additional_tag='')

# # ##############
# # # MASKS
# # ##############
#
# # country mask
# COU.create_mask_country('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','tmp',overwrite=False)
#
# # masks for administrative regions
# COU.create_mask_admin('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','tmp',overwrite=False)

# # population weighted
# #COU.create_mask_country('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','tmp',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
# #COU.create_mask_country('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','tmp',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
#
# # ##############
# # # Zoom - save small netcdf files covering the country
# # ##############
# COU.country_zoom('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.tmp.dat.nc',var_name='tmp',given_var_name='tas',data_type='CRU',time_format='monthly',overwrite=True)
# COU.country_zoom('/Users/peterpfleiderer/Desktop/peter/Documents/data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc',var_name='pre',given_var_name='pr',data_type='CRU',time_format='monthly',overwrite=True)
#
# COU.classify_ensemble(filters=['CRU'])
# COU.save_data()
COU.load_data()
COU.period_statistic([1940,1960],'test')
COU.period_events_above_thresh([1940,1960],'thresh_test',names=[u'pr_CRU_monthly])
asdasd
# merges historical with rcp scenarios
COU.hist_merge()

# creates ensemble mean files
COU.ensemble_mean()

# show loaded data
COU.summary()

# calculates area averages
COU.area_average('lat_weighted',overwrite=True)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)

# compresses the data package into a zip file
COU.zip_it()
