# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')
try:del sys.modules['country_analysis']
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')


country_iso='BFA'
overwrite=False

print country_iso
COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag='')

os.chdir('/p/projects/tumble/carls/shared_folder/country_analysis/')
COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag='')


# ##############
# # MASKS
# ##############

COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',overwrite=overwrite)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=overwrite)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=overwrite)
COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',overwrite=overwrite)


# ##############
# #CMIP5_BC
# ##############
#tas
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_*rcp45*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	COU.country_zoom(in_file,name='tas_'+model+'_rcp45',var_name='tas',overwrite=overwrite)

#tas
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_*historical*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	COU.country_zoom(in_file,name='tas_'+model+'_rcp45',var_name='tas',overwrite=overwrite)
#
# #pr
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/pr/mon_pr_*rcp45*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)
#
# #SPI
# for months in [3]:
# 	all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPI/SPI_*'+str(months)+'m*')
# 	for in_file in all_files:
# 		model=in_file.split('/')[-1].split('_')[1]
# 		rcp=in_file.split('/')[-1].split('_')[2]
# 		print rcp,model,in_file
# 		COU.country_zoom(in_file,var_name='spi',given_var_name='SPI_'+str(months)+'m',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)
#

COU.save_data()
