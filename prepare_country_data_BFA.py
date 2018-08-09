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
COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',COU._working_directory+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=overwrite)


# ##############
# #CMIP5_BC
# ##############
#tas
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	rcp=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='tas',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)

#pr
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/pr/mon_pr_*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	rcp=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='pr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)

#SPI
for months in [3]:
	all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPI/SPI_*'+str(months)+'m_*')
	for in_file in all_files:
		model=in_file.split('/')[-1].split('_')[1]
		rcp=in_file.split('/')[-1].split('_')[2]
		print rcp,model,in_file
		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)

# #SPEI
# for months in [1,3,6,12]:
# 	all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPEI/SPEI_*'+str(months)+'m_*penman*')
# 	for in_file in all_files:
# 		model=in_file.split('/')[-1].split('_')[-2]
# 		rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
# 		print rcp,model,in_file
# 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m_penman',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=overwrite)



COU.hist_merge()

COU.ensemble_mean()

COU.summary()


COU.area_average('lat_weighted',overwrite=overwrite)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)


COU.zip_it()

# except Exception,e:
# 	print str(e)
# 	print 'issues with '+country_iso
