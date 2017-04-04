import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis_obj'] 
except:pass
from country_analysis_obj import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')

SEN=country_analysis('SEN','/p/projects/tumble/carls/shared_folder/country_analysis/')


SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
SEN.create_mask_admin('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/SEN_adm_shp/SEN_adm1')


SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
SEN.create_mask_country('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
SEN.create_mask_admin('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/SEN_adm_shp/SEN_adm1')


SEN.create_mask_country('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/CCLM4_HADGEM2/Hist/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-19501230.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',lat_name='rlat',lon_name='rlon')
SEN.create_mask_country('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/CCLM4_HADGEM2/Hist/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-19501230.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',lat_name='rlat',lon_name='rlon')
SEN.create_mask_country('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/CCLM4_HADGEM2/Hist/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-19501230.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',lat_name='rlat',lon_name='rlon')
SEN.create_mask_admin('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/CCLM4_HADGEM2/Hist/tas_AFR-44_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM4-8-17_v1_day_19491201-19501230.nc','tas','/home/pepflei/CA/masks/shapefiles/SEN_adm_shp/SEN_adm1',lat_name='rlat',lon_name='rlon')

##############
# CORDEX
##############

#pr 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/monthly/*/*/mon_pr_*')
for in_file in all_files:
	rcp=in_file.split('/')[-2]
	model=in_file.split('/')[-3]
	print rcp,model,in_file
	SEN.country_zoom(in_file,'pr',tags={'type':'CORDEX','scenario':rcp,'model':model,'var':'pr'},overwrite=False)

#tas 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/Org_Data/RCM/monthly/*/*/mon_tas_*')
for in_file in all_files:
	rcp=in_file.split('/')[-2]
	model=in_file.split('/')[-3]
	print rcp,model,in_file
	SEN.country_zoom(in_file,'tas',tags={'type':'CORDEX','scenario':rcp,'model':model,'var':'tas'},overwrite=False)

# ##############
# # CORDEX BC
# ##############

# #tas 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/*/*/mon_tas_*')
# for in_file in all_files:
# 	rcp=in_file.split('/')[-2]
# 	model=in_file.split('/')[-3]
# 	print rcp,model,in_file
# 	SEN.country_zoom(in_file,'tas',tags={'type':'CORDEX_BC','scenario':rcp,'model':model,'var':'tas'},overwrite=False)

##############
# CMIP5
##############

#tas 
all_files=glob.glob('/p/projects/climber3/knaus/Global/Input_data/CMIP5/*/*/mon_tas_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	SEN.country_zoom(in_file,'tas',tags={'type':'CMIP5','scenario':rcp,'model':model,'var':'tas'})

#pr 
all_files=glob.glob('/p/projects/climber3/knaus/Global/Input_data/CMIP5/*/*/mon_pr_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	SEN.country_zoom(in_file,'pr',tags={'type':'CMIP5','scenario':rcp,'model':model,'var':'pr'})

# #pr new
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/*/*/mon_pr_*')
# for in_file in all_files:
# 	rcp=in_file.split('_')[-2]
# 	model=in_file.split('_')[-3]
# 	SEN.country_zoom(in_file,'pr',meta_data=['pr','CMIP5',rcp,model])

#rx5 
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/rx5/mon_rx5_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	SEN.country_zoom(in_file,'rx5',tags={'type':'CMIP5','scenario':rcp,'model':model,'var':'rx5'})

# #rx1 new
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/*/*/mon_rx1_*')
# for in_file in all_files:
# 	rcp=in_file.split('_')[-2]
# 	model=in_file.split('_')[-3]
# 	SEN.country_zoom(in_file,'rx1',meta_data=['rx1','CMIP5',rcp,model])

##############
# NCEP
##############
SEN.country_zoom('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas',tags={'type':'NCEP','var':'tas'})
SEN.country_zoom('/p/projects/climber3/knaus/Global/Input_data/NCEP/pr_ncep_1948-2014.nc','pr',tags={'type':'NCEP','var':'pr'})
all_files=glob.glob('/p/projects/climber3/knaus/Global/Indices/SPEI/NCEP/*')
for in_file in all_files:
	months=in_file.split('/')[-1].split('.')[0].split('_')[-1]
	SEN.country_zoom(in_file,'SPEI',tags={'type':'NCEP','var':'SPEI'})


# ##############
# # CRU
# ##############
# SEN.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',tags={'type':'CRU','var':'tas'})
# SEN.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.pre.dat.nc','pre',tags={'type':'CRU','var':'pr'})
# all_files=glob.glob('/home/pepflei/CA/data/SPEI/CRU/*')
# for in_file in all_files:
# 	months=in_file.split('/')[-1].split('.')[0].split('i')[1].replace('0','')
# 	SEN.country_zoom(in_file,'spei',tags={'type':'CRU','var':'SPEI'})


# SEN.hist_merge()

# SEN.average('lat_weighted')
# SEN.average('pop2015_weighted')
# SEN.average('pop1990_weighted')

SEN.hist_merge()

SEN.display()

adas

SEN.ensemble_mean()

SEN.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]})

SEN.prepare_for_download()




