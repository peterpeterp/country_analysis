import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/')

GHA=country_analysis('GHA','/p/projects/tumble/carls/shared_folder/country_anaylsis/')


GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',0.0,mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',0.0,mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',0.0)


GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',-180.0)
GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',-180.0,mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',-180.0,mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')

##############
# CMIP5
##############

#tas 
all_files=glob.glob('/p/projects/climber3/knaus/Global/Input_data/CMIP5/*/*/mon_tas_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'tas',meta_data=['tas','CMIP5',rcp,model])

#pr 
all_files=glob.glob('/p/projects/climber3/knaus/Global/Input_data/CMIP5/*/*/mon_pr_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'pr',meta_data=['pr','CMIP5',rcp,model])

#spei 
all_files=glob.glob('/p/projects/climber3/knaus/Global/Indices/SPEI/CMIP5/spei_*_12m*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	GHA.country_zoom(in_file,'SPEI',meta_data=['SPEI_12m','CMIP5',rcp,model])

all_files=glob.glob('/p/projects/climber3/knaus/Global/Indices/SPEI/CMIP5/spei_*_1m*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	GHA.country_zoom(in_file,'SPEI',meta_data=['SPEI_1m','CMIP5',rcp,model])

all_files=glob.glob('/p/projects/climber3/knaus/Global/Indices/SPEI/CMIP5/spei_*_3m*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	GHA.country_zoom(in_file,'SPEI',meta_data=['SPEI_3m','CMIP5',rcp,model])

#rx5 
all_files=glob.glob('/p/projects/tumble/carls/shared_folder/rx5/mon_rx5_*')
for in_file in all_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'rx5',meta_data=['rx5','CMIP5',rcp,model])


##############
# NCEP
##############
GHA.country_zoom('/p/projects/climber3/knaus/Global/Input_data/NCEP/tas_ncep_1948-2014.nc','tas',meta_data=['tas','NCEP'])
GHA.country_zoom('/p/projects/climber3/knaus/Global/Input_data/NCEP/pr_ncep_1948-2014.nc','pr',meta_data=['pr','NCEP'])
all_files=glob.glob('/p/projects/climber3/knaus/Global/Indices/SPEI/NCEP/*')
for in_file in all_files:
	months=in_file.split('/')[-1].split('.')[0].split('_')[-1]
	GHA.country_zoom(in_file,'SPEI',meta_data=['SPEI_'+months,'NCEP'])


##############
# CRU
##############
GHA.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',meta_data=['tas','CRU'])
GHA.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.pre.dat.nc','pre',meta_data=['pr','CRU'])
all_files=glob.glob('/home/pepflei/CA/data/SPEI/CRU/*')
for in_file in all_files:
	months=in_file.split('/')[-1].split('.')[0].split('i')[1].replace('0','')
	GHA.country_zoom(in_file,'spei',meta_data=['SPEI_'+months+'m','CRU'])



# output = open(GHA._working_directory+GHA._iso+'/'+GHA._iso+'_masks.pkl', 'wb')
# pickle.dump(GHA._masks, output)	;	output.close()

# output = open(GHA._working_directory+GHA._iso+'/'+GHA._iso+'_data.pkl', 'wb')
# pickle.dump(GHA._data, output)	;	output.close()

# output = open(GHA._working_directory+GHA._iso+'/'+GHA._iso+'_meta.pkl', 'wb')
# pickle.dump(GHA._meta, output)	;	output.close()





