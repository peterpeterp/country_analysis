import sys,glob,os
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/')

GHA=country_analysis('GHA','SPEI_1m','/p/projects/tumble/carls/shared_folder/country_anaylsis/')


GHA.create_mask('/p/projects/climber3/knaus/Global/Input_data/CMIP5/RCP2.6/hadgem2-es/mon_tas_hadgem2-es_rcp2.6_1950-2099.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',0.0)

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



GHA.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['SPEI_1m','OBS','NCEP'])


