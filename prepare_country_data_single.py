# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')


country_iso='PHL'
additional_tag='_complete'

os.chdir('/p/projects/tumble/carls/shared_folder/country_analysis/')
COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag=additional_tag)


# # ##############
# # # MASKS
# # ##############

# COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',overwrite=False)
# COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
# COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
# COU.create_mask_admin('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4','pr',COU._working_directory+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=False)

# COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',overwrite=False)
# COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
# COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
# COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas',COU._working_directory+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=False)

COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',overwrite=False)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
#COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',COU._working_directory+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=False)


# # ##############
# # #CMIP5_BC
# # ##############

# #tas 
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tas',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

#tas 
all_files=glob.glob('/home/pepflei/CA/data/year_tas/RCP26/year_tas_*nc4')
for in_file in all_files:
	if in_file.split('/')[-1].split('_')[2]!='anom':
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='tas',given_var_name='tas_yr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

#tas 
all_files=glob.glob('/home/pepflei/CA/data/year_tas/RCP26/year_tas_*nc4')
for in_file in all_files:
	if in_file.split('/')[-1].split('_')[2]=='anom':
			model=in_file.split('/')[-1].split('_')[3]
			rcp=in_file.split('/')[-1].split('_')[4]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='tas',given_var_name='tas_yr_anom',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)


# #tasmax
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tasmax/mon_tasmax_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

# #TXx 
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/TXx/mon_TXx_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #pr 
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/pr/mon_pr_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

# #RX5 
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/RX5/mon_RX5_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX5',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #RX1 
# all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/RX1/mon_RX1_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# # #year_CDD 
# # all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/CDD/year_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[4]
# # 	rcp=in_file.split('/')[-1].split('_')[5]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Apr_Jul_cdd
# # all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/CDD/Apr-Jul_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[4]
# # 	rcp=in_file.split('/')[-1].split('_')[5]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Jun_Sep_cdd 
# # all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/CDD/Jun-Sep_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[4]
# # 	rcp=in_file.split('/')[-1].split('_')[5]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #May_Oct_cdd 
# # all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/CDD/May-Oct_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[4]
# # 	rcp=in_file.split('/')[-1].split('_')[5]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #10day_pr
# # all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/10_day_pr/10day_mean_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[4]
# # 	rcp=in_file.split('/')[-1].split('_')[5]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='pr',given_var_name='10day_pr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='10day',overwrite=True)

# # #SPEI
# # for months in [1,3,6,12]:
# # 	all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPEI/SPEI_*'+str(months)+'m_*thronthwaite*')
# # 	for in_file in all_files:
# # 		model=in_file.split('/')[-1].split('_')[-2]
# # 		rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
# # 		print rcp,model,in_file
# # 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# # #SPEI
# # for months in [1,3,6,12]:
# # 	all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPEI/SPEI_*'+str(months)+'m_*penman*')
# # 	for in_file in all_files:
# # 		model=in_file.split('/')[-1].split('_')[-2]
# # 		rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
# # 		print rcp,model,in_file
# # 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m_penman',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)



# # ##############
# # # CORDEX
# # ##############

# #tasmax 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/tasmax/mon_tasmax_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #tas 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/tas/mon_tas_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tas',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #TXx 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/TXx/mon_TXx_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #pr 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #RX5 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/RX5/mon_RX5_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX5',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #RX1 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/RX1/mon_RX1_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# # #year_CDD 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/CDD/year_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Apr_Jul_cdd
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/CDD/Apr-Jul_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Jun_Sep_cdd 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/CDD/Jun-Sep_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #May_Oct_cdd 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/CDD/May-Oct_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #10day_pr
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/10_day_pr/10day_mean_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='pr',given_var_name='10day_pr',data_type='CORDEX',scenario=rcp,model=model,time_format='10day',overwrite=True)

# # #SPEI
# # for months in [1,3,6,12]:
# # 	all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/SPEI/SPEI_*'+str(months)+'m_*')
# # 	for in_file in all_files:
# # 		model=in_file.split('/')[-1].split('_')[-2]
# # 		rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
# # 		print rcp,model,in_file
# # 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)


# # ##############
# # #CORDEX BC
# # ##############

# #tasmax 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/tasmax/mon_tasmax_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #tas 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/tas/mon_tas_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tas',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #TXx 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/TXx/mon_TXx_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #pr 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/pr/mon_pr_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #RX5 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/RX5/mon_RX5_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX5',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #RX1 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/RX1/mon_RX1_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# # #year_CDD 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/year_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Apr_Jul_cdd
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/Apr-Jul_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #Jun_Sep_cdd 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/Jun-Sep_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #May_Oct_cdd 
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/May-Oct_cdd_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# # #10day_pr
# # all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/10_day_pr/10day_mean_*')
# # for in_file in all_files:
# # 	model=in_file.split('/')[-1].split('_')[2]
# # 	rcp=in_file.split('/')[-1].split('_')[3]
# # 	print rcp,model,in_file
# # 	COU.country_zoom(in_file,var_name='pr',given_var_name='10day_pr',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='10day',overwrite=True)

# # # SPEI
# # for months in [1,3,6,12]:
# # 	all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/SPEI/SPEI_*'+str(months)+'m_*')
# # 	for in_file in all_files:
# # 		model=in_file.split('/')[-1].split('_')[-2]
# # 		rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
# # 		print rcp,model,in_file
# # 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)


# ##############
# # EWEMBI
# ##############

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4',var_name='tas',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tasmax_EWEMBI_1979-2014.nc4',var_name='tasmax',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_TXx_EWEMBI_1979-2014.nc4',var_name='tasmax',given_var_name='TXx',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_pr_EWEMBI_1979-2014.nc4',var_name='pr',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX5',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX1_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX1',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/10day_pr/pr_10day_EWEMBI_1979-2013.nc4',var_name='pr',given_var_name='10day_pr',data_type='EWEMBI',time_format='10day',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/year_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/Jun-Sep_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/May-Oct_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/Apr-Jul_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

# for months in [1,3,6,12]:
# 	COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/PET/SPEI_'+str(months)+'m_thornthwaite_EWEMBI_1979-2014.nc4',var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='EWEMBI',time_format='monthly',overwrite=True)



COU.hist_merge()

COU.ensemble_mean()

COU.summary()

# if len(COU.selection(['year_RX5']))<2:
# 	for data in COU.selection(['RX5']): data.year_max('year_RX5')
# 	for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
# 	COU.ensemble_mean()

# COU.summary()	

COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag=additional_tag)
COU.load_data(quiet=False,load_region_polygons=False) 

COU.area_average('lat_weighted',overwrite=False)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)


COU.zip_it()







