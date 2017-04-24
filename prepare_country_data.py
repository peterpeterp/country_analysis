import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis_obj'] 
except:pass
from country_analysis_obj import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')


country_iso='BEN'
COU=country_analysis(country_iso,'/p/projects/tumble/carls/shared_folder/country_analysis/'+country_iso+'/')


# ##############
# # CORDEX
# ##############

COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
COU.create_mask_admin('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/'+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=True)


##############
# CORDEX BC
##############

#tas 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/mon_tas_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='tas',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


#tasmax 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/mon_tasmax_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


#pr 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/mon_pr_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='pr',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


#mon_RX5 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/mon_RX5_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='pr',given_var_name='mon_RX5',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


#year_RX5 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/year_RX5_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='pr',given_var_name='year_RX5',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


#RX1 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/RX1_*')
for in_file in all_files:
	rcp=in_file.split('_')[-3]
	model=in_file.split('_')[-4]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)

#year_CDD 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/year_CDD_*')
for in_file in all_files:
	rcp=in_file.split('/')[-1].split('_')[4]
	model=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_CDD',data_type='CORDEX_BC',scenario=rcp,model=model,overwrite=True)


##############
# EWEMBI
##############
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/mon_tas_EWEMBI_1979-2014.nc4',var_name='tas',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/mon_TXx_EWEMBI_1979-2014.nc4',var_name='tasmax',given_var_name='TXx',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/mon_pr_EWEMBI_1979-2014.nc4',var_name='pr',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/year_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='year_RX5',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/mon_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='mon_RX5',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/mon_RX1_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX1',data_type='EWEMBI',overwrite=True)
COU.country_zoom('/p/projects/tumble/carls/shared_folder/monthly_mean_and_extremes/EWEMBI/year_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_CDD',data_type='EWEMBI',overwrite=True)


COU.hist_merge()

COU.ensemble_mean()

COU.data_summary()

COU.zip_it()


