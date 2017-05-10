import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')


country_iso='SEN'
COU=country_analysis(country_iso,'/p/projects/tumble/carls/shared_folder/country_analysis/'+country_iso+'/',additional_tag='_cdd')

# ##############
# # CORDEX
# ##############

COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
COU.create_mask_country('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
COU.create_mask_admin('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_ECEARTH_historical_19510101-20101231.nc4','pr','/home/pepflei/CA/masks/shapefiles/'+country_iso+'_adm_shp/'+country_iso+'_adm1')

COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc')
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc')
COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas','/home/pepflei/CA/masks/shapefiles/'+country_iso+'_adm_shp/'+country_iso+'_adm1')

# #tas 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_tas_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tas',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #tasmax 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_tasmax_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #pr 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/mon_pr_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',data_type='CORDEX',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #year_RX5 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/year_RX5_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='year_RX5',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# #RX1 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/RX1_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[1]
# 	rcp=in_file.split('/')[-1].split('_')[2]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# #year_CDD 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/year_CDD_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_CDD',data_type='CORDEX',scenario=rcp,model=model,time_format='yearly',overwrite=True)



##############
# CORDEX BC
##############

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
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/RX1/RX1_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[1]
# 	rcp=in_file.split('/')[-1].split('_')[2]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)

# #year_CDD 
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/year_cdd_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

#Apr_Jul_cdd
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/Apr-Jul_cdd_*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	rcp=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

#Jun_Sep_cdd 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/Jun-Sep_cdd_*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	rcp=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

#May_Oct_cdd 
all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/CDD/May-Oct_cdd_*')
for in_file in all_files:
	model=in_file.split('/')[-1].split('_')[2]
	rcp=in_file.split('/')[-1].split('_')[3]
	print rcp,model,in_file
	COU.country_zoom(in_file,var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='yearly',overwrite=True)

# #10day_pr
# all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/10_day_pr/10day_mean_*')
# for in_file in all_files:
# 	model=in_file.split('/')[-1].split('_')[2]
# 	rcp=in_file.split('/')[-1].split('_')[3]
# 	print rcp,model,in_file
# 	COU.country_zoom(in_file,var_name='pr',given_var_name='10day_pr',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='10day',overwrite=True)

# SPEI
# for months in [1,3,6,12]:
# 	all_files=glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/PET/SPEI_thornthwaite_*'+str(months)+'m_*')
# 	for in_file in all_files:
# 		model=in_file.split('/')[-1].split('_')[-4]
# 		rcp=in_file.split('/')[-1].split('_')[-3]
# 		print rcp,model,in_file
# 		COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='CORDEX_BC',scenario=rcp,model=model,time_format='monthly',overwrite=True)


##############
# EWEMBI
##############
# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4',var_name='tas',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_TXx_EWEMBI_1979-2014.nc4',var_name='tasmax',given_var_name='TXx',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_pr_EWEMBI_1979-2014.nc4',var_name='pr',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX5',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX1_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX1',data_type='EWEMBI',time_format='monthly',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/10day_pr/pr_10day_EWEMBI_1979-2013.nc4',var_name='pr',given_var_name='10day_pr',data_type='EWEMBI',time_format='10day',overwrite=True)

# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/year_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='year_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/Jun-Sep_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='Jun_Sep_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/May-Oct_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='May_Oct_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/Apr-Jul_CDD_EWEMBI_1979-2014.nc4',var_name='consecutive_dry_days_index_per_time_period',given_var_name='Apr_Jul_cdd',data_type='EWEMBI',time_format='yearly',overwrite=True)

# for months in [1,3,6,12]:
# 	COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/PET/SPEI_'+str(months)+'m_thornthwaite_EWEMBI_1979-2014.nc4',var_name='spei',given_var_name='SPEI_'+str(months)+'m',data_type='EWEMBI',time_format='monthly',overwrite=True)



COU.hist_merge()

COU.ensemble_mean()

COU.data_summary()

COU.zip_it('raw_cdd')


