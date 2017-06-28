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


all_isos=['GRD']

# ##############
# # Pre-Prepare
# ##############
for country_iso in all_isos:
	print country_iso
	COU=country_analysis(country_iso,'/p/projects/tumble/carls/shared_folder/country_analysis/data/'+country_iso+'/',additional_tag='')

	os.chdir(COU._working_directory)
	os.system('wget biogeo.ucdavis.edu/data/gadm2.8/shp/'+country_iso+'_adm_shp.zip')
	os.system('mkdir '+country_iso+'_adm_shp')
	os.system('ls')
	os.chdir(COU._working_directory+country_iso+'_adm_shp')
	os.system('unzip ../'+country_iso+'_adm_shp.zip')


try:
	job_id=int(os.environ["SLURM_ARRAY_TASK_ID"])
	print 'job_id=',job_id
	print job_id
	isos=all_isos[(job_id*6):((job_id+1)*6)]
	print isos

except:
	isos=['GRD']	



for country_iso in isos:
	print country_iso
	if True:
		COU=country_analysis(country_iso,'/p/projects/tumble/carls/shared_folder/country_analysis/data/'+country_iso+'/',additional_tag='')

		# ##############
		# # MASKS
		# ##############

		COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',overwrite=False)
		COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
		COU.create_mask_country('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
		COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_HadGEM2-ES_historical_0-2004_BC.nc4','tas',COU._working_directory+country_iso+'_adm_shp/'+country_iso+'_adm1',overwrite=False)


		# ##############
		# #CMIP5_BC
		# ##############

		#tn90p 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/TN90P/mon_tn90p_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='warm_nights_percent_wrt_90th_percentile_of_reference_period',given_var_name='tn90p',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)


		#tx90p 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/TX90P/mon_tx90p_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='very_warm_days_percent_wrt_90th_percentile_of_reference_period',given_var_name='tx90p',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#r95p 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/R95P/mon_r95p_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='very_wet_days_wrt_95th_percentile_of_reference_period',given_var_name='r95p',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#tas 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/tas/mon_tas_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='tas',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#TXx 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/TXx/mon_TXx_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='tasmax',given_var_name='TXx',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#pr 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/pr/mon_pr_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='pr',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#RX5 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/RX5/mon_RX5_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='pr',given_var_name='RX5',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)

		#RX1 
		all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/RX1/mon_RX1_*')
		for in_file in all_files:
			model=in_file.split('/')[-1].split('_')[2]
			rcp=in_file.split('/')[-1].split('_')[3]
			print rcp,model,in_file
			COU.country_zoom(in_file,var_name='pr',given_var_name='RX1',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)


		#SPEI
		for months in [1,3,6,12]:
			all_files=glob.glob('/p/projects/tumble/carls/shared_folder/CMIP5_monthly/SPEI/SPEI_*'+str(months)+'m_*penman*')
			for in_file in all_files:
				model=in_file.split('/')[-1].split('_')[-2]
				rcp=in_file.split('/')[-1].split('_')[-1].split('.')[0]
				print rcp,model,in_file
				COU.country_zoom(in_file,var_name='spei',given_var_name='SPEI_'+str(months)+'m_penman',data_type='CMIP5_BC',scenario=rcp,model=model,time_format='monthly',overwrite=False)




		# ##############
		# # EWEMBI
		# ##############
		# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4',var_name='tas',data_type='EWEMBI',time_format='monthly',overwrite=True)

		# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_TXx_EWEMBI_1979-2014.nc4',var_name='tasmax',given_var_name='TXx',data_type='EWEMBI',time_format='monthly',overwrite=True)

		# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_pr_EWEMBI_1979-2014.nc4',var_name='pr',data_type='EWEMBI',time_format='monthly',overwrite=True)

		# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX5',data_type='EWEMBI',time_format='monthly',overwrite=True)

		# COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX1_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX1',data_type='EWEMBI',time_format='monthly',overwrite=True)


		COU.hist_merge()

		COU.ensemble_mean()

		COU.summary()
		
		if COU.selection(['year_RX5'])<2:
			for data in COU.selection(['RX5']): data.year_max('year_RX5')
			for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
			COU.ensemble_mean()

		# COU.summary()

		# COU.area_average('lat_weighted',overwrite=True)

		# COU.zip_it()
		#COU.zip_it('raw'+COU._additional_tag)

	# except:
	# 	print 'issues with '+country_iso






