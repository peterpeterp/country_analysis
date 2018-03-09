# -*- coding: utf-8 -*-


import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

try:
	basepath='/Users/peterpfleiderer/Documents/Projects/'
	os.chdir(basepath)
except:
	basepath='/p/projects/tumble/carls/shared_folder/'
	os.chdir(basepath)


sys.path.append(basepath+'country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis']
except:pass
from country_analysis import country_analysis
sys.path.append(basepath)

os.chdir(basepath)


all_isos=['AGO', 'DZA', 'EGY', 'GNQ', 'BEN', 'NGA', 'NER', 'ZWE', 'NAM', 'GNB', 'SWZ', 'GHA', 'COG', 'SLE', 'ETH', 'COM', 'ERI', 'CPV', 'LBR', 'LBY', 'LSO', 'UGA', 'RWA', 'SOM', 'MDG', 'CMR', 'TZA', 'BWA', 'SEN', 'TCD', 'GAB', 'BFA', 'MWI', 'MOZ', 'MRT', 'GMB', 'MLI', 'BDI', 'STP', 'DJI', 'GIN', 'ESH', 'KEN', 'MAR', 'COD', 'ZMB', 'ZAF', 'TGO', 'TUN', 'CAF', 'SSD', 'SDN', 'CIV']

try:
	job_id=int(os.environ["SLURM_ARRAY_TASK_ID"])
	print 'job_id=',job_id
	print job_id
	isos=all_isos[(job_id*6):((job_id+1)*6)]
	print isos

except:
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--country",'-c', help="country iso",required=False)
	parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
	args = parser.parse_args()

	if args.overwrite:
	    overwrite=True
	else:
	    overwrite=False

	if args.country is not None:
		isos=[args.country]
	else:
		isos=['CPV']


for country_iso in isos:
	# try:
	print country_iso

	COU=country_analysis(country_iso,'country_analysis/data/'+country_iso+'/')
	COU.load_data(quiet=False,load_region_polygons=False)



	if len(COU.selection(['year_RX5']))<2:
		for data in COU.selection(['RX5']): data.year_max('year_RX5')
		for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
		COU.ensemble_mean()

	COU.summary()


	COU.area_average('lat_weighted',overwrite=True)
	# COU.area_average('pop2015_weighted',overwrite=False)
	# COU.area_average('pop1990_weighted',overwrite=False)


	COU.zip_it()

	# except Exception,e:
	# 	print str(e)
	# 	print 'issues with '+country_iso
