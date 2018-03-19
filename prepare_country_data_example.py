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
	iso=args.country
else:
	return 0

print iso
os.chdir('/p/projects/tumble/carls/shared_folder/country_analysis/')
COU=country_analysis(iso,'data/'+iso+'/',additional_tag='')

# ##############
# # MASKS
# ##############
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas',overwrite=True)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=True)
COU.create_mask_country('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=True)
COU.create_mask_admin('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4','tas',overwrite=True)

# ##############
# # EWEMBI
# ##############
COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4',var_name='tas',data_type='EWEMBI',time_format='monthly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_TXx_EWEMBI_1979-2014.nc4',var_name='tasmax',given_var_name='TXx',data_type='EWEMBI',time_format='monthly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_pr_EWEMBI_1979-2014.nc4',var_name='pr',data_type='EWEMBI',time_format='monthly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX5_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX5',data_type='EWEMBI',time_format='monthly',overwrite=True)

COU.country_zoom('/p/projects/tumble/carls/shared_folder/EWEMBI/mon_year/mon_RX1_EWEMBI_1979-2014.nc4',var_name='pr',given_var_name='RX1',data_type='EWEMBI',time_format='monthly',overwrite=True)


COU.hist_merge()

COU.ensemble_mean()

COU.summary()

if len(COU.selection(['year_RX5']))<2:
	for data in COU.selection(['RX5']): data.year_max('year_RX5')
	for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
	COU.ensemble_mean()

COU.summary()


COU.area_average('lat_weighted',overwrite=True)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)

COU.zip_it()
