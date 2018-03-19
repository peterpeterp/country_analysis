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
parser.add_argument("--country",'-c', help="country iso",required=True)
parser.add_argument("--overwrite",'-o', help="overwrite output files",action="store_true")
args = parser.parse_args()

if args.overwrite:
    overwrite=True
else:
    overwrite=False


iso=args.country


print iso
os.chdir('/p/projects/tumble/carls/shared_folder/country_analysis/')
COU=country_analysis(iso,'data/'+iso+'/',additional_tag='')

# ##############
# # MASKS
# ##############
COU.create_mask_country('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',overwrite=False)
COU.create_mask_country('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',mask_style='pop1990_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_1990_incrLat.nc',overwrite=False)
COU.create_mask_country('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',mask_style='pop2015_weighted',pop_mask_file='/home/pepflei/CA/masks/population/population_2015_incrLat.nc',overwrite=False)
COU.create_mask_admin('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc','tmp',overwrite=False)

# ##############
# # Zoom
# ##############
COU.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.tmp.dat.nc',var_name='tmp',given_var_name='tas',data_type='CRU',time_format='monthly',overwrite=False)
COU.country_zoom('/p/projects/elis/CRUDATA_TS3_23/cru_ts3.23.1901.2014.pre.dat.nc',var_name='pre',given_var_name='pr',data_type='CRU',time_format='monthly',overwrite=False)


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
