import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/allgemeine_scripte/')
from plot_functions import *
sys.path.append('/Users/peterpfleiderer/Documents/')

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

GHA=country_analysis('GHA','Projects/country_analysis/')

GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries')
GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')
GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')

GHA.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries')
GHA.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')


in_files=glob.glob('data/raw/mon_rx5/CMIP5/*/mon_rx5_*_1950-2099.nc4')
for in_file in in_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'mon_rx5',meta_data=['rx5','CMIP5',rcp,model])

GHA.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['SPEI_1m','NCEP'])

GHA.country_zoom('data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','pre',meta_data=['pr','CRU'])

#GHA.plot_map(meta_data=['pr','CRU'])
GHA.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]},meta_data=['rx5','CMIP5','rcp2.6'])




#GHA.plot_map(meta_data=['94x192','lat_weighted'],source='_masks',limits=[-4,2,3,12])
#GHA.plot_map(meta_data=['pr','CRU'],period='ref')
#GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],period='ref',source='_data')
#GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],source='_data')
#GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],period='ref',source='_data',limits=[-4,2,4,12])


#GHA.average('pop2015_weighted',meta_data=['rx5','CMIP5','rcp2.6'])
#GHA.plot_transcient(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],mask_style='pop2015_weighted')
#GHA.plot_transcient(meta_data=['rx5','CMIP5','rcp8.5','ensemble_mean'],mask_style='pop2015_weighted')


year=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['year']
time=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['year']
month=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['month']
lon=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['lon']
lat=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['lat']
data=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['data']


