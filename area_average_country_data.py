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

country_iso='SEN'

COU=country_analysis(country_iso,'country_analysis/data/'+country_iso+'/')
COU.load_data()



if COU.selection(['year_RX5'])<2:
	for data in COU.selection(['RX5']): data.year_max('year_RX5')
	for data in COU.selection(['year_RX5','ensemble_mean']): COU._DATA.remove(data)
	COU.ensemble_mean()

COU.summary()	


COU.area_average('lat_weighted',overwrite=False)
# COU.area_average('pop2015_weighted',overwrite=False)
# COU.area_average('pop1990_weighted',overwrite=False)


COU.zip_it()


