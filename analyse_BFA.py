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


country_iso='BFA'
overwrite=False

print country_iso
COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag='')

os.chdir('/p/projects/tumble/carls/shared_folder/country_analysis/')
COU=country_analysis(country_iso,'data/'+country_iso+'/',additional_tag='')

COU.load_data()
