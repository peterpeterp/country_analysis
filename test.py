# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')
try:del sys.modules['country_analysis']
except:pass
import country_analysis; reload(country_analysis)
sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')


iso='BEN'
os.chdir('/Users/peterpfleiderer/Documents/Projects/country_analysis/')
# data will be stored in working_directory
COU=country_analysis.country_analysis(iso,working_directory='data/'+iso+'/',additional_tag='')
