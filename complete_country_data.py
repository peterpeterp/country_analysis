import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis_obj'] 
except:pass
from country_analysis_obj import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_analysis/')


country_iso='SEN'
COU=country_analysis(country_iso,'/p/projects/tumble/carls/shared_folder/country_analysis/'+country_iso+'/')
COU.load_data()

#COU.hist_merge()

COU.ensemble_mean()

COU.data_summary()

COU.zip_it()


