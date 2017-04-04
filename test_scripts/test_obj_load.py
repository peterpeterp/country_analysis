import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/allgemeine_scripte/')
from plot_functions import *
sys.path.append('/Users/peterpfleiderer/Documents/')

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis_obj'] 
except:pass
from country_analysis_obj import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')



os.chdir('/Users/peterpfleiderer/Documents/')

SEN=country_analysis('SEN','Projects/country_analysis/')
SEN.load_from_tar('Projects/country_analysis/SEN.tar.gz')


SEN.display()

SEN.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]})

# SEN.prepare_for_download()

SEN._DATA[42].plot_map(period='ref')


#SEN._DATA[44].plot_map(period='2030s-ref')


fig,axes=plt.subplots(rows=2,cols=5,figsize=(10,10))
axes=axes.flatten()
count=0

for 