import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/p/projects/tumble/carls/shared_folder/country_anaylsis/')

GRL=country_analysis('GRL','/p/projects/tumble/carls/shared_folder/country_anaylsis/')


GRL.create_mask('/home/pepflei/CA/data/txx_erainterim.1979-2016.nc','t2m','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')
GRL.create_mask('/home/pepflei/CA/data/txx.20CR.1851-2014.nc','air','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')

GRL.create_mask('/home/pepflei/CA/data/H2_TXx_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc','Ann','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')

GRL.create_mask('/home/pepflei/CA/data/GHCND_TXx_1951-2016_RegularGrid_global_2.5x2.5deg_LSmask.nc','Ann','/home/pepflei/CA/masks/shapefiles/world/ne_50m_admin_0_countries')



