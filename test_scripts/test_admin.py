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

GHA.create_mask_admin('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/GHA_adm_shp/GHA_adm1')
GHA.create_mask_country('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries')
GHA.create_mask_country('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')


GHA.plot_map(meta_data=['94x192','pop1990_weighted','GHA'],source='_masks')
