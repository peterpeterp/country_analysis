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





SEN=country_analysis('SEN','Projects/country_analysis/')
SEN.load_from_tar('Projects/country_analysis/SEN.tar.gz')

#GHA.load_from_tar('Projects/country_analysis/SEN.tar.gz')







# GHA=country_analysis('GHA','Projects/country_analysis/')
# GHA.create_mask_admin('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/GHA_adm_shp/GHA_adm1',overwrite=False)

# GHA.create_mask_country('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',overwrite=False)
# GHA.create_mask_country('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc',overwrite=False)
# GHA.prepare_for_download()

#GHA.create_mask_admin('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/GHA_adm_shp/GHA_adm1')

# # GHA.plot_map(meta_data=['94x192','pop1990_weighted','GHA'],source='_masks')


# # GHA.plot_map(meta_data=['94x192','lat_weighted','Ashanti'],source='_masks')
# GHA.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['SPEI_1m','NCEP'])


# GHA.average('lat_weighted',meta_data=['SPEI_1m','NCEP'],overwrite=False)

# # # plot
# # fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(4,4))
# # for name in GHA._regions:
# #     GHA.plot_transient(meta_data=['SPEI_1m','NCEP'],mask_style='lat_weighted',region=name,ax=axes,running_mean=240,show=False,ylabel='SPEI',label=name)

# # plt.show()





