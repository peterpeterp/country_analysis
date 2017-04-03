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

import seaborn as sns   




SEN=country_analysis('SEN','Projects/country_analysis/')
SEN.load_from_tar('Projects/country_analysis/SEN.tar.gz')

#SEN.ensemble_mean()


# # plot
# fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(6,5))
# axes.set_color_cycle(sns.color_palette("Set3",len(SEN._regions)))
# axes.set_axis_bgcolor('white')
# for name in SEN._regions:
# 	SEN.plot_transient(meta_data=['pr','CMIP5','rcp8.5','hadgem2-es'],mask_style='lat_weighted',region=name.encode('utf-8'),ax=axes,running_mean=240,show=False,ylabel='precipitation [mm]',label=name)

# box = axes.get_position()
# axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.show()

# compute period average
# SEN.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]})

# var='pr'
# for rcp in ['rcp2.6','rcp8.5']:
# 	SEN._data[var]['CMIP5'][rcp]['ensemble_mean']={'period':{}}
# 	for period in ['ref','2030s','2040s','2030s-ref','2040s-ref']:
# 		SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period]=SEN._data[var]['CMIP5'][rcp]['hadgem2-es']['period']['ref'].copy()*0		
# 		SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['lat']=SEN._data[var]['CMIP5'][rcp]['hadgem2-es']['lat'].copy()
# 		SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['lon']=SEN._data[var]['CMIP5'][rcp]['hadgem2-es']['lon'].copy()
# 		for model in SEN._data[var]['CMIP5'][rcp].keys():
# 			SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period]+=SEN._data[var]['CMIP5'][rcp][model]['period'][period]
# 		SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period]/=5


# 	SEN._data[var]['CMIP5'][rcp]['agreement']={}
# 	for period in ['2030s-ref','2040s-ref']:
# 		if period!='ref':
# 			SEN._data[var]['CMIP5'][rcp]['agreement'][period]=SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period].copy()*0
# 			for model in SEN._data[var]['CMIP5'][rcp].keys():
# 				if model not in ['ensemble_mean','agreement']:
# 					SEN._data[var]['CMIP5'][rcp]['agreement'][period][np.where(np.sign(SEN._data[var]['CMIP5'][rcp][model]['period'][period])==np.sign(SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period]))]+=1
# 			SEN._data[var]['CMIP5'][rcp]['agreement'][period][SEN._data[var]['CMIP5'][rcp]['agreement'][period]>3]=np.nan
# 			SEN._data[var]['CMIP5'][rcp]['agreement'][period][np.isnan(SEN._data[var]['CMIP5'][rcp]['agreement'][period])==False]=0.5
# 			SEN._data[var]['CMIP5'][rcp]['agreement'][period][np.ma.getmask(SEN._data[var]['CMIP5'][rcp]['ensemble_mean']['period'][period])]=np.nan


# fig, axes = plt.subplots(nrows=2, ncols=2,figsize=(7,7))
# axes=axes.flatten()
# count=0

# for period in ['2030s-ref','2040s-ref']:
# 	for rcp in ['rcp2.6','rcp8.5']:
# 		im=SEN.plot_map(meta_data=['pr','CMIP5','rcp8.5','ensemble_mean'],period=period,source='_data',ax=axes[count],title=period+rcp,show=False,color_bar=False,color_palette=plt.cm.RdYlBu,color_range=[-10,10])
# 		count+=1
    
    
# cbar_ax=fig.add_axes([0.1,0.05,0.8,0.1])
# cbar_ax.axis('off')
# cb=fig.colorbar(im,orientation='horizontal',label='monthly rx5 [mm]')






# # custom plots
# fig, axes = plt.subplots(nrows=3, ncols=5,figsize=(10,7))
# axes=axes.flatten()
# count=0
# for model in SEN._data['pr']['CMIP5']['rcp8.5'].keys():
#     im=SEN.plot_map(meta_data=['pr','CMIP5','rcp8.5',model],period='ref',source='_data',ax=axes[count],title=model,show=False,color_bar=False,color_palette=plt.cm.YlGnBu,color_range=[25,65])
#     if model == 'ipsl-cm5a-lr': axes[count].set_ylabel('ref')
#     count+=1

# for model in SEN._data['pr']['CMIP5']['rcp8.5'].keys():
#     im=SEN.plot_map(meta_data=['pr','CMIP5','rcp8.5',model],period='2030s-ref',source='_data',ax=axes[count],title=model,show=False,color_bar=False,color_palette=plt.cm.YlGnBu,color_range=[25,65])
#     if model == 'ipsl-cm5a-lr': axes[count].set_ylabel('2030s')
#     count+=1
    
# for model in SEN._data['pr']['CMIP5']['rcp8.5'].keys():
#     im=SEN.plot_map(meta_data=['pr','CMIP5','rcp8.5',model],period='2040s-ref',source='_data',ax=axes[count],title=model,show=False,color_bar=False,color_palette=plt.cm.YlGnBu,color_range=[25,65])
#     if model == 'ipsl-cm5a-lr': axes[count].set_ylabel('2040s')
#     count+=1
    
# cbar_ax=fig.add_axes([0.1,0.05,0.8,0.1])
# cbar_ax.axis('off')
# cb=fig.colorbar(im,orientation='horizontal',label='monthly rx5 [mm]')

# #plt.suptitle('\n\n\nReference Period')
# #plt.tight_layout()
# #plt.show()


# fig, axes = plt.subplots(nrows=2, ncols=5,figsize=(7,7))
# axes=axes.flatten()
# count=0

# for model in SEN._data['pr']['CMIP5']['rcp2.6'].keys():
#     im=SEN.plot_map(meta_data=['pr','CMIP5','rcp2.6',model],period='2030s-ref',source='_data',ax=axes[count],title=model,show=False,color_bar=False,color_palette=plt.cm.RdYlBu,color_range=[-10,10])
#     if model == 'ipsl-cm5a-lr': axes[count].set_ylabel('2030s')
#     count+=1
    
# for model in SEN._data['pr']['CMIP5']['rcp2.6'].keys():
#     im=SEN.plot_map(meta_data=['pr','CMIP5','rcp2.6',model],period='2040s-ref',source='_data',ax=axes[count],title=model,show=False,color_bar=False,color_palette=plt.cm.RdYlBu,color_range=[-10,10])
#     if model == 'ipsl-cm5a-lr': axes[count].set_ylabel('2040s')
#     count+=1
    
# cbar_ax=fig.add_axes([0.1,0.05,0.8,0.1])
# cbar_ax.axis('off')
# cb=fig.colorbar(im,orientation='horizontal',label='monthly rx5 [mm]')

# plt.savefig('SEN_2.6')




