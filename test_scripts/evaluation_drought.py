import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

#################
# plot settings
#################
from mpl_toolkits.basemap import Basemap, cm
import mpl_toolkits.basemap
import matplotlib.pylab as plt
import matplotlib as mpl
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from matplotlib.colors import ListedColormap

sys.path.append('/Users/peterpfleiderer/Documents/Scripts/allgemeine_scripte/')
from plot_functions import *
sys.path.append('/Users/peterpfleiderer/Documents/')

risk = make_colormap([col_conv('green'), col_conv('white'), 0.2, col_conv('white'), col_conv('yellow'), 0.4, col_conv('yellow'), col_conv('orange'), 0.6, col_conv('orange'), col_conv('red'), 0.8, col_conv('red'), col_conv('violet')])
month_color = mpl.colors.ListedColormap(sns.color_palette("cubehelix", 12))

#################
# load object
#################

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/country_analysis_scripts/')
try:del sys.modules['country_analysis'] 
except:pass
from country_analysis import country_analysis
sys.path.append('/Users/peterpfleiderer/Documents/')

GHA=country_analysis('GHA','Projects/country_analysis/')

pkl_file = open(GHA._working_directory+'/'+GHA._iso+'_masks.pkl', 'rb')
GHA._masks = pickle.load(pkl_file)	;	pkl_file.close()  

pkl_file = open(GHA._working_directory+'/'+GHA._iso+'_data.pkl', 'rb')
GHA._data = pickle.load(pkl_file)	;	pkl_file.close()  

pkl_file = open(GHA._working_directory+'/'+GHA._iso+'_meta.pkl', 'rb')
GHA._meta = pickle.load(pkl_file)	;	pkl_file.close()  

#################
# extreme drought
#################
periods={'ref':[1985,2005],'2030s':[2024,2044],'2040s':[2034,2054]}

# ensemble mean
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['SPEI_12m']['CMIP5'][rcp]
	tmp['ensemble_mean']=tmp[tmp.keys()[0]].copy()
	tmp['ensemble_mean']['data']=tmp[tmp.keys()[0]]['data'].copy()*0
	for model in tmp.keys():
		if model!='ensemble_mean':
			tmp['ensemble_mean']['data']+=tmp[model]['data'].copy()
	tmp['ensemble_mean']['data']/=5
	GHA._meta.append(['SPEI_12m','CMIP5',rcp,'ensemble_mean'])

GHA.period_averages(periods=periods)

# exposure
GHA._data['SPEI_12m_expo-2']={}
for meta_list in GHA._meta:
	if meta_list[0]=='SPEI_12m':
		tmp_in = GHA._data['SPEI_12m']
		tmp_out = GHA._data['SPEI_12m_expo-2']
		for meta_info in meta_list[1:]:
			tmp_in=tmp_in[meta_info]
			if meta_info not in tmp_out:tmp_out[meta_info]={}
			tmp_out=tmp_out[meta_info]

		for period in periods:
			relevant_data_indices=np.where((tmp_in['year']>periods[period][0]) & (tmp_in['year']<=periods[period][1]))[0]
			tmp_out[period]=tmp_in['period']['ref'][:,:].copy()*np.nan

			for y in range(tmp_out[period].shape[0]):
				for x in range(tmp_out[period].shape[1]):
					if len(relevant_data_indices)>0: tmp_out[period][y,x]=float(len(np.where(tmp_in['data'][relevant_data_indices,y,x].flatten()<=-2)[0]))/len(relevant_data_indices)


# ensemble mean
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
	tmp['ensemble_mean']={}
	for period in periods:
		tmp['ensemble_mean'][period]=tmp['hadgem2-es'][period].copy()*0
		for model in tmp.keys():
			if model!='ensemble_mean':
				tmp['ensemble_mean'][period]+=tmp[model][period].copy()
		tmp['ensemble_mean'][period]/=5



# period diff
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
	for model in tmp.keys():
		if model not in ['agreement']:
			for period in periods:
				if period!='ref':
					tmp[model][period+'-'+'ref']=tmp[model][period]-tmp[model]['ref']


# model agremment
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
	tmp['agreement']={}
	for period in tmp['ensemble_mean'].keys():
		if '-' in period:	
			tmp['agreement'][period]=tmp['ensemble_mean'][period].copy()*0
			for model in tmp.keys():
				if model not in ['ensemble_mean','agreement']:
					tmp['agreement'][period][np.where(np.sign(tmp[model][period])==np.sign(tmp['ensemble_mean'][period]))]+=1
			tmp['agreement'][period][tmp['agreement'][period]>3]=np.nan
			tmp['agreement'][period][np.isnan(tmp['agreement'][period])==False]=0.5
			tmp['agreement'][period][np.ma.getmask(tmp['ensemble_mean'][period])]=np.nan



##############
# plot final result
##############
fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(4,4))
count=0
lon=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lon']
lat=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lat']
for rcp in ['rcp2.6','rcp8.5']:
	for period in ['2030s-ref','2040s-ref']:
		tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
		ax,im=plot_map(axes.flatten()[count],lon,lat,tmp['ensemble_mean'][period]*100,color_type=risk,color_range=[-7,28],color_label=None,subtitle='',grey_area=tmp['agreement'][period])
		if rcp=='rcp2.6':ax.set_title(period)
		if period=='2030s-ref':ax.set_ylabel(rcp)
		count+=1
		if period=='2040s-ref':
			axes.flatten()[count].axis('off')
			count+=1
cbar_ax=fig.add_axes([0.45,0.2,0.3,0.6])
cbar_ax.axis('off')
cb=fig.colorbar(im,orientation='vertical',label='change in frequency of extreme dry events \n [percentage of affected months]')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
plt.savefig(GHA._working_directory+'/plots/spei_X_sum.png')

####################
# plot reference period
####################
tmp=GHA._data['SPEI_12m_expo-2']
fig,axes=plt.subplots(nrows=1,ncols=8,figsize=(8,3))
ax,im=plot_map(axes.flatten()[0],GHA._data['SPEI_12m']['CRU']['lon'],GHA._data['SPEI_12m']['CRU']['lat'],GHA._data['SPEI_12m_expo-2']['CRU']['ref']*100,color_type=plt.cm.YlOrBr,color_range=[0,15],color_label=None,subtitle='CRU')
ax,im=plot_map(axes.flatten()[1],GHA._data['SPEI_12m']['NCEP']['lon'],GHA._data['SPEI_12m']['NCEP']['lat'],GHA._data['SPEI_12m_expo-2']['NCEP']['ref']*100,color_type=plt.cm.YlOrBr,color_range=[0,15],color_label=None,subtitle='NCEP',limits=[-3.25,1.25,4.75,11.25])

lon=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lon']
lat=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lat']
ax,im=plot_map(axes.flatten()[2],lon,lat,GHA._data['SPEI_12m_expo-2']['CMIP5']['rcp2.6']['ensemble_mean']['ref']*100,color_type=plt.cm.YlOrBr,color_range=[0,15],color_label=None,subtitle='ens. mn.')
count=3
tmp=GHA._data['SPEI_12m_expo-2']['CMIP5']['rcp2.6']
for model in tmp.keys():
	if model not in ['ensemble_mean','agreement']:
		ax,im=plot_map(axes.flatten()[count],lon,lat,GHA._data['SPEI_12m_expo-2']['CMIP5']['rcp2.6'][model]['ref']*100,color_type=plt.cm.YlOrBr,color_range=[0,15],color_label=None,subtitle='')
		ax.set_title(model,fontsize=8)
		count+=1
cbar_ax=fig.add_axes([0,0.17,1,0.3])
cbar_ax.axis('off')
cb=fig.colorbar(im,orientation='horizontal',label='ratio of months affected by SPEI -2 [%]')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

plt.suptitle('\n\nRefernce Period 1986-2005')
plt.savefig(GHA._working_directory+'/plots/spei_X_ref.png')

###############
# plot exposure
###############
lon=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lon']
lat=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lat']

for rcp in ['rcp2.6','rcp8.5']:
	fig,axes=plt.subplots(nrows=3,ncols=7,figsize=(10,6))
	count=0
	tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
	for period in ['ref','2030s','2040s']:
		ax,im=plot_map(axes.flatten()[count],lon,lat,tmp['ensemble_mean'][period]*100,color_type=plt.cm.YlOrBr,color_range=[0,20],color_label=None,subtitle='')
		if period=='ref':ax.set_title('ens. mean')
		ax.set_ylabel(period)
		count+=1
		for model in tmp.keys():
			if model not in ['agreement','ensemble_mean']:
				ax,im=plot_map(axes.flatten()[count],lon,lat,tmp[model][period]*100,color_type=plt.cm.YlOrBr,color_range=[0,20],color_label=None,subtitle='')
				if period=='ref':ax.set_title(model)
				count+=1

		ax=axes.flatten()[count]
		ax.axis('off')
		count+=1 

	cbar_ax=fig.add_axes([0.8,0.2,0.1,0.6])
	cbar_ax.axis('off')
	cb=fig.colorbar(im,orientation='vertical',label='ratio of months affected by heat extremes [%]')
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()
	plt.savefig(GHA._working_directory+'/plots/spei_X_'+rcp+'.png')

###############
# plot exposure diff
###############
lon=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lon']
lat=GHA._data['SPEI_12m']['CMIP5'][rcp]['ensemble_mean']['lat']

for rcp in ['rcp2.6','rcp8.5']:
	fig,axes=plt.subplots(nrows=2,ncols=8,figsize=(10,4))
	count=0
	tmp=GHA._data['SPEI_12m_expo-2']['CMIP5'][rcp]
	for period in ['2030s-ref','2040s-ref']:
		ax,im=plot_map(axes.flatten()[count],lon,lat,tmp['ensemble_mean'][period]*100,color_type=risk,color_range=[-7,28],color_label=None,subtitle='')
		if period=='2030s-ref':ax.set_title('ens. mean')
		ax.set_ylabel(period)
		count+=1
		ax,im=plot_map(axes.flatten()[count],lon,lat,tmp['ensemble_mean'][period]*100,color_type=risk,color_range=[-7,28],color_label=None,subtitle='',grey_area=tmp['agreement'][period])
		if period=='2030s-ref':ax.set_title('ens. mean')
		count+=1
		for model in tmp.keys():
			if model not in ['agreement','ensemble_mean']:
				ax,im=plot_map(axes.flatten()[count],lon,lat,tmp[model][period]*100,color_type=risk,color_range=[-7,28],color_label=None,subtitle='')
				if period=='2030s-ref':ax.set_title(model)
				count+=1

		ax=axes.flatten()[count]
		ax.axis('off')
		count+=1 

	cbar_ax=fig.add_axes([0.8,0.2,0.1,0.6])
	cbar_ax.axis('off')
	cb=fig.colorbar(im,orientation='vertical',label='change in frequency of extreme dry events \n [percentage of affected months]')
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()
	plt.savefig(GHA._working_directory+'/plots/spei_X_diff_'+rcp+'.png')



