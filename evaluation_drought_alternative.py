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
# drought risk
#################

# annual min
GHA._data['SPEI_12m_yrMin']={}
for meta_list in GHA._meta:
	if meta_list[0]=='SPEI_12m':
		tmp_in = GHA._data['SPEI_12m']
		tmp_out = GHA._data['SPEI_12m_yrMin']
		for meta_info in meta_list[1:]:
			tmp_in=tmp_in[meta_info]
			if meta_info not in tmp_out:tmp_out[meta_info]={}
			tmp_out=tmp_out[meta_info]

		GHA._meta.append(['SPEI_12m_yrMin']+meta_list[1:])
		for key in tmp_in.keys():
			if key not in ['data','year','month','time']:
				tmp_out[key]=tmp_in[key]

		tmp_out['year']=np.array(sorted(list(set(tmp_in['year']))))	;	nYr=len(tmp_out['year'])
		tmp_out['month']=np.zeros([nYr])
		tmp_out['time']=np.zeros([nYr])
		tmp_out['data']=tmp_in['data'][range(nYr),:,:].copy()
		for yr in range(nYr):
			tmp_out['data'][yr,:,:]=np.min(tmp_in['data'][np.where(tmp_in['year']==tmp_out['year'][yr])[0],:,:],axis=0)

		

periods={'ref':[1985,2005],'2030s':[2024,2044],'2040s':[2034,2054]}

# ensemble mean
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['SPEI_12m_yrMin']['CMIP5'][rcp]
	tmp['ensemble_mean']=tmp[tmp.keys()[0]].copy()
	tmp['ensemble_mean']['data']=tmp[tmp.keys()[0]]['data'].copy()*0
	for model in tmp.keys():
		if model!='ensemble_mean':
			tmp['ensemble_mean']['data']+=tmp[model]['data'].copy()
	tmp['ensemble_mean']['data']/=5
	GHA._meta.append(['SPEI_12m_yrMin','CMIP5',rcp,'ensemble_mean'])


GHA.period_averages(periods=periods)



####################
# plot reference period
####################
tmp=GHA._data['SPEI_12m_yrMin']
fig,axes=plt.subplots(nrows=1,ncols=8,figsize=(8,3))
ax,im=plot_map(axes.flatten()[0],GHA._data['SPEI_12m']['CRU']['lon'],GHA._data['SPEI_12m']['CRU']['lat'],tmp['CRU']['period']['ref'],color_type=plt.cm.spring,color_range=[-2,1],color_label=None,subtitle='CRU')
ax,im=plot_map(axes.flatten()[1],GHA._data['SPEI_12m']['NCEP']['lon'],GHA._data['SPEI_12m']['NCEP']['lat'],tmp['NCEP']['period']['ref'],color_type=plt.cm.spring,color_range=[-2,1],color_label=None,subtitle='NCEP',limits=[-3.25,1.25,4.75,11.25])

lon=GHA._data['SPEI_12m']['CMIP5'][rcp]['hadgem2-es']['lon']
lat=GHA._data['SPEI_12m']['CMIP5'][rcp]['hadgem2-es']['lat']
ax,im=plot_map(axes.flatten()[2],lon,lat,tmp['CMIP5']['rcp2.6']['ensemble_mean']['period']['ref'],color_type=plt.cm.spring,color_range=[-2,1],color_label=None,subtitle='ens. mn.')
count=3
tmp=tmp['CMIP5']['rcp2.6']
for model in tmp.keys():
	if model not in ['ensemble_mean','agreement']:
		ax,im=plot_map(axes.flatten()[count],lon,lat,tmp[model]['period']['ref'],color_type=plt.cm.spring,color_range=[-2,1],color_label=None,subtitle='')
		ax.set_title(model,fontsize=8)
		count+=1
cbar_ax=fig.add_axes([0,0.17,1,0.3])
cbar_ax.axis('off')
cb=fig.colorbar(im,orientation='horizontal',label='annual minimal SPEI value')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

plt.suptitle('\n\nRefernce Period 1986-2005')
plt.savefig(GHA._working_directory+'/plots/spei_yrM_ref.png')




GHA.average('pop2015_weighted')
for rcp in GHA._transcient['SPEI_12m_yrMin']['CMIP5'].keys():
	tmp=GHA._transcient['SPEI_12m_yrMin']['CMIP5'][rcp]

	fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(4,4))
	relevant_years=range(30,100)
	for model in tmp.keys():

		if model!='ensemble_mean':plt.plot(tmp[model]['year'][relevant_years],tmp[model]['pop2015_weighted'][relevant_years],linestyle='--',label=model)
		if model=='ensemble_mean':plt.plot(tmp[model]['year'][relevant_years],tmp[model]['pop2015_weighted'][relevant_years],linestyle='-',label=model)

	plt.legend(loc='best')
	plt.savefig(GHA._working_directory+'/plots/spei_yrM_'+rcp+'.png')





