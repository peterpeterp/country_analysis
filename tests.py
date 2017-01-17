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

GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries')
GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop1990_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')
GHA.create_mask('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')

GHA.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries')
GHA.create_mask('data/raw/SPEI/CMIP5/spei_hadgem2-es_rcp2.6_1950-2099_1m.nc','SPEI','masks/shapefiles/world/ne_50m_admin_0_countries',mask_style='pop2015_weighted',pop_mask_file='masks/population/population_1990_incrLat.nc')


in_files=glob.glob('data/raw/mon_rx5/CMIP5/*/mon_rx5_*_1950-2099.nc4')
for in_file in in_files:
	rcp=in_file.split('_')[-2]
	model=in_file.split('_')[-3]
	GHA.country_zoom(in_file,'mon_rx5',meta_data=['rx5','CMIP5',rcp,model])

GHA.country_zoom('data/raw/SPEI/NCEP/SPEI_ncep_1948-2014_1m.nc','SPEI',meta_data=['SPEI_1m','NCEP'])

GHA.country_zoom('data/raw/cru/cru_ts3.23.1901.2014.pre.dat.nc','pre',meta_data=['pr','CRU'])

#GHA.plot_map(meta_data=['pr','CRU'])
GHA.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]})






fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,3))

Z=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['period']['ref']
lat=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['lat']
lon=GHA._data['rx5']['CMIP5']['rcp8.5']['gfdl-esm2m']['lon']
limits=[np.min(lon)-0.25,np.max(lon)+0.25,np.min(lat)-0.25,np.max(lat)+0.25]

print limits
if lat[0]>lat[1]:Z=Z[::-1,:]

m = Basemap(ax=axes,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
m.drawmapboundary(fill_color='1.')

im1 = m.imshow(Z,cmap=plt.cm.plasma,interpolation='none',extent=limits)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
plt.show()



















#GHA.plot_map(meta_data=['94x192','lat_weighted'],source='_masks',limits=[-4,2,3,12])
#GHA.plot_map(meta_data=['pr','CRU'],period='ref')
#GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],period='ref',source='_data')
GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],source='_data')
#GHA.plot_map(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],period='ref',source='_data',limits=[-4,2,4,12])


#GHA.average('pop2015_weighted')
#GHA.plot_transcient(meta_data=['rx5','CMIP5','rcp8.5','gfdl-esm2m'],mask_style='pop2015_weighted')
#GHA.plot_transcient(meta_data=['rx5','CMIP5','rcp8.5','ensemble_mean'],mask_style='pop2015_weighted')



# ensemble mean
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['rx5']['CMIP5'][rcp]
	tmp['ensemble_mean']=tmp[tmp.keys()[0]]
	tmp['ensemble_mean']['data']=tmp[tmp.keys()[0]]['data'].copy()*0
	for model in tmp.keys():
		if model!='ensemble_mean':
			tmp['ensemble_mean']['data']+=tmp[model]['data']
	tmp['ensemble_mean']['data']/=5
	GHA._meta.append(['rx5','CMIP5',rcp,'ensemble_mean'])

GHA.period_averages(periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]})

# model agremment
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['rx5']['CMIP5'][rcp]
	tmp['ensemble_mean']['agreement']={}
	for period in tmp['ensemble_mean']['period']:
		if '-' in period:			
			tmp['ensemble_mean']['agreement'][period]=tmp['ensemble_mean']['period'][period].copy()*0
			for model in tmp.keys():
				if model!='ensemble_mean':
					tmp['ensemble_mean']['agreement'][period][np.where(np.sign(tmp[model]['period'][period])==np.sign(tmp['ensemble_mean']['period'][period]))]+=1
			tmp['ensemble_mean']['agreement'][period][tmp['ensemble_mean']['agreement'][period]>3]=np.nan
			tmp['ensemble_mean']['agreement'][period][np.isnan(tmp['ensemble_mean']['agreement'][period])==False]=0.5
			tmp['ensemble_mean']['agreement'][period][np.ma.getmask(tmp['ensemble_mean']['period'][period])]=np.nan



fig, axes = plt.subplots(nrows=2, ncols=6,figsize=(12,4))
count=0
for rcp in ['rcp2.6','rcp8.5']:
	tmp=GHA._data['rx5']['CMIP5'][rcp]
	ax,im=plot_map(axes.flatten()[count],tmp['ensemble_mean']['lon'],tmp['ensemble_mean']['lat'],tmp['ensemble_mean']['period']['2040s-ref'],color_range=[-10,10],color_label=None,subtitle='',grey_area=tmp['ensemble_mean']['agreement']['2040s-ref'])
	ax.set_ylabel(rcp)
	if rcp=='rcp2.6':ax.set_title(model)
	count+=1
	for model in tmp.keys():
		if model!='ensemble_mean':
			ax,im=plot_map(axes.flatten()[count],tmp[model]['lon'],tmp[model]['lat'],tmp[model]['period']['2040s-ref'],color_range=[-10,10],color_label=None,subtitle='')
			if rcp=='rcp2.6':ax.set_title(model)
			count+=1
cbar_ax=fig.add_axes([0.85,0.2,0.1,0.6])
cbar_ax.axis('off')
cb=fig.colorbar(im,orientation='vertical',label='increase in RX5 (2040s-ref) [mm]')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
plt.savefig(GHA._working_directory+GHA._iso+'/plots/test.png')




