# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,num2date
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import cartopy.crs as ccrs
import cartopy

iso = 'BFA'
print(iso)
os.chdir('/Users/peterpfleiderer/Projects/country_analysis/')
import country_analysis as country_analysis; reload(country_analysis)

# data will be stored in working_directory
COU=country_analysis.country_analysis(iso,working_directory='data/'+iso+'/',additional_tag='')
grid='13x17_lat_15.25_9.25_lon_-5.75_2.25'

# Load all the data
COU.load_data('BFA.nc4')
plot_mask=np.ma.getdata(COU._masks[grid]['lat_weighted'][iso])

# #This has already been done and saved. No need to run again
# for rcp in ['rcp2p6','rcp4p5','rcp8p5']:
# 	members_tas=COU.classify_ensemble(['tas',rcp])
# 	COU.ensemble_statistic(members_tas,'tas_ensmedian_'+rcp)
# 	members_tas=COU.classify_ensemble(['pr',rcp])
# 	COU.ensemble_statistic(members_tas,'pr_ensmedian_'+rcp)
# 	members_tas=COU.classify_ensemble(['spi3m',rcp])
# 	COU.ensemble_statistic(members_tas,'spi3m_ensmedian_'+rcp)
# COU.save_data('BFA.nc4')


# Here you can define some periods for which a mean is calculated
COU.period_statistic([1986,2006],'ref_mean')
COU.period_statistic([2025,2045],'2030s_mean')
COU.period_statistic([2035,2055],'2040s_mean')

# Make a difference of these period means
COU.period_diff('diff_ref_2040s',ref_period='ref_mean',target_period='2040s_mean')
COU.period_diff('rel_diff_ref_2040s',ref_period='ref_mean',target_period='2040s_mean',relative=True)

for rcp,RCP in zip(['rcp2p6','rcp4p5','rcp8p5'],['RCP26','RCP45','RCP85']):
	###############
	# SPI
	###############

	# define ensemble and specify ensemble median
	# the fields are already loaded, here it's just about the names given to the fields
	all_spi=COU.classify_ensemble(['spi3m',rcp])
	members_spi=all_spi[:-1]
	ensmedian_spi=all_spi[-1]
	print(members_spi,ensmedian_spi)

	# get the frequency of events below a threshold
	COU.period_events_above_thresh([1986,2006],'ref_droughtRisk',threshold=-1.5,names=all_spi,below=True)
	COU.period_events_above_thresh([2025,2045],'2040s_droughtRisk',threshold=-1.5,names=all_spi,below=True)

	# compute the difference
	COU.period_diff('diff_droughtRisk_2040s',ref_period='ref_droughtRisk',target_period='2040s_droughtRisk')

	# get model agreement for this difference
	COU.model_agreement(ensmedian_spi,members_spi,'diff_droughtRisk_2040s','spi3m_ensmedian_agree_'+rcp)

	# plot a map for drought risk
	ret=COU.plot_map(COU._periods[grid]['diff_droughtRisk_2040s'][ensmedian_spi],
				grey_area=COU._periods[grid]['diff_droughtRisk_2040s']['spi3m_ensmedian_agree_'+rcp],
				out_file=COU._working_directory+'/plots/spi_diff_droughtRisk_2040s_ensmedian_'+rcp+'.png',
				add_mask=plot_mask,
				color_palette=plt.cm.RdBu_r,
				color_range=[-2,2],
				color_label='change in precipitation \n deficit risk [%]',
				title='2040s vs 1986-2005 '+RCP)

	# get the frequency of events above a threshold
	COU.period_events_above_thresh([1986,2006],'ref_wetRisk',threshold=1.5,names=all_spi)
	COU.period_events_above_thresh([2025,2045],'2040s_wetRisk',threshold=1.5,names=all_spi)
	COU.period_diff('diff_wetRisk_2040s',ref_period='ref_wetRisk',target_period='2040s_wetRisk')
	COU.model_agreement(ensmedian_spi,members_spi,'diff_wetRisk_2040s','spi3m_ensmedian_agree_'+rcp)

	# plot map again
	ret=COU.plot_map(COU._periods[grid]['diff_wetRisk_2040s'][ensmedian_spi],
				grey_area=COU._periods[grid]['diff_wetRisk_2040s']['spi3m_ensmedian_agree_'+rcp],
				out_file=COU._working_directory+'/plots/spi_diff_wetRisk_2040s_ensmedian_'+rcp+'.png',
				add_mask=plot_mask,
				color_palette=plt.cm.RdBu,
				color_range=[-10,10],
				color_label='change in precipitation \n  surplus risk [%]',
				title='2040s vs 1986-2005 '+RCP)


	# this plots 10 panaels as an overview of the signals for each model
	fig, axes = plt.subplots(nrows=3, ncols=5,figsize=(12,5),subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw = {'height_ratios':[3,3,1]})
	for row,period in zip(range(2),['ref_droughtRisk','2040s_droughtRisk']):
		for member,ax in zip(members_spi,axes[row,:]):
			ret=COU.plot_map(COU._periods[grid][period][member],
						ax=ax,
						add_mask=plot_mask,
						color_palette=plt.cm.Reds,
						color_bar=False,
						color_range=[0,20],
						title=member.split('_')[1])

	axes[0,0].annotate('ref', xy=(-0.2, 0.6),rotation=90, xycoords='axes fraction', fontsize=12)
	axes[1,0].annotate('2040s', xy=(-0.2, 0.6),rotation=90, xycoords='axes fraction', fontsize=12)

	for ax in axes[2,:]:
		ax.outline_patch.set_edgecolor('white');	ax.axis('off')
	cbar_ax=fig.add_axes([0,0.1,1,0.3]);	cbar_ax.axis('off')
	cb=fig.colorbar(ret[1],orientation='horizontal',label='risk of SPI3m < -1.5',ax=cbar_ax)
	# plt.tight_layout()
	plt.savefig(COU._working_directory+'/plots/spi_overwiew_drought_risk_'+rcp+'.png',)

	# this plots 5 panels showing differences for each model
	fig, axes = plt.subplots(nrows=2, ncols=5,figsize=(12,4),subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw = {'height_ratios':[3,1]})
	for member,ax in zip(members_spi,axes[0,:]):
		ret=COU.plot_map(COU._periods[grid]['diff_droughtRisk_2040s'][member],
					ax=ax,
					add_mask=plot_mask,
					color_palette=plt.cm.RdBu,
					color_bar=False,
					color_range=[-5,5],
					title=member.split('_')[1])
	for ax in axes[1,:]:
		ax.outline_patch.set_edgecolor('white');	ax.axis('off')
	cbar_ax=fig.add_axes([0,0.2,1,0.5]);	cbar_ax.axis('off')
	cb=fig.colorbar(ret[1],orientation='horizontal',label='change in drought risk [%]',ax=cbar_ax)
	plt.tight_layout()
	plt.savefig(COU._working_directory+'/plots/spi_risk-1.5_diff_2040s_overview_'+rcp+'.png',)

	###############
	# Tas - similar procedure as for SPI but a bit simpler
	###############

	members_tas=COU.classify_ensemble(['tas',rcp])
	members_tas=members_tas[:-1]
	ensmedian_tas=members_tas[-1]
	COU.model_agreement(ensmedian_tas,members_tas,'diff_ref_2040s','tas_ensmedian_agree_'+rcp)

	ret=COU.plot_map(COU._periods[grid]['diff_ref_2040s'][ensmedian_tas],
				grey_area=COU._periods[grid]['diff_ref_2040s']['tas_ensmedian_agree_'+rcp],
				out_file=COU._working_directory+'/plots/tas_diff_2040s_ensmedian_'+rcp+'.png',
				add_mask=plot_mask,
				color_palette=plt.cm.Reds,
				color_range=[0.5,2],
				color_label='difference in mean \ntemperature [$^\circ C$]',
				title='2040s vs 1986-2005 '+RCP)

	###############
	# Pr - similar procedure as for SPI but a bit simpler
	###############
	members_pr=COU.classify_ensemble(['pr',rcp])
	members_pr=members_pr[:-1]
	ensmedian_pr=members_pr[-1]
	COU.model_agreement(ensmedian_pr,members_pr,'diff_ref_2040s','pr_ensmedian_agree_'+rcp)

	ret=COU.plot_map(COU._periods[grid]['rel_diff_ref_2040s'][ensmedian_pr],
				grey_area=COU._periods[grid]['diff_ref_2040s']['pr_ensmedian_agree_'+rcp],
				out_file=COU._working_directory+'/plots/pr_diff_2040s_ensmedian_'+rcp+'.png',
				add_mask=plot_mask,
				color_palette=plt.cm.RdBu,
				color_range=[-20,20],
				color_label='relative difference \npreciptiation [%]',
				title='2040s vs 1986-2005 '+RCP)

	fig, axes = plt.subplots(nrows=2, ncols=5,figsize=(12,4),subplot_kw={'projection': ccrs.PlateCarree()},gridspec_kw = {'height_ratios':[3,1]})
	for member,ax in zip(members_pr,axes[0,:]):
		ret=COU.plot_map(COU._periods[grid]['rel_diff_ref_2040s'][member],
					ax=ax,
					add_mask=plot_mask,
					color_palette=plt.cm.RdBu,
					color_bar=False,
					color_range=[-10,10],
					title=member.split('_')[1])
	for ax in axes[1,:]:
		ax.outline_patch.set_edgecolor('white');	ax.axis('off')
	cbar_ax=fig.add_axes([0,0.2,1,0.5]);	cbar_ax.axis('off')
	cb=fig.colorbar(ret[1],orientation='horizontal',label='change in precipitation [%]',ax=cbar_ax)
	plt.tight_layout()
	plt.savefig(COU._working_directory+'/plots/pr_diff_2040s_overview_'+rcp+'.png',)












#
