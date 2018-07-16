# -*- coding: utf-8 -*-

import sys,glob,os,itertools,datetime,pickle,subprocess,time
import numpy as np
from netCDF4 import Dataset,num2date
import pandas as pd
from shapely.geometry import mapping, Polygon, MultiPolygon
import matplotlib.pylab as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from unidecode import unidecode
import dimarray as da

sys.path.append('/Users/peterpfleiderer/Projects/country_analysis/')
try:del sys.modules['country_analysis']
except:pass
import country_analysis; reload(country_analysis)
sys.path.append('/Users/peterpfleiderer/Projects/country_analysis/')

os.chdir('/Users/peterpfleiderer/Projects/isimip_fast-track_processed')

cmap_hist = matplotlib.colors.LinearSegmentedColormap.from_list("", ["saddlebrown","yellow","green"])
cmap_change = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"])
cmap_adaptation_potential = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","yellow","green"])

all_isos=['NGA','MOZ','DZA','AGO',  'EGY', 'GNQ', 'BEN',  'NER', 'ZWE', 'NAM', 'GNB', 'SWZ', 'GHA', 'COG', 'SLE', 'ETH', 'COM', 'ERI', 'CPV', 'LBR',\
            'LBY', 'LSO', 'UGA', 'RWA', 'SOM', 'MDG', 'CMR', 'TZA', 'BWA', 'SEN', 'TCD', 'GAB', 'BFA', 'MWI',  'MRT', 'GMB', 'MLI', 'BDI', \
            'STP', 'DJI', 'GIN', 'ESH', 'KEN', 'MAR', 'COD', 'ZMB', 'ZAF', 'TGO', 'TUN', 'CAF', 'SSD', 'SDN', 'CIV','SYC','MUS']


for crop,crop_short in zip(['rice','wheat','soybean','maize'],['ric','whe','soy','mai']):

    try:
        cult_frac=da.read_nc('masks/'+crop+'_ha_0.5.nc')['cropdata']
    except:
        cult_frac=da.read_nc('masks/'+crop+'_0.5.nc')['surta'].ix[0,0,:,:].squeeze()

    ds=da.Dataset({'area':cult_frac})
    ds.write_nc('masks/'+crop+'_cultivated_frac.nc')

    irr_frac=da.read_nc('masks/irrig_'+crop+'_0.5.nc')['area'].ix[2,:,:]
    ds=da.Dataset({'area':irr_frac})
    ds.write_nc('masks/'+crop+'_irrigated_frac.nc')

    rain_frac=cult_frac*(1-irr_frac)
    ds=da.Dataset({'area':rain_frac})
    ds.write_nc('masks/'+crop+'_rainfed_frac.nc')

    cult_mask_01=irr_frac.copy()*0.0+1; cult_mask_01[cult_frac<0.001]=0; cult_mask_01[np.isfinite(cult_frac)==False]=0
    irr_mask_01=irr_frac.copy()*0.0+1; irr_mask_01[irr_frac<0.001]=0; irr_mask_01[np.isfinite(irr_frac)==False]=0
    rain_mask_01=irr_frac.copy()*0.0+1; rain_mask_01[irr_mask_01==1]=0; rain_mask_01[cult_frac<0.001]=0; rain_mask_01[np.isfinite(rain_frac)==False]=0

    mask=cult_mask_01.copy()
    mask[mask==0]=np.nan


    for iso in all_isos:
        print(iso,crop,crop_short)
        if os.path.isfile('maps/'+'_'.join([iso,'firr','2p0'])+'.png')==False or True:
            # data will be stored in working_directory
            COU=country_analysis.country_analysis(iso,working_directory='cou_data/'+iso+'/',additional_tag='')
            COU.create_mask_country('isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_mai',overwrite=False)
            #COU.create_mask_admin('masked_isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_mai',overwrite=False)

            print(np.nanmean(COU._masks['360x720_lat_89.75_-89.75_lon_-179.75_179.75']['lat_weighted'][iso]*mask))
            if np.isfinite(np.nanmean(COU._masks['360x720_lat_89.75_-89.75_lon_-179.75_179.75']['lat_weighted'][iso]*mask)):

                #lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_yield__+crop_short2p0.nc4
                for filename in glob.glob('isimip_output_warming/*/*_rcp8p5_ssp2_co2_noirr_*yield_*'+crop_short+'*'):
                    model=filename.split('/')[1]
                    gcm=filename.split('/')[-1].split('_')[1]
                    rcp=filename.split('/')[-1].split('_')[2]
                    ssp=filename.split('/')[-1].split('_')[3]
                    co2=filename.split('/')[-1].split('_')[4]
                    wlvl=filename.split('/')[-1].split('_')[8].split('.')[0]
                    COU.country_zoom(filename,var_name='yield_'+crop_short,given_var_name='yield_'+crop_short,data_type='_'.join([ssp,co2,'noirr',rcp]),scenario=wlvl,model=model+'_'+gcm,time_format='snapshot',overwrite=False)

                for filename in glob.glob('isimip_output_warming/*/*_rcp8p5_ssp2_co2_firr_*yield_*'+crop_short+'*'):
                    model=filename.split('/')[1]
                    gcm=filename.split('/')[-1].split('_')[1]
                    rcp=filename.split('/')[-1].split('_')[2]
                    ssp=filename.split('/')[-1].split('_')[3]
                    co2=filename.split('/')[-1].split('_')[4]
                    wlvl=filename.split('/')[-1].split('_')[8].split('.')[0]

                    COU.country_zoom(filename,var_name='yield_'+crop_short,given_var_name='yield_'+crop_short,data_type='_'.join([ssp,co2,'firr',rcp]),scenario=wlvl,model=model+'_'+gcm,time_format='snapshot',overwrite=False)
                    COU.country_zoom(filename,var_name='yield_'+crop_short,given_var_name='yield_'+crop_short,data_type='_'.join([ssp,co2,'total',rcp]),scenario=wlvl,model=model+'_'+gcm,time_format='snapshot',overwrite=False)
                    data_noirr=COU.selection([wlvl,model+'_'+gcm,'_'.join([ssp,co2,'noirr',rcp])])[0]
                    data_firr=COU.selection([wlvl,model+'_'+gcm,'_'.join([ssp,co2,'firr',rcp])])[0]
                    data_total=COU.selection([wlvl,model+'_'+gcm,'_'.join([ssp,co2,'total',rcp])])[0]
                    data_total.raw=data_noirr.raw.copy()*rain_mask_01[data_noirr.lat,data_noirr.lon].values+data_firr.raw.copy()*irr_mask_01[data_noirr.lat,data_noirr.lon].values


                    # ds=da.Dataset({'yield':da.DimArray(data_noirr.raw.squeeze()*rain_mask_01[data_noirr.lat,data_noirr.lon].values,axes=[data_noirr.lat,data_noirr.lon],dims=['lat','lon'])}).write_nc('noirr.nc')
                    # ds=da.Dataset({'yield':da.DimArray(data_firr.raw.squeeze()*irr_mask_01[data_noirr.lat,data_noirr.lon].values,axes=[data_firr.lat,data_firr.lon],dims=['lat','lon'])}).write_nc('firr.nc')
                    # ds=da.Dataset({'yield':da.DimArray(data_total.raw.squeeze(),axes=[data_total.lat,data_total.lon],dims=['lat','lon'])}).write_nc('total.nc')


                COU.ensemble_statistic('median')

                for data in COU.selection(['hist','_'.join([ssp,co2,'total',rcp])]):
                    data.period={'mean':{'year':{}}}
                    for wlvl in ['hist','1p0','1p5','2p0','2p5','3p0']:
                        data.period['mean']['year'][wlvl]=COU.selection([data.model,data.data_type,wlvl])[0].raw.squeeze()
                    COU.period_statistic_diff(data,'mean','year',ref_name='hist')
                    for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                        data.period['mean']['year'][wlvl+'_firr']=COU.selection([data.model,data.data_type.replace('total','firr'),wlvl])[0].raw.squeeze()
                        data.period['mean']['year'][wlvl+'_noirr']=COU.selection([data.model,data.data_type.replace('total','noirr'),wlvl])[0].raw.squeeze()
                        COU.period_statistic_diff(data,'mean','year',ref_name=wlvl+'_noirr',proj_period_names=[wlvl+'_firr'])

                for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                    COU.period_model_agreement(ref_name=wlvl+'_noirr',ens_statistic='median',proj_period_names=[wlvl+'_firr'],relChangeThresh=5)

                COU.period_model_agreement(ref_name='hist',ens_statistic='median',relChangeThresh=5)



                for co2 in ['co2']:
                    ens_mean=COU.find_ensemble(['hist','ssp2_'+co2+'_total_rcp8p5'])['median']
                    title=' '.join([iso,crop])
                    for wlvl,wlvl_name in zip(['1p0','1p5','2p0','2p5','3p0'],['+1.0$^\circ$C','+1.5$^\circ$C','+2.0$^\circ$C','+2.5$^\circ$C','+3.0$^\circ$C']):
                        ens_mean.display_map(period='diff_relative_'+wlvl+'-hist',
                                                out_file='maps/'+'_'.join([iso,crop,'total',wlvl])+'.png',
                                                title=title+'  '+wlvl_name,
                                                color_label='relative change in yield [%]',
                                                add_mask=mask[ens_mean.lat,ens_mean.lon].values,
                                                color_palette=cmap_change,
                                                color_range=[-20,20])
                        ens_mean.display_map(period='diff_relative_'+wlvl+'_firr-'+wlvl+'_noirr',
                                                out_file='maps/'+'_'.join([iso,crop,'irr-added-value',wlvl])+'.png',
                                                title=title+'  '+wlvl_name,
                                                color_label='irrigation potential [%]',
                                                add_mask=mask[ens_mean.lat,ens_mean.lon].values,
                                                color_palette=cmap_adaptation_potential,
                                                color_range=[0,20])
                    ens_mean.display_map(period='hist',
                                out_file='maps/'+'_'.join([iso,crop,'total','hist'])+'.png',
                                title=title+'  historical',
                                color_label='yield [t/ha/yr]',
                                add_mask=mask[ens_mean.lat,ens_mean.lon].values,
                                color_palette=cmap_hist)

                os.system('rm '+COU._working_directory+'/raw/*')


# cult_mask=irr_frac.copy()*0.0+1; cult_mask[cult_frac<0.001]=np.nan
# irr_mask=irr_frac.copy()*0.0+1; irr_mask[irr_frac<0.001]=np.nan
# rain_mask=irr_frac.copy()*0.0+1; rain_mask[rain_frac<0.001]=np.nan
#
# masks={'noirr':rain_mask,'firr':irr_mask}
#

#
# fig,axes=plt.subplots(nrows=3,figsize=(4,6),subplot_kw={'projection': ccrs.PlateCarree()})
# for ax,toplot,mask in zip(axes.flatten(),[cult_frac,irr_frac,rain_frac],[cult_mask,irr_mask,rain_mask]):
#     ax.pcolormesh(cult_frac.longitude,cult_frac.latitude,toplot*mask,cmap=cmap_hist,vmin=0,vmax=1)
# plt.tight_layout()
# plt.savefig('masks.png')
