# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.colors
import cartopy.crs as ccrs
import cartopy

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')
try:del sys.modules['country_analysis']
except:pass
import country_analysis; reload(country_analysis)
sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')

os.chdir('/Users/peterpfleiderer/Documents/Projects/isimip_fast-track_processed')



cult_frac=da.read_nc('maize_ha_0.5.nc')['cropdata']
irr_frac=da.read_nc('irrig_maize_0.5.nc')['area'].ix[2,:,:]
rain_frac=cult_frac*(1-irr_frac)

cult_mask=irr_frac.copy()*0.0+1; cult_mask[cult_frac<0.001]=np.nan
irr_mask=irr_frac.copy()*0.0+1; irr_mask[irr_frac<0.001]=np.nan
rain_mask=irr_frac.copy()*0.0+1; rain_mask[rain_frac<0.001]=np.nan

cmap_hist = matplotlib.colors.LinearSegmentedColormap.from_list("", ["saddlebrown","yellow","green"])
cmap_change = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"])

fig,axes=plt.subplots(nrows=3,figsize=(4,6),subplot_kw={'projection': ccrs.PlateCarree()})
for ax,toplot,mask in zip(axes.flatten(),[cult_frac,irr_frac,rain_frac],[cult_mask,irr_mask,rain_mask]):
    ax.pcolormesh(cult_frac.longitude,cult_frac.latitude,toplot*mask,cmap=cmap_hist,vmin=0,vmax=1)
plt.tight_layout()
plt.savefig('masks.png')


all_isos=['DZA','AGO',  'EGY', 'GNQ', 'BEN', 'NGA', 'NER', 'ZWE', 'NAM', 'GNB', 'SWZ', 'GHA', 'COG', 'SLE', 'ETH', 'COM', 'ERI', 'CPV', 'LBR',\
            'LBY', 'LSO', 'UGA', 'RWA', 'SOM', 'MDG', 'CMR', 'TZA', 'BWA', 'SEN', 'TCD', 'GAB', 'BFA', 'MWI', 'MOZ', 'MRT', 'GMB', 'MLI', 'BDI', \
            'STP', 'DJI', 'GIN', 'ESH', 'KEN', 'MAR', 'COD', 'ZMB', 'ZAF', 'TGO', 'TUN', 'CAF', 'SSD', 'SDN', 'CIV','SYC','MUS']


for iso in all_isos:
    print(iso)
    if os.path.isfile('maps/'+'_'.join([iso,'firr','2p0'])+'.pngas')==False:
        # data will be stored in working_directory
        COU=country_analysis.country_analysis(iso,working_directory='cou_data/'+iso+'/',additional_tag='')
        COU.create_mask_country('isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_mai',overwrite=False)
        #COU.create_mask_admin('masked_isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_mai',overwrite=False)

        #lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_yield_mai_2p0.nc4
        for filename in glob.glob('isimip_output_warming/*/*_rcp8p5_ssp2_co2_*yield_mai*'):
            model=filename.split('/')[1]
            gcm=filename.split('/')[-1].split('_')[1]
            rcp=filename.split('/')[-1].split('_')[2]
            ssp=filename.split('/')[-1].split('_')[3]
            co2=filename.split('/')[-1].split('_')[4]
            irr=filename.split('/')[-1].split('_')[5]
            crop=filename.split('/')[-1].split('_')[7]
            wlvl=filename.split('/')[-1].split('_')[8].split('.')[0]

            COU.country_zoom(filename,var_name='yield_mai',given_var_name='yield_mai',data_type='_'.join([ssp,co2,irr,rcp]),scenario=wlvl,model=model+'_'+gcm,time_format='snapshot',overwrite=False)

        COU.ensemble_statistic('median')

        for data in COU.selection(['hist']):
            data.period={'mean':{'year':{}}}
            for wlvl in ['hist','1p0','1p5','2p0','2p5','3p0']:
                data.period['mean']['year'][wlvl]=COU.selection([data.model,data.data_type,wlvl])[0].raw.squeeze()
            COU.period_statistic_diff(data,'mean','year',ref_name='hist')

        COU.period_model_agreement(ref_name='hist',ens_statistic='median')


        for irr,mask,title in zip(['firr','noirr'],[irr_mask,rain_mask],['irrigated','rainfed']):
            ensemble=COU.find_ensemble(['hist','ssp2_co2_'+irr+'_rcp8p5'])
            for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                ensemble['median'].display_map(period='diff_relative_'+wlvl+'-hist',out_file='maps/'+'_'.join([iso,irr,wlvl])+'.png',title=title+' '+wlvl,color_label='relative change in yield [%]',add_mask=mask[ensemble['median'].lat,ensemble['median'].lon].values,color_palette=cmap_change)
            ensemble['median'].display_map(period='hist',out_file='maps/'+'_'.join([iso,irr,'hist'])+'.png',title=title+' historical',color_label='yield [t ha-1 yr-1]',add_mask=mask[ensemble['median'].lat,ensemble['median'].lon].values,color_palette=cmap_hist)


        os.system('rm '+COU._working_directory+'/raw/*')
    # except:
    #     print('***************')



    # COU.period_model_agreement(ref_name='hist')

    # plt.close('all')
    # fig,axes = plt.subplots(nrows=7,ncols=5,figsize=(20,20),subplot_kw={'projection': plate_carree})
    #
    # ensemble['median'].display_map(period='diff_2p0-hist',ax=axes[0,0],color_range=[-0.5,0.5],color_bar=False,title='ensemble median')
    # ensemble['mean'].display_map(period='diff_2p0-hist',ax=axes[0,1],color_range=[-0.5,0.5],color_bar=False,title='ensemble mean')
    # for i,agm in enumerate(['lpjml','lpj-guess','pegasus','pdssat','epic','gepic']):
    #     for j,gcm in enumerate(['gfdl-esm2m','noresm1-m','hadgem2-es','ipsl-cm5a-lr','miroc-esm-chem']):
    #         model=agm+'_'+gcm
    #         tmp=ensemble['models'][model]
    #         tmp.display_map(period='diff_2p0-hist',ax=axes[i+1,j],color_range=[-0.5,0.5],color_bar=False,title=model)
    #
    # #plt.tight_layout()
    # plt.savefig('test__.png')










#
