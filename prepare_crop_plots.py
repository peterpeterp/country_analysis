# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd

sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')
try:del sys.modules['country_analysis']
except:pass
import country_analysis; reload(country_analysis)
sys.path.append('/Users/peterpfleiderer/Documents/Projects/country_analysis/')

os.chdir('/Users/peterpfleiderer/Documents/Projects/isimip_fast-track_processed')
#
# #mask agricultural area
# for filename in glob.glob('isimip_output_warming/*/*'):
#     model=filename.split('/')[1]
#     if os.path.isdir('masked_isimip_output_warming/'+model)==False:
#         os.system('mkdir masked_isimip_output_warming/'+model)
#     if filename.split('_')[-2]=='mai':
#         os.system('cdo mul '+filename+' mask_maize_area.nc'+' '+filename.replace('isimip_output_warming','masked_isimip_output_warming'))
#

all_isos=['AGO', 'DZA', 'EGY', 'GNQ', 'BEN', 'NGA', 'NER', 'ZWE', 'NAM', 'GNB', 'SWZ', 'GHA', 'COG', 'SLE', 'ETH', 'COM', 'ERI', 'CPV', 'LBR',\
            'LBY', 'LSO', 'UGA', 'RWA', 'SOM', 'MDG', 'CMR', 'TZA', 'BWA', 'SEN', 'TCD', 'GAB', 'BFA', 'MWI', 'MOZ', 'MRT', 'GMB', 'MLI', 'BDI', \
            'STP', 'DJI', 'GIN', 'ESH', 'KEN', 'MAR', 'COD', 'ZMB', 'ZAF', 'TGO', 'TUN', 'CAF', 'SSD', 'SDN', 'CIV','SYC','MUS']

all_isos=['GNQ']

for iso in all_isos:
    print(iso)
    if os.path.isfile('maps/'+'_'.join([iso,'firr','2p0'])+'.png_saads')==False:
        # data will be stored in working_directory
        COU=country_analysis.country_analysis(iso,working_directory='cou_data/'+iso+'/',additional_tag='')
        COU.create_mask_country('masked_isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_mai',overwrite=False)
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

        for irr in ['firr','noirr']:
            ensemble=COU.find_ensemble(['hist','ssp2_co2_'+irr+'_rcp8p5'])
            for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                ensemble['median'].display_map(period='diff_'+wlvl+'-hist',color_range=[-0.25,0.25],out_file='maps/'+'_'.join([iso,irr,wlvl])+'.png',title=irr+' '+wlvl,color_label='change in yield')

        #os.system('rm '+COU._working_directory+'/raw/*')
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
