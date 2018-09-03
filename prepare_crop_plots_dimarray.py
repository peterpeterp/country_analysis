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
cmap_hist.set_under('saddlebrown');cmap_hist.set_over('green')
cmap_change = matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"],N=5)
cmap_change.set_under('red');cmap_change.set_over('green')
cmap_adaptation_potential = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","yellow","green"])
cmap_adaptation_potential.set_over('green')

all_isos=['NGA','MOZ','DZA','AGO',  'EGY', 'GNQ', 'BEN',  'NER', 'ZWE', 'NAM', 'GNB', 'SWZ', 'GHA', 'COG', 'SLE', 'ETH', 'COM', 'ERI', 'CPV', 'LBR',\
            'LBY', 'LSO', 'UGA', 'RWA', 'SOM', 'MDG', 'CMR', 'TZA', 'BWA', 'SEN', 'TCD', 'GAB', 'BFA', 'MWI',  'MRT', 'GMB', 'MLI', 'BDI', \
             'DJI', 'GIN', 'ESH', 'KEN', 'MAR', 'COD', 'ZMB', 'ZAF', 'TGO', 'TUN', 'CAF', 'SSD', 'SDN', 'CIV','SYC','MUS','STP']


for crop,crop_short in zip(['rice','wheat','soybean','maize'][:],['ric','whe','soy','mai'][:]):

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

    cult_mask_01=irr_frac.copy()*0.0+1; cult_mask_01[cult_frac<0.01]=0; cult_mask_01[np.isfinite(cult_frac)==False]=0
    irr_mask_01=irr_frac.copy()*0.0+1; irr_mask_01[irr_frac<0.01]=0; irr_mask_01[cult_frac<0.01]=0; irr_mask_01[np.isfinite(irr_frac)==False]=0
    rain_mask_01=irr_frac.copy()*0.0+1; rain_mask_01[irr_mask_01==1]=0; rain_mask_01[cult_frac<0.01]=0; rain_mask_01[np.isfinite(rain_frac)==False]=0

    mask=cult_mask_01.copy()
    mask[mask==0]=np.nan

    no_cult_mask=cult_mask_01.copy()
    no_cult_mask[no_cult_mask!=0]=np.nan
    no_cult_mask[no_cult_mask==0]=1


    for iso in all_isos:
        print(iso,crop,crop_short)

        if os.path.isfile('checks/'+iso+'_'+crop_short+'_check_ensmedian.png')==False:

            # data will be stored in working_directory
            COU=country_analysis.country_analysis(iso,working_directory='cou_data/'+iso+'/',additional_tag='',time_axis=np.array(['hist','1p0','1p5','2p0','2p5','3p0','3p5','4p0']))
            if os.path.isfile('cou_data/'+iso+'/'+iso+'_'+crop_short+'.nc'):
                COU.load_data(iso+'_'+crop_short+'.nc')
                grid=COU._DATA.keys()[0]

            else:
                COU.create_mask_country('isimip_output_warming/lpjml/lpjml_gfdl-esm2m_rcp8p5_ssp2_co2_firr_aet_mai_1p0.nc4','aet_'+crop_short,overwrite=False)
                grid=COU._DATA.keys()[0]
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
                        COU.country_zoom(filename,var_name='yield_'+crop_short,name='_'.join([crop_short,model,gcm,'noirr']),in_time=np.array([wlvl]),overwrite=False)

                    for filename in glob.glob('isimip_output_warming/*/*_rcp8p5_ssp2_co2_firr_*yield_*'+crop_short+'*'):
                        model=filename.split('/')[1]
                        gcm=filename.split('/')[-1].split('_')[1]
                        rcp=filename.split('/')[-1].split('_')[2]
                        ssp=filename.split('/')[-1].split('_')[3]
                        co2=filename.split('/')[-1].split('_')[4]
                        wlvl=filename.split('/')[-1].split('_')[8].split('.')[0]

                        COU.country_zoom(filename,var_name='yield_'+crop_short,name='_'.join([crop_short,model,gcm,'firr']),in_time=np.array([wlvl]),overwrite=False)
                        COU.country_zoom(filename,var_name='yield_'+crop_short,name='_'.join([crop_short,model,gcm,'total']),in_time=np.array([wlvl]),overwrite=False)

                        data_noirr=COU._DATA[grid]['_'.join([crop_short,model,gcm,'noirr'])]
                        data_firr=COU._DATA[grid]['_'.join([crop_short,model,gcm,'firr'])]

                        COU._DATA[grid]['_'.join([crop_short,model,gcm,'total'])]=data_noirr.copy()*rain_mask_01[data_noirr.lat,data_noirr.lon].values+data_firr.copy()*irr_mask_01[data_noirr.lat,data_noirr.lon].values



            if np.isfinite(np.nanmean(COU._masks['360x720_lat_89.75_-89.75_lon_-179.75_179.75']['lat_weighted'][iso]*mask)):
                for irr in ['firr','noirr','total']:
                    members=COU.classify_ensemble([crop_short,irr])
                    COU.ensemble_statistic(members,'ensmedian_'+'_'.join([crop_short,irr]))

                for wlvl in ['hist','1p0','1p5','2p0','2p5','3p0','3p5']:
                    COU.period_statistic([wlvl,wlvl],wlvl)

                for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                    COU.period_diff('rel_diff_hist_'+wlvl,ref_period='hist',target_period=wlvl,relative=True)

                for wlvl in ['1p0','1p5','2p0','2p5','3p0']:
                    for irr in ['firr','noirr','total']:
                        members=COU.classify_ensemble([crop_short,irr])
                        members=[memb for memb in members if memb != 'ensmedian_'+'_'.join([crop_short,irr])]
                        ensmedian='ensmedian_'+'_'.join([crop_short,irr])
                        COU.model_agreement(ensmedian,members,'rel_diff_hist_'+wlvl,ensmedian+'_agree')
                        COU.small_change(ensmedian,members,'rel_diff_hist_'+wlvl,ensmedian+'_smallChange',ChangeThresh=5)


                cou_mask=np.ma.getdata(COU._masks[grid]['lat_weighted'][iso])

                title=' '.join([iso,crop])
                for wlvl,wlvl_name in zip(['1p0','1p5','2p0','2p5','3p0'],['+1.0$^\circ$C','+1.5$^\circ$C','+2.0$^\circ$C','+2.5$^\circ$C','+3.0$^\circ$C']):
                    ensmedian='ensmedian_'+'_'.join([crop_short,'total'])
                    plot_mask=mask[COU._periods[grid]['rel_diff_hist_'+wlvl].lat,COU._periods[grid]['rel_diff_hist_'+wlvl].lon].values*cou_mask
                    ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid]['rel_diff_hist_'+wlvl][ensmedian],
                    			grey_area=COU._periods[grid]['rel_diff_hist_'+wlvl][ensmedian+'_agree'],
                    			add_mask=plot_mask,
                    			color_palette=cmap_change,
                    			color_range=[-25,25],
                    			color_label='relative change in yield [%]',
                    			title=title+'  '+wlvl_name)

                    smallChange=COU._periods[grid]['rel_diff_hist_'+wlvl][ensmedian+'_smallChange']*plot_mask
                    smallChange[smallChange!=1]=np.nan
                    try:
                        ax.pcolormesh(x,y,smallChange,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow"]),vmin=0,vmax=1)
                    except:
                        pass
                    cbar.ax.get_yaxis().set_ticks([])
                    for j, lab in enumerate(['<-15','-5;-15','5;-5','15;5','>15']):
                        cbar.ax.text(.5, (2 * j + 1) / 10.0, lab, ha='center', va='center',rotation=90,fontsize=8)
                    cbar.ax.get_yaxis().labelpad = 15
                    cbar.ax.set_aspect(10)
                    cbar.ax.set_ylabel('relative change in yield [%]', rotation=90)
                    plt.tight_layout(); plt.savefig('maps/'+'_'.join([iso,crop,'total',wlvl])+'.png')
                    plt.savefig('maps/'+'_'.join([iso,crop,'total',wlvl])+'.pdf')

                    ensmedian_firr='ensmedian_'+'_'.join([crop_short,'firr'])
                    ensmedian_noirr='ensmedian_'+'_'.join([crop_short,'noirr'])
                    to_plot=(COU._periods[grid][wlvl][ensmedian_firr] - COU._periods[grid][wlvl][ensmedian_noirr] ) / COU._periods[grid][wlvl][ensmedian_noirr] *100
                    ax,im,range,x,y,cbar=COU.plot_map(to_plot,
                    			add_mask=plot_mask,
                    			color_palette=cmap_adaptation_potential,
                    			color_range=[0,20],
                    			color_label='irrigation potential [%]',
                    			title=title+'  '+wlvl_name)
                    plt.tight_layout(); plt.savefig('maps/'+'_'.join([iso,crop,'irr-added-value',wlvl])+'.png')
                    plt.savefig('maps/'+'_'.join([iso,crop,'irr-added-value',wlvl])+'.pdf')

                ensmedian='ensmedian_'+'_'.join([crop_short,'total'])
                ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid]['hist'][ensmedian],
                			add_mask=plot_mask,
                			color_palette=cmap_hist,
                            color_range=np.nanpercentile(COU._periods[grid]['hist'][ensmedian],[10,90]),
                			color_label='yield [t/ha/yr]',
                			title=title+'  historical')
                plt.tight_layout(); plt.savefig('maps/'+'_'.join([iso,crop,'total','hist'])+'.png')
                plt.savefig('maps/'+'_'.join([iso,crop,'total','hist'])+'.pdf')

                ensmedian='ensmedian_'+'_'.join([crop_short,'total'])
                ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid]['hist'][ensmedian],
                			add_mask=plot_mask,
                			color_palette=cmap_hist,
                            color_range=np.nanpercentile(COU._periods[grid]['hist'][ensmedian],[10,90]),
                			color_label='yield [t/ha/yr]',
                			title=title+'  historical')
                plt.tight_layout(); plt.savefig('maps/'+'_'.join([iso,crop,'total','hist'])+'.png')
                plt.savefig('maps/'+'_'.join([iso,crop,'total','hist'])+'.pdf')

                COU.save_data(iso+'_'+crop_short+'.nc')

                if True:
                    fig,axes = plt.subplots(nrows=6,ncols=6 ,figsize=(20,20),subplot_kw={'projection': ccrs.PlateCarree()})
                    for wlvl,row in zip(['hist','1p0','1p5','2p0','2p5','3p0'],np.arange(0,6)):
                        ensmedian_total='ensmedian_'+'_'.join([crop_short,'total'])
                        to_plot=(COU._periods[grid][wlvl][ensmedian_total]-COU._periods[grid]['hist'][ensmedian_total]) / COU._periods[grid]['hist'][ensmedian_total] *100
                        ax,im,range,x,y,cbar=COU.plot_map(to_plot,
                        			add_mask=plot_mask,
                        			color_palette=cmap_change,
                        			color_range=[-25,25],
                        			color_label='rel change yield',
                        			title='total diff  '+wlvl,
                                    ax=axes[row,0])
                        for irr,ax in zip(['total','noirr','firr'],axes[row,1:4]):
                            ensmedian='ensmedian_'+'_'.join([crop_short,irr])
                            ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid][wlvl][ensmedian],
                            			add_mask=plot_mask,
                            			color_palette=cmap_hist,
                            			color_range=[1,2.2],
                            			color_label='yield',
                            			title=irr+'  '+wlvl,
                                        ax=ax)
                        ensmedian_firr='ensmedian_'+'_'.join([crop_short,'firr'])
                        ensmedian_noirr='ensmedian_'+'_'.join([crop_short,'noirr'])
                        ensmedian_total='ensmedian_'+'_'.join([crop_short,'total'])
                        ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid][wlvl][ensmedian_firr]-COU._periods[grid][wlvl][ensmedian_noirr],
                        			add_mask=plot_mask,
                        			color_palette=cmap_adaptation_potential,
                        			color_range=[0,1],
                        			color_label='yield',
                        			title='firr - noirr '+wlvl,
                                    ax=axes[row,4])
                        ax,im,range,x,y,cbar=COU.plot_map(COU._periods[grid][wlvl][ensmedian_total]-COU._periods[grid][wlvl][ensmedian_noirr],
                        			add_mask=plot_mask,
                        			color_palette=cmap_adaptation_potential,
                        			color_range=[0,1],
                        			color_label='yield',
                        			title='total - noirr '+wlvl,
                                    ax=axes[row,5])
                    plt.tight_layout(); plt.savefig('checks/'+iso+'_'+crop_short+'_check_ensmedian.png')


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
