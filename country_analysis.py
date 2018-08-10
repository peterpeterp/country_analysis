# -*- coding: utf-8 -*-
'''
Class to analyze climate data on national (& sub national) scale
Peter Pfleiderer
peter.pfleiderer@climateanalytics.org
'''

import sys,glob,os,itertools,datetime,pickle,subprocess,time
import numpy as np
from netCDF4 import Dataset,num2date
import pandas as pd
from shapely.geometry import mapping, Polygon, MultiPolygon
import matplotlib.pylab as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
from unidecode import unidecode
import dimarray as da

'''
more elegant with subprocess?
'''


from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "sans-serif"
plt.style.use('classic')

def depth(d, level=1):
    if not isinstance(d, dict) or not d:
        return level
    return max(depth(d[k], level + 1) for k in d)

def running_mean_func(xx,N):
    if N==1:
        return xx
    if N!=1:
        x=np.ma.masked_invalid(xx.copy())
        ru_mean=x.copy()*np.nan
        for t in range(int(N/2),len(x)-int(N/2)):
            ru_mean[t]=np.nanmean(x[t-int(N/2):t+int(N/2)])
        return ru_mean


class country_analysis(object):

    def __init__(self,iso,working_directory,seasons={'year':range(1,13)},additional_tag=''):
        '''
        Prepare directories and meta-data
        iso: str: Country name or iso
        working_directory: path: Path where files will be stored
        seasons: dict: seasons relevant for the country. 'season name':{months in season as int 1-12}
        additional_tag: str: when specified raw data will be stored in a separate directory with additional tag name.
        '''

        self._iso=iso
        #self._working_directory=working_directory
        self._working_directory=os.getcwd()+'/'+working_directory

        self._additional_tag=additional_tag
        self._working_directory_raw=self._working_directory+'/raw'+additional_tag

        self._seasons=seasons

        self._masks={}
        self._grid_dict={}
        self._DATA={}
        self._periods={}
        self._region_names={}

    	year=np.array([[yr]*12 for yr in range(1900,2101)]).flatten()
    	month=np.array([[mon]*len(range(1900,2101)) for mon in range(12)]).T.flatten()
        self._time_axis=np.array([yr+float(mn)/12. for yr,mn in zip(year,month)])

        if os.path.isdir(self._working_directory)==False:os.system('mkdir '+self._working_directory)
        if os.path.isdir(self._working_directory+'/masks')==False:os.system('mkdir '+self._working_directory+'/masks')
        if os.path.isdir(self._working_directory+'/plots')==False:os.system('mkdir '+self._working_directory+'/plots')
        if os.path.isdir(self._working_directory_raw)==False:os.system('mkdir '+self._working_directory_raw)
        if os.path.isdir(self._working_directory+'/area_average')==False:os.system('mkdir '+self._working_directory+'/area_average')

        # get shapefiles for country
        if os.path.isdir(self._working_directory+'/'+iso+'_adm_shp')==False:
            current_dir=os.getcwd()
            os.chdir(self._working_directory)
            os.system('wget biogeo.ucdavis.edu/data/gadm2.8/shp/'+iso+'_adm_shp.zip')
            os.system('mkdir '+iso+'_adm_shp')
            os.system('ls')
            os.chdir(self._working_directory+iso+'_adm_shp')
            os.system('unzip ../'+iso+'_adm_shp.zip')
            os.chdir(current_dir)

        # load required shapefiles
        print 'load regions'
        start_time=time.time()
        adm_shapefiles=shapereader.Reader(self._working_directory+self._iso+'_adm_shp/'+self._iso+'_adm1').records()

        # collect all shapes of region
        self._adm_polygons={}
        for item in adm_shapefiles:
            shape,region=item.geometry,item.attributes
            region = {k.lower():v for k,v in region.items()}
            name_full = region['name_1']
            name=unidecode(name_full.decode('utf8')).replace(' ','_')
            self._region_names[name]=name_full
            # simplify could be added here to speed up things
            self._adm_polygons[name]=MultiPolygon(shape)

        # for region_name in self._region_names.keys():
        # 	if '+' in region_name:
        # 		sub_regs=region_name.split('+')
        # 		self._adm_polygons[region_name]=self._adm_polygons[sub_regs[0]]
        # 		for region in sub_regs[1:]:
        # 			self._adm_polygons[region_name] = \
        # 			self._adm_polygons[region_name].symmetric_difference(self._adm_polygons[region])

        adm_shapefiles=shapereader.Reader(self._working_directory+self._iso+'_adm_shp/'+self._iso+'_adm0').records()
        name=self._iso
        self._region_names[name]=name
        self._adm_polygons[self._iso]=MultiPolygon(next(adm_shapefiles).geometry)

        print self._adm_polygons.keys()
        print 'regions loaded '+str(time.time()-start_time)

    def save_data(self):
        da.Dataset(self._DATA).write_nc(self._working_directory+'/data.nc')
        os.system('rm '+self._working_directory+'/raw/*')

    def load_data(self):
        tmp=da.read_nc(self._working_directory+'/data.nc')
        self._DATA={}
        for grid in tmp.keys():
            self._DATA[grid]=tmp[grid]


        for file in glob.glob(self._working_directory+'/masks/'+self._iso+'*.nc*'):
            # need the country mask first
            if 'admin' not in file.split('_') and len(file.split('+'))==1:
                file_new=self._working_directory+'/masks'+file.split('masks')[-1]
                self.load_masks(file_new)

        for file in glob.glob(self._working_directory+'/masks/'+self._iso+'*.nc*'):
            if 'admin' in file.split('_'):
                if len(file.split('+'))==1 or load_merged_regions: # exclude merged regions
                    file_new=self._working_directory+'/masks'+file.split('masks')[-1]
                    self.load_masks(file_new)

    ###########
    # masks
    ###########

    def identify_grid(self,input_file,lat_name,lon_name):
        '''
        get information about grid of input data
        input_file: file_path: file to be analyzed
        lat_name: str: name of latitude variable
        lon_name: str: name of longitude variable
        '''
        nc_in=Dataset(input_file,'r')
        lat = nc_in.variables[lat_name][:]
        lon = nc_in.variables[lon_name][:].squeeze()

        if len(lat.shape)==2:
            lat=lat[:,0]
            lon=lon[0,:]

        # formerly shift_lon
        if max(lon)>200:	lon_shift=-180.0
        else:				lon_shift=0.0
        lon+=lon_shift

        nx = len(lon)	;	ny = len(lat)
        grid=str(ny)+'x'+str(nx)+'_lat_'+str(lat[0])+'_'+str(lat[-1])+'_lon_'+str(lon[0])+'_'+str(lon[-1])
        nc_in.close()

        return lon,lat,grid,lon_shift

    def load_masks(self,mask_file):
        '''
        load existing mask
        mask_file: file_path: mask file to load
        '''
        nc_mask=Dataset(mask_file,'r')

        grid=nc_mask.getncattr('original_grid')
        mask_style=nc_mask.getncattr('mask_style')

        if grid not in self._masks.keys():
            self._masks[grid]={}
        if mask_style not in self._masks[grid].keys():
            self._masks[grid][mask_style]={}

        self._masks[grid]['lat_mask'] = nc_mask.variables['lat'][:]
        self._masks[grid]['lon_mask'] = nc_mask.variables['lon'][:]

        # get all variables (regions)
        for name in nc_mask.variables.keys():
            if name not in ['lat','lon']:
                self._masks[grid][mask_style][name] = nc_mask.variables[name][:,:]
                self.zoom_mask(grid,mask_style,name)

    def get_grid_polygons(self,grid,lon,lat,lon_shift):
        '''
        create polygons for each grid-cell
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        lon_shift: float: deg longitudes that have to be shifted to be on a -180 to 180 grid (computed in identify_grid)
        '''
        # loop over the grid to get grid polygons
        nx = len(lon)	;	ny = len(lat)

        grid_polygons = np.empty((nx,ny),dtype=Polygon)
        dx = np.zeros((nx))
        dy = np.zeros((ny))
        dx[1:] = np.abs(np.diff(lon,1))
        dx[0] = dx[1]
        dy[1:] = np.abs(np.diff(lat,1))
        dy[0] = dy[1]
        for i in range(nx):
            x1 = lon[i]-dx[i]/2.
            x2 = lon[i]+dx[i]/2.
            for j in range(ny):
                y1 = lat[j]-dy[j]/2.
                y2 = lat[j]+dy[j]/2.
                grid_polygons[i,j] = Polygon([(x1,y1),(x1,y2),(x2,y2),(x2,y1)])
                #grid_polygons[i,j] = Polygon([(y1,x1),(y1,x2),(y2,x2),(y2,x1)])

        # since the lon axis has been shifted, masks and outputs will have to be shifted as well. This shift is computed here
        lon-=lon_shift
        shift = len(lon)-np.where(lon==lon[0]-lon_shift)[0][0]

        self._masks[grid]['lat_mask']=lat
        self._masks[grid]['lon_mask']=lon

        return grid_polygons,shift

    def merge_adm_regions(self,region_names,new_region_name=None):
        if new_region_name is None:
            single_regions=[]
            for region in region_names:
                for split in region.split('+'):
                    single_regions.append(split)
            new_region_name='+'.join(sorted(single_regions))
            self._region_names[new_region_name]='+'.join([self._region_names[reg] for reg in sorted(single_regions)])
        self._adm_polygons[new_region_name]=self._adm_polygons[region_names[0]]
        for region in region_names[1:]:
            self._adm_polygons[new_region_name] = \
            self._adm_polygons[new_region_name].symmetric_difference(self._adm_polygons[region])
        return new_region_name

    def get_region_area(self,region):
        poly=self._adm_polygons[region]
        lat=poly.centroid.xy[1][0]
        return({'km2':poly.area*(12742./360.)**2*np.cos(np.radians(lat))*10,'latxlon':poly.area})

    def regrid_pop_mask(self,grid,lon,lat,shift,pop_mask_file,mask_style):
        '''
        regrid population masks
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        shift: int: number of elements to roll around longitude axis. Has to be considered since the longitude axis might be shifted (see get_grid_polygons)
        pop_mask_file: file_path: path of used population mask
        mask_style: str: name of the created mask
        '''
        mygrid=open(self._working_directory+'/masks/'+grid+'.txt','w')
        mygrid.write('gridtype=lonlat\nxsize='+str(len(lon))+'\nysize='+str(len(lat))+'\nxfirst='+str(lon[0])+'\nxinc='+str(np.mean(np.diff(lon,1)))+'\nyfirst='+str(lat[0])+'\nyinc='+str(np.mean(np.diff(lat,1))))
        mygrid.close()
        os.system('cdo remapbil,'+self._working_directory+'/masks/'+grid+'.txt '+pop_mask_file+' '+self._working_directory+'/masks/'+mask_style+'_'+grid+'.nc')
        nc_pop_mask = Dataset(self._working_directory+'/masks/'+mask_style+'_'+grid+'.nc')
        pop_mask = np.array(nc_pop_mask.variables['mask'][:,:]).squeeze()
        pop_mask = np.roll(pop_mask,shift,axis=1)

        return pop_mask

    def grid_polygon_overlap(self,grid,lon,lat,grid_polygons,country_polygons,shift,mask_style,ext_poly,name,pop_mask=None):
        '''
        Compute overlap betwwen grid polygons (get_grid_polygons) and country polygons
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        grid_polygons: list: List of polygons created in get_grid_polygons
        country_polgons: list: list of polygons representing the country
        shift: int: number of elements to roll around longitude axis. Has to be considered since the longitude axis might be shifted (see get_grid_polygons)
        mask_style: str: Can be 'lat_weighted' or population weighted. If population weighted, mask_style is a given name
        est_poly: Polygon: Polygon limiting the are where overlaps are computed
        name: str: country-name or region name
        pop_mask: np.array: population mask from regrid_pop_mask
        '''
        nx = len(lon)	;	ny = len(lat)

        overlap = np.zeros((ny,nx))
        for i in range(nx):
            for j in range(ny):
                # check gridcell is relevant
                if grid_polygons[i,j].intersects(ext_poly):
                    # get fraction of grid-cell covered by polygon
                    intersect = grid_polygons[i,j].intersection(country_polygons).area/grid_polygons[i,j].area*country_polygons.area
                    if pop_mask is not None:
                        # population weighting
                        overlap[j,i] = intersect*pop_mask[j,i]
                    if mask_style=='lat_weighted':
                        # multiply overlap with latitude weighting
                        overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

        # renormalize overlap to get sum(mask)=1
        overlap_sum=sum(overlap.copy().flatten())
        if overlap_sum!=0:
            output=overlap.copy()/overlap_sum
            # mask zeros
            output[output==0]=np.nan
            output=np.ma.masked_invalid(output)
            # shift back to original longitudes
            self._masks[grid][mask_style][name]=np.roll(output,shift,axis=1)
            return True
        else:
            print 'something went wrong with the mask'
            return False

    def create_mask_country(self,input_file,var_name,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon'):
        '''
        create country mask
        input_file: str: location of example input data (required for the identification of the grid)
        var_name: str: variable name of input file
        shape_file: str: location of the shape_file used to identify country borders
        mask_style: str: name under which the mask will be stored (important for further analysis)
        pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
        overwrite: bool: if True, old files is deleted, new mask created
        lat_name: str: name of latitude variable in netcdf file
        lon_name: str: name of longitude variable in netcdf file
        '''

        lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

        if grid not in self._masks.keys():
            self._masks[grid]={}
        if mask_style not in self._masks[grid].keys():
            self._masks[grid][mask_style]={}

        mask_file=self._working_directory+'/masks/'+self._iso+'_'+grid+'_'+mask_style+'.nc4'

        if os.path.isfile(mask_file) and overwrite==False:
                self.load_masks(mask_file)

        if os.path.isfile(mask_file) and overwrite:
            os.system('rm '+mask_file)

        if os.path.isfile(mask_file)==False:
            grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)

            country_polygons = self._adm_polygons[self._iso]

            # get boundaries for faster computation
            x1, y1, x2, y2 = country_polygons.bounds
            xmin, xmax, ymin, ymax = min([x1,x2]), max([x1,x2]), min([y1,y2]), max([y1,y2])
            ext = [(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin),(xmin,ymin)]
            ext_poly = Polygon(ext)

            # load population mask
            if pop_mask_file=='':
                pop_mask = np.ones((len(lat),len(lon)))
            else:
                pop_mask = self.regrid_pop_mask(grid,lon,lat,shift,pop_mask_file,mask_style)

            # compute overlap
            self.grid_polygon_overlap(grid,lon, lat, grid_polygons, country_polygons, shift, mask_style, ext_poly, self._iso, pop_mask)

            # save mask
            print mask_file
            nc_mask=Dataset(mask_file,'w')
            nc_mask.createDimension('lat', len(lat))
            nc_mask.createDimension('lon', len(lon))
            outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
            outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')
            outVar = nc_mask.createVariable(self._iso, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][self._iso][:,:]

            nc_mask.setncattr('original_grid',grid)
            nc_mask.setncattr('mask_style',mask_style)
            nc_mask.close()

            self.zoom_mask(grid,mask_style,self._iso)

    def create_mask_admin(self,input_file,var_name,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon',regions=None):
        '''
        create country mask
        input_file: str: location of example input data (required for the identification of the grid)
        var_name: str: variable name of input file
        shape_file: str: location of the shape_file used to identify country borders
        mask_style: str: name under which the mask will be stored (important for further analysis)
        pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
        overwrite: bool: if True, old files is deleted, new mask created
        lat_name: str: name of latitude variable in netcdf file
        lon_name: str: name of longitude variable in netcdf file
        '''

        lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

        if grid not in self._masks.keys():
            self._masks[grid]={}
        if mask_style not in self._masks[grid].keys():
            self._masks[grid][mask_style]={}

        if regions is None:
            mask_file=self._working_directory+'/masks/'+self._iso+'_admin_'+grid+'_'+mask_style+'_all.nc4'
        if regions is not None:
            mask_file=self._working_directory+'/masks/'+self._iso+'_admin_'+grid+'_'+mask_style+'_'+'_'.join(regions)+'.nc4'

        if os.path.isfile(mask_file) and overwrite==False:
            self.load_masks(mask_file)

        if os.path.isfile(mask_file) and overwrite:
            os.system('rm '+mask_file)

        if os.path.isfile(mask_file)==False:
            grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)

            region_polygons=self._adm_polygons

            # get boundaries for faster computation
            xs, ys = [], []
            for name in region_polygons.keys():
                bounds=region_polygons[name].bounds
                xs.append(bounds[0])
                xs.append(bounds[2])
                ys.append(bounds[1])
                ys.append(bounds[3])
            xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
            ext = [(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin),(xmin,ymin)]
            ext_poly = Polygon(ext)

            # load population mask
            if pop_mask_file=='':
                pop_mask = np.ones((len(lat),len(lon)))
            else:
                pop_mask = self.regrid_pop_mask(grid,lon,lat,shift,pop_mask_file,mask_style)

            # prepare outputfile
            nc_mask=Dataset(mask_file,'w')
            nc_mask.createDimension('lat', len(lat))
            nc_mask.createDimension('lon', len(lon))
            outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
            outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')

            if regions is None:
                selected_regions=region_polygons.keys()
            if regions is not None:
                selected_regions=regions

            for name in selected_regions:
                #print name,region_polygons.keys()
                region_shape = region_polygons[name]
                if self.grid_polygon_overlap(grid,lon, lat, grid_polygons, region_shape, shift, mask_style, ext_poly, name, pop_mask):
                    outVar = nc_mask.createVariable(name, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][name][:,:]
                self.zoom_mask(grid,mask_style,name)

            nc_mask.setncattr('original_grid',grid)
            nc_mask.setncattr('mask_style',mask_style)
            nc_mask.close()

    def zoom_mask(self,grid,mask_style,region):
        '''
        store a mask file restricted to the country relevant rectangle
        grid: str: name of the grid
        mask_style: str: Can be 'lat_weighted' or population weighted. If population weighted, mask_style is a given name
        region: str: country-name or region name
        '''
        mask=self._masks[grid][mask_style][region]
        lat_mask=self._masks[grid]['lat_mask']
        lon_mask=self._masks[grid]['lon_mask']

        cou_mask=self._masks[grid][mask_style][self._iso]
        cou_mask=np.ma.getdata(cou_mask)

        lon_mean=np.nanmean(cou_mask,0)
        #lons=np.where(lon_mean!=0)[0]
        lons=sorted(np.where(np.isfinite(lon_mean))[0])
        lon_=lon_mask[lons[0]:lons[-1]+1]

        lat_mean=np.nanmean(cou_mask,1)
        #lats=np.where(lat_mean!=0)[0]
        lats=sorted(np.where(np.isfinite(lat_mean))[0])
        lat_=lat_mask[lats[0]:lats[-1]+1]

        small_grid=str(len(lat_))+'x'+str(len(lon_))+'_lat_'+str(lat_[0])+'_'+str(lat_[-1])+'_lon_'+str(lon_[0])+'_'+str(lon_[-1])
        if small_grid not in self._masks.keys():	self._masks[small_grid]={}
        if mask_style not in self._masks[small_grid].keys():	self._masks[small_grid][mask_style]={}

        print(small_grid,len(lons),len(lats))
        self._masks[small_grid][mask_style][region]=mask[lats[0]:lats[-1]+1,lons[0]:lons[-1]+1]
        self._grid_dict[grid]=small_grid

        if small_grid not in self._DATA.keys():
            self._DATA[small_grid] = da.DimArray(axes=[['dummy'],self._time_axis,np.ma.getdata(lat_),np.ma.getdata(lon_)],dims=['name','time','lat','lon'])

        if small_grid not in self._periods.keys():
            self._periods[small_grid] = {}

    ###########
    # raw data treatment
    ###########

    def country_zoom(self,input_file,name,var_name,mask_style='lat_weighted',lat_name='lat',lon_name='lon',overwrite=False,**kwargs):
        '''
        zoom input_file to area relevant for the country
        input_file: str: file to be processed
        var_name: str: name of the variable of interest
        mask_style: str: name of the mask used to load the data (see create_mask_country and create_mask_admin)
        time_units: str: time format units
        time_calendar: str: calendar information
        lat_name: str: name of latitude variable
        lon_name: str: name of longitude variable
        overwrite: bool: if True, old file is deleted, new one created
        **kwargs: dict: tags given to the dataset. These tags are going to be interpreted in country_data_object() __init__()
        '''

        #print kwargs
        #out_file=self._working_directory_raw+'/'+input_file.split('/')[-1].replace('.nc','_'+self._iso+'.nc')
        out_file=self._working_directory_raw+'/'+self._iso+'_'+var_name+'_'+name+'.nc4'
        print out_file

        if os.path.isfile(out_file.replace('.nc','_merged.nc'))==False and os.path.isfile(out_file)==False:
            # open file to get information
            print input_file
            lon_in,lat_in,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

            country_mask=self._masks[grid][mask_style][self._iso]
            country_mask=np.ma.getdata(country_mask)
            lat_mask=self._masks[grid]['lat_mask']
            lon_mask=self._masks[grid]['lon_mask']

            # find relevant area (as rectangle)
            lon_mean=np.nanmean(country_mask,0)
            #lons=sorted(np.where(lon_mean!=0)[0])
            lons=sorted(np.where(np.isfinite(lon_mean))[0])

            lat_mean=np.nanmean(country_mask,1)
            #lats=sorted(np.where(lat_mean!=0)[0])
            lats=sorted(np.where(np.isfinite(lat_mean))[0])

            nx,ny=len(lons),len(lats)

            lon=lon_mask[list(lons)]
            lat=lat_mask[list(lats)]

            # zoom to relevant area
            os.system('cdo -O sellonlatbox,'+str(min(lon))+','+str(max(lon))+','+str(min(lat))+','+str(max(lat))+' '+input_file+' '+out_file)

            grid=self._grid_dict[grid]
            in_nc=da.read_nc(out_file)

            datevar=[num2date(tt,units = in_nc['time'].units) for tt in in_nc['time']]
            year=[date.year for date in datevar]
            month=[date.month for date in datevar]
            in_time=np.array([yr+float(mn)/12. for yr,mn in zip(year,month)])
            relevant_steps=np.where(np.isfinite(np.nanmean(in_nc[var_name],axis=(-1,-2))))[0]
            in_time=in_time[relevant_steps]

            print(name)
            if name in self._DATA[grid]:
                print(in_nc[var_name][in_time,:,:].shape)
                self._DATA[grid][name,in_time,:,:]=in_nc[var_name][in_time,:,:]

            else:
                tmp = da.DimArray(axes=[[name],self._time_axis,in_nc[lat_name].values,in_nc[lon_name].values],dims=['name','time','lat','lon'])
                tmp[name,in_time,:,:]=in_nc[var_name][in_time,:,:]
                self._DATA[grid] = da.concatenate((self._DATA[grid],tmp))

    def classify_ensemble(self,filters,grid=None):
        if grid is None and len(self._DATA)==1:
            grid=self._DATA.keys()[0]
        names=self._DATA[grid].name[:]
        for filter in filters:
            names=[xx for xx in names if filter in xx.split('_')]
        return names

    def ensemble_statistic(self,members,ens_name,stat='median',grid=None,write=True):
        if grid is None and len(self._DATA)==1:
            grid=self._DATA.keys()[0]

        member0=self._DATA[grid][members[0],:,:,:]
        tmp = da.DimArray(axes=[[ens_name],self._time_axis,member0.lat,member0.lon],dims=['name','time','lat','lon'])
        tmp[ens_name,:,:,:]=np.nanpercentile(self._DATA[grid][members,:,:,:],[50],axis=0)

        self._DATA[grid] = da.concatenate((self._DATA[grid],tmp))


    def period_statistic(self,period,period_name,stat='mean'):
        for grid in self._DATA:

            self._periods[grid][period_name] = da.DimArray(np.nanmean(self._DATA[grid][:,period[0]:period[1],:,:],axis=1),axes=[self._DATA[grid].name,self._DATA[grid].lat,self._DATA[grid].lon],dims=['name','lat','lon'])

    def period_diff(self,period1,period2,period_name):
        for grid in self._DATA:
            self._periods[grid][period_name] = self._periods[grid][period1] - self._periods[grid][period2]


    def period_events_above_thresh(self,period,period_name,names=None,threshold=0,below=False):
        for grid in self._DATA:
            self._periods[grid][period_name] = np.nanmean(self._DATA[grid][names,period[0]:period[1],:,:],axis=1)*np.nan
            for name in names:
                for y in self._DATA[grid].lat:
                    for x in self._DATA[grid].lon:
                        diff = self._DATA[grid][name,period[0]:period[1],y,x] - threshold
                        if below:
                            self._periods[grid][period_name][name,y,x] = len(np.where(diff<0)[0])/float(diff.shape[1])*100
                        else:
                            self._periods[grid][period_name][name,y,x] = len(np.where(diff>0)[0])/float(diff.shape[1])*100

    def model_agreement(self,ens_name,members,period_name,agreement_name):
        for grid in self._DATA:
            ens=self._periods[grid][period_name][ens_name,:,:]

            agreement=ens.copy()*0
            for member in members:
                agreement+=np.sign(self._periods[grid][period_name][member,:,:])==np.sign(ens)
            agreement[agreement<2./3.*len(members)]=1
            agreement[agreement>=2./3.*len(members)]=np.nan
            tmp = da.DimArray(agreement ,axes=[[agreement_name],self.lat,self.lon],dims=['name','lat','lon'])
            self._periods[grid][period_name] = da.concatenate((self._periods[grid][period_name],tmp))




















    ###########
    # analysis tools
    ###########
