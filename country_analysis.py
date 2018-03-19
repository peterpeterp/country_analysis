# -*- coding: utf-8 -*-
'''
Class to analyze climate data on national (& sub national) scale
Peter Pfleiderer
peter.pfleiderer@climateanalytics.org
'''

import sys,glob,os,itertools,datetime,pickle,subprocess,time
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
from mpl_toolkits.basemap import Basemap
from shapely.geometry import mapping, Polygon, MultiPolygon
import matplotlib.pylab as plt
import matplotlib as mpl

'''
This is a stable but outdated version used in RegioClim
'''


from matplotlib import rc
rc('text', usetex=True)
plt.rcParams["font.family"] = "sans-serif"
plt.style.use('classic')


try:
	from unidecode import unidecode
except:
	pass

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
		self._DATA=[]
		self._regions={}

		if os.path.isdir(self._working_directory)==False:os.system('mkdir '+self._working_directory)
		if os.path.isdir(self._working_directory+'/masks')==False:os.system('mkdir '+self._working_directory+'/masks')
		if os.path.isdir(self._working_directory+'/plots')==False:os.system('mkdir '+self._working_directory+'/plots')
		if os.path.isdir(self._working_directory_raw)==False:os.system('mkdir '+self._working_directory_raw)
		if os.path.isdir(self._working_directory+'/area_average')==False:os.system('mkdir '+self._working_directory+'/area_average')

	def zip_it(self,subfolder=None):
		'''
		Prepares a compressed tar
		subfolder: str: name of the sub-folder that has to be compressed. If none, all the data for the country is compressed
		'''
		actual_path=os.getcwd()
		os.chdir(self._working_directory)
		os.chdir('../')
		if subfolder is None:
			os.system('tar -zcf '+self._iso+'.tar.gz '+self._iso)
		else:
			os.system('tar -zcf '+self._working_directory[0:-1]+'_'+subfolder+'.tar.gz '+self._working_directory.split('/')[-2]+'/'+subfolder)
		os.chdir(actual_path)

	def load_data(self,quiet=True,filename_filter='',load_mask=True,load_raw=True,load_area_averages=True,load_region_polygons=True,load_merged_regions=False):
		'''
		Loads data from existing country_analysis project
		quiet: bool: if quiet, no file-names are printed
		'''
		if load_mask:
			for file in glob.glob(self._working_directory+'/masks/'+self._iso+'*.nc*'):
				# need the country mask first
				if 'admin' not in file.split('_'):
					file_new=self._working_directory+'/masks'+file.split('masks')[-1]
					if quiet==False:print file_new
					self.load_masks(file_new)

			for file in glob.glob(self._working_directory+'/masks/'+self._iso+'*.nc*'):
				if 'admin' in file.split('_'):
					if len(file.split('+'))==1 or load_merged_regions: # exclude merged regions
						file_new=self._working_directory+'/masks'+file.split('masks')[-1]
						if quiet==False:print file_new
						self.load_masks(file_new)


		if load_raw:
			for file in glob.glob(self._working_directory_raw+'/*'+filename_filter+'*'):
				file_new=self._working_directory_raw+file.split('raw'+self._additional_tag)[-1]
				if file_new not in [data.raw_file for data in self._DATA]:
					if quiet==False:print file_new
					nc_out=Dataset(file_new,"r")
					tags={}
					for key,val in zip(nc_out.getncattr('tags_keys').split('**'),nc_out.getncattr('tags_values').split('**')):
						tags[key]=val
					try:
						var_name=tags['original_var_name']
					except:
						var_name=tags['var_name']

					new_data=country_data_object(outer_self=self,**tags)
					new_data.raw_file=file_new
					new_data.add_data(raw=nc_out.variables[var_name][:,:,:],lat=nc_out.variables['lat'][:],lon=nc_out.variables['lon'][:],time=nc_out.variables['time'][:],year=nc_out.variables['year'][:],month=nc_out.variables['month'][:])
					try:
						new_data.add_data(day=nc_out.variables['day'][:])
					except:
						print 'no days'
					new_data.create_time_stamp()


		if load_area_averages:
			for file in glob.glob(self._working_directory+'/area_average/*'):
				mask_style=file.split('-')[-1].split('.')[0]
				name=file.split('-')[-2]

				file_new=self._working_directory+'/area_average'+file.split('area_average')[-1]
				if quiet==False:print file_new
				for data in self._DATA:
					if sorted(data.name.split('_'))==sorted(name.split('_')):
						table=pd.read_csv(file_new,sep=';')
						for key in table.keys():
							if key not in ['time','year','month','index']:
								if mask_style not in data.area_average.keys():	data.area_average[mask_style]={}
								if 'unidecode' in sys.modules:
									data.area_average[mask_style][unidecode(key)]=np.array(table[key])
								if 'unidecode' not in sys.modules:
									data.area_average[mask_style][key]=np.array(table[key])




		# try to load polygons of adm regions for plotting
		if load_region_polygons:
			print 'load regions'
			start_time=time.time()
			self._adm_polygons=self.simplify_adm_polygons(self._working_directory+self._iso+'_adm_shp/'+self._iso+'_adm1',tolerance=0.02)
			for region_name in self._regions.keys():
				if '+' in region_name:
					sub_regs=region_name.split('+')
					self._adm_polygons[region_name]=self._adm_polygons[sub_regs[0]]
					for region in sub_regs[1:]:
						self._adm_polygons[region_name] = \
						self._adm_polygons[region_name].symmetric_difference(self._adm_polygons[region])

			m = Basemap()
			m.readshapefile(self._working_directory+self._iso+'_adm_shp/'+self._iso+'_adm0', 'admin', drawbounds=False)
			name=self._iso
			for shape in m.admin:
				if name in self._adm_polygons.keys():
					self._adm_polygons[name] = \
					self._adm_polygons[name].symmetric_difference(Polygon(shape))
				else:
					self._adm_polygons[name] = Polygon(shape)
			print self._adm_polygons.keys()
			print 'regions loaded '+str(time.time()-start_time)

	def get_historical_extreme_events(self,path):
		'''
		Load csv table with historical extreme events.
		Columns should include: ID, country, start, end, region, description, type
		'''
		table=pd.read_csv(path,sep=';').to_dict()
		event_groups={'flood':['flash floods','floods'],'drought':['drought']}
		extreme_events={}
		for event in ['flood','drought']:
			extreme_events[event]=[]
			for i in table['ID']:
				if table['country'][i]=='SEN':
					if len(set(table['type'][i].replace(' ','').split(',')).intersection(event_groups[event]))>=1:
						extreme_events[event].append({
							'start':table['start'][i],
							'end':table['end'][i],
							'region':table['region'][i],
							'description':table['type'][i],
							})

		self._extreme_events=extreme_events

	def variables(self):
		print '\n'.join(sorted(set([dd.var_name for dd in self._DATA])))

	def datasets(self):
		print '\n'.join(sorted(set([dd.data_type for dd in self._DATA])))

	def short_summary(self):
		'''
		display available variables, datasets and scenarios
		'''
		types=sorted(set([dd.data_type for dd in self._DATA]))
		var_names=sorted(set([dd.var_name for dd in self._DATA]))

		scenarios=[]
		for dd in self._DATA:
			if hasattr(dd,'scenario'):
				scenarios.append(dd.scenario)
		scenarios=sorted(list(set(scenarios)))

		print '--------------Variables------------\n'+'\t'.join(var_names)
		print '**************Datasets*************\n'+'\t'.join(types)
		print '______________Scenarios____________\n'+'\t'.join(scenarios)

	def summary(self):
		'''
		display available variables, datasets and scenarios
		'''
		types=sorted(set([dd.data_type for dd in self._DATA]))
		var_names=sorted(set([dd.var_name for dd in self._DATA]))

		scenarios=[]
		for dd in self._DATA:
			if hasattr(dd,'scenario'):
				scenarios.append(dd.scenario)
		scenarios=sorted(list(set(scenarios)))

		for var_name in var_names:
			print var_name+': '
			for data_type in types:
				if len(self.selection([data_type,var_name],show_selection=False))>0:
					out='\t-'+data_type+': '
					if hasattr(self.selection([data_type,var_name],show_selection=False)[0],'scenario'):
						for scenario in scenarios:
							if len(self.selection([data_type,var_name,scenario],show_selection=False))>0:
								out+=scenario+', '
					print out

	def complete_summary(self,detailed=True):
		'''
		display available variables, datasets and scenarios
		'''
		types=set([dd.data_type for dd in self._DATA])
		var_names=set([dd.var_name for dd in self._DATA])

		for data_type in types:
			print '\n\n***********',data_type,'*********************************'
			for var_name in var_names:
				if len(self.selection([data_type,var_name],show_selection=False))>0:
					print '\n____________',var_name,'___________'
					self.selection([data_type,var_name],detailed)

	def selection(self,filters,show_selection=False,detailed=True):
		'''
		select datasets for further analysis. returns a list of country_data_objects
		filters: list: keywords that all selected elements should share
		show_selection: bool: if True, information about the selected data is shown
		detailed: bool: if True, detailed information about selected data is shown
		'''
		selection=[]
		for data in self._DATA:
			selected=True
			for key in filters:
				if key not in data.all_tags:
					selected=False
			if selected:
				selection.append(data)
		if show_selection==True:
			for count,data in zip(range(len(selection)),selection):
				print count,data.details(detailed)
		return selection

	def find_ensemble(self,filters,quiet=True):
		'''
		selects an ensemble for further analysis. returns a dict containing country_data_objects
		filters: list: keywords that all selected elements should share
		quiet: bool: if False, information about the selected data is shown
		'''
		ensemble={}
		ensemble_mean=None
		ensemble_median=None
		ensemble_min=None
		ensemble_max=None
		ensemble_25=None
		ensemble_75=None
		for data in self._DATA:
			selected=True
			for key in filters:
				if key not in data.all_tags:
					selected=False
			if selected and data.model not in ['ensemble_mean','ensemble_median','ensemble_min','ensemble_max','ensemble_25','ensemble_75']:
				ensemble[data.model]=data
				if quiet==False:
					print data.model+': '+data.name+' '+str(min(data.year))+'-'+str(max(data.year))
			if selected and data.model=='ensemble_mean':
				ensemble_mean=data
			if selected and data.model=='ensemble_median':
				ensemble_median=data
			if selected and data.model=='ensemble_min':
				ensemble_min=data
			if selected and data.model=='ensemble_max':
				ensemble_max=data
			if selected and data.model=='ensemble_25':
				ensemble_25=data
			if selected and data.model=='ensemble_75':
				ensemble_75=data

		return {'models':ensemble,'mean':ensemble_mean,'median':ensemble_median,'min':ensemble_min,'max':ensemble_max,'25':ensemble_25,'75':ensemble_75}

	def unit_conversions(self,selection=None):
		'''
		convert K to deg C
		'''
		if selection is None:
			selection=self._DATA[:]

		for data in selection:
			if data.var_name in ['tas','TXx','tasmax'] or 'tas' in data.var_name.split('_'):
				if np.nanmean(data.raw)>100:
					data.raw-=273.15
				for mask_style in data.area_average.keys():
					for region in data.area_average[mask_style].keys():
						if region in self._regions.values()+[self._iso]:
							if np.nanmean(data.area_average[mask_style][region])>100:
								data.area_average[mask_style][region]-=273.15

	def clean_above_or_below(self,above=999,below=-999):
		'''
		delete weird features
		'''
		for data in self._DATA:
			data.raw[data.raw>above]=np.nan
			data.raw[data.raw<below]=np.nan
			for mask_style in data.area_average.keys():
				for region in data.area_average[mask_style].keys():
					if region in self._regions.keys():
						x=data.area_average[mask_style][region].copy()
						x[x>above]=np.nan
						x[x<below]=np.nan

	def get_warming_slices(self,warming_lvls=[1.5,2],ref_period=[1986,2006],window=21,ref_period_name='ref',model_real_names=None,wlcalculator_path='/Users/peterpfleiderer/Documents/Projects/wlcalculator/app/'):
		'''
		get model specific periods corresponding to global mean temperature warming levels
		GMT_path: path: Path to GMT files (wlcalculator)
		warming_levels: list: GMT warming levels in deg C above preindustrial
		ref_period: list(start year, end year): reference period
		warming_of_ref_period: float: GMT of reference period in deg C above preindustrial
		model_real_names: dict: exact names of models {'given model name in country_analysis class':'corresponding model name in GMT files'}
		'''

		current_path=os.getcwd()
		os.chdir(wlcalculator_path)
		sys.path.append(wlcalculator_path)
		os.system('ls')
		import wacalc.CmipData as CmipData; reload(CmipData)
		import wacalc.hadcrut_warming as hadcrut_warming; reload(hadcrut_warming) # I think I don't need this

		self._warming_slices={}
		for data in self._DATA:
			if hasattr(data,'model'):
				if data.model not in self._warming_slices.keys() and data.model!='ensemble_mean':
					self._warming_slices[data.model]={}
					if data.scenario not in self._warming_slices[data.model].keys():
						self._warming_slices[data.model][data.scenario]={ref_period_name:ref_period}

						# model names from cordex are not explicit!
						if model_real_names is not None:		model_name=model_real_names[data.model]
						if model_real_names is None:		model_name=data.model.lower()

						scenario=data.scenario.replace('4p5','45').replace('2p6','26').replace('6p0','60').replace('8p5','85')

						cmipdata = CmipData.CmipData('CMIP5',[model_name],[scenario])
						cmipdata.get_cmip()
						cmipdata.compute_period( ref_period, [1850,1900], warming_lvls, window=window)
						lvls=cmipdata.exceedance_tm

						for wlvl in lvls.wlevel:
							self._warming_slices[data.model][data.scenario][str(wlvl)]=[round(lvls[scenario][wlvl]-10.),round(lvls[scenario][wlvl]+10.)]

		os.chdir(current_path)

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
				if 'unidecode' in sys.modules:
					self._masks[grid][mask_style][unidecode(name)] = nc_mask.variables[name][:,:]
					if unidecode(name) not in self._regions.keys(): self._regions[unidecode(name).replace(' ','_')]=name
					self.zoom_mask(grid,mask_style,unidecode(name))
				if 'unidecode' not in sys.modules:
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

	def get_admin_polygons(self,shape_file):
		# load shape file
		m = Basemap()
		m.readshapefile(shape_file, 'admin', drawbounds=False)

		# collect all shapes of region
		region_polygons={}
		for shape, region in zip(m.admin, m.admin_info):
			region = {k.lower():v for k,v in region.items()}
			if 'unidecode' in sys.modules: name = unidecode(region['name_1'].decode('utf-8')).replace(' ','_')
			if 'unidecode' not in sys.modules: name = region['name_1'].replace(' ','_')

			if name in region_polygons.keys():
				region_polygons[name] = \
				region_polygons[name].symmetric_difference(Polygon(shape))
			else:
				region_polygons[name] = Polygon(shape)

		return region_polygons

	def simplify_adm_polygons(self,shape_file,tolerance=0.05):
		# load shape file
		m = Basemap()
		m.readshapefile(shape_file, 'admin', drawbounds=False)

		# collect all shapes of region
		region_polygons={}
		for shape, region in zip(m.admin, m.admin_info):
			region = {k.lower():v for k,v in region.items()}
			if 'unidecode' in sys.modules: name = unidecode(region['name_1'].decode('utf-8')).replace(' ','_')
			if 'unidecode' not in sys.modules: name = region['name_1'].replace(' ','_')

			if name in region_polygons.keys():
				region_polygons[name] = \
				region_polygons[name].symmetric_difference(Polygon(shape).simplify(tolerance))
			else:
				region_polygons[name] = Polygon(shape).simplify(tolerance)

		import fiona

		# Define a polygon feature geometry with one attribute
		schema = {
		    'geometry': 'Polygon',
		    'properties': {'name_1': 'str'},
		}

		# Write a new Shapefile
		with fiona.open(shape_file+'_simplified.shp', 'w', 'ESRI Shapefile', schema) as c:
		    for region_name,shape in zip(region_polygons.keys(),region_polygons.values()):
			    c.write({
			        'geometry': mapping(shape),
			        'properties': {'name_1': region_name},
			    })

		return region_polygons

	def merge_adm_regions(self,region_names,new_region_name=None):
		if new_region_name is None:
			single_regions=[]
			for region in region_names:
				for split in region.split('+'):
					single_regions.append(split)
			new_region_name='+'.join(sorted(single_regions))
			self._regions[new_region_name]='+'.join([self._regions[reg] for reg in sorted(single_regions)])
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

	def create_mask_country(self,input_file,var_name,shape_file,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon'):
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

		else:
			grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)

			# load shape file
			m = Basemap()
			m.readshapefile(shape_file, 'admin', drawbounds=False)

			# collect all shapes of country
			for shape, country in zip(m.admin, m.admin_info):
				country = {k.lower():v for k,v in country.items()}
				if (country['iso_a3']==self._iso) | (country['adm0_a3']==self._iso):
					if 'country_polygons' in locals():
						country_polygons = \
						country_polygons.symmetric_difference(Polygon(shape))
					else:
						country_polygons = Polygon(shape)

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
			os.system('rm '+mask_file)
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

	def create_mask_admin(self,input_file,var_name,shape_file=None,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon',regions=None):
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

		else:
			grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)

			'''
			not testes yet!!!
			'''
			# get admin polygons
			if hasattr(self,'_adm_polygons')==False:
				if shape_file is not None:	self._adm_polygons=self.get_admin_polygons(shape_file)
				if shape_file is None: 		return('no shape_file specified')
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
			os.system('rm '+mask_file)
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

			nc_mask.setncattr('original_grid',grid)
			nc_mask.setncattr('mask_style',mask_style)
			nc_mask.close()

			self.zoom_mask(grid,mask_style,name)

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
		lon_=lon_mask[lons]

		lat_mean=np.nanmean(cou_mask,1)
		#lats=np.where(lat_mean!=0)[0]
		lats=sorted(np.where(np.isfinite(lat_mean))[0])
		lat_=lat_mask[lats]

		small_grid=str(len(lat_))+'x'+str(len(lon_))+'_lat_'+str(lat_[0])+'_'+str(lat_[-1])+'_lon_'+str(lon_[0])+'_'+str(lon_[-1])
		if small_grid not in self._masks.keys():	self._masks[small_grid]={}
		if mask_style not in self._masks[small_grid].keys():	self._masks[small_grid][mask_style]={}

		self._masks[small_grid][mask_style][region]=mask[np.ix_(list(lats),list(lons))]
		self._grid_dict[grid]=small_grid

	###########
	# raw data treatment
	###########

	def understand_time_format(self,nc_in=None,time=None,time_units=None,time_calendar=None):
		'''
		interpret time format and store additional information
		nc_in: file_path: input file. If specified time will be read from the file
		time: array: if given, this time array will be treated
		time_units: str: time format units
		time_calendar: str: calendar information
		'''
		if time is None:
			time=nc_in.variables['time'][:]
		# issue with hadgem2 file
		#time[time<0]=-999999
		time=np.delete(time,np.where(time<0))

		datevar = []
		# if specified units and calendar
		if time_units is not None and time_calendar is not None:
			datevar.append(num2date(time,units = time_units,calendar= time_calendar))
		# if no specification
		if time_units is None and time_calendar is None:
			time_unit=nc_in.variables['time'].units
			try:
				cal_temps = nc_in.variables['time'].calendar
				datevar.append(num2date(time,units = time_unit,calendar = cal_temps))
			except:
				datevar.append(num2date(time,units = time_unit))
		# create index variable
		year=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])
		month=np.array([int(str(date).split("-")[1])	for date in datevar[0][:]])
		day=np.array([int(str(date).split("-")[2].split(' ')[0])	for date in datevar[0][:]])

		return(time,year,month,day)

	def create_complete_time_axis(self,data_object,yearmin=None,yearmax=None):
		'''
		create homogeneous time axis
		data_object: country_data_object: data_object for which the axis is created
		yearmin: int: start-year
		yearmax: int: end-year
		'''
		if yearmin is None:	np.nanmin(data_object.year)
		if yearmax is None:	np.nanmax(data_object.year)

		if len(set(data_object.month))==1:
			time_axis=[(year,6,15) for year in np.arange(yearmin,yearmax+1,1)]
		elif len(set(data_object.day))==1:
			time_axis=[]
			for year in np.arange(yearmin,yearmax+1,1):
				for month in sorted(set(data_object.month)):
					time_axis.append((year,month,15))
		elif len(set(data_object.day))>1:
			time_axis=[]
			for year in np.arange(yearmin,yearmax+1,1):
				for month in sorted(set(data_object.month)):
					for day in sorted(set(data_object.day)):
						time_axis.append((year,month,day))

		return(time_axis)

	def fill_gaps_in_time_axis(self,data_object,in_file,out_file):
		'''
		transform data to data with homogeneous time axis. missing time steps will be filled with np.nan
		data_object: country_data_object: data_object to be transformed
		in_file: file path: file to be treated
		out_file: file path: where to store new file
		'''

		# write additional information in copied file
		nc_in=Dataset(in_file,"r")
		raw_data=nc_in.variables[data_object.original_var_name][:,:,:]

		# understand time format and create continuous time axis
		time,year,month,day=self.understand_time_format(nc_in)
		data_object.add_data(time=time,year=year,month=month,day=day)
		data_object.create_time_stamp()
		data_object.add_data(lon=nc_in.variables['lon'][:],lat=nc_in.variables['lat'][:])

		time_axis=self.create_complete_time_axis(data_object,yearmin=np.nanmin(data_object.year),yearmax=np.nanmax(data_object.year))

		out_data=np.zeros([len(time_axis),len(data_object.lat),len(data_object.lon)])*np.nan

		for t in time_axis:
			if t in data_object.time_stamp:
				# some files have more than one value per time step????
				#print t,time_axis[np.where(time_axis==t)[0]],data_object.time_stamp[np.where(data_object.time_stamp==t)[0]]
				out_data[time_axis.index(t),:,:]=raw_data[data_object.time_stamp.index(t),:,:]

		data_object.time_stamp=time_axis
		data_object.convert_time_stamp()
		data_object.add_data(raw=out_data)

		os.system('rm '+out_file)
		nc_out=Dataset(out_file,"w")

		nc_out.createDimension('lat', len(data_object.lat))
		nc_out.createDimension('lon', len(data_object.lon))
		nc_out.createDimension('time', len(time_axis))

		varin=nc_in.variables['lat']
		outVar = nc_out.createVariable('lat', varin.datatype, varin.dimensions)
		outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
		outVar[:] = varin[:]

		varin=nc_in.variables['lon']
		outVar = nc_out.createVariable('lon', varin.datatype, varin.dimensions)
		outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
		outVar[:] = varin[:]

		outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=data_object.time
		outVar.setncatts({'calendar':'proleptic_gregorian','units':'days since 1950-1-1 00:00:00'})

		outVar = nc_out.createVariable('year', 'f', ('time',)) ; outVar[:]=data_object.year
		outVar = nc_out.createVariable('month', 'f', ('time',)) ; outVar[:]=data_object.month
		outVar = nc_out.createVariable('day', 'f', ('time',)) ; outVar[:]=data_object.day

		varin=nc_in.variables[data_object.original_var_name]
		outVar = nc_out.createVariable(data_object.original_var_name, varin.datatype, varin.dimensions)
		outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
		outVar[:] = out_data

		nc_out.setncattr('original_grid',data_object.grid)
		nc_out.setncattr('tags_keys','**'.join(data_object.all_tags_dict.keys()))
		nc_out.setncattr('tags_values','**'.join(data_object.all_tags_dict.values()))

		nc_out.close()
		nc_in.close()
		os.system('rm '+in_file)

	def country_zoom(self,input_file,var_name,mask_style='lat_weighted',time_units=None,time_calendar=None,lat_name='lat',lon_name='lon',overwrite=False,**kwargs):
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
		if 'given_var_name' in kwargs.keys():
			out_file=self._working_directory_raw+'/'+self._iso+'_'+kwargs['given_var_name']
		if 'given_var_name' not in kwargs.keys():
			out_file=self._working_directory_raw+'/'+self._iso+'_'+var_name
		for key in sorted(kwargs.keys()):
			if key not in ['time_format','var_name','given_var_name']:
				out_file+='_'+kwargs[key]
		out_file+='.nc4'
		print out_file

		if overwrite:
			os.system('rm '+out_file)
			os.system('rm '+out_file.replace('.nc','_merged.nc'))

		if os.path.isfile(out_file):
			nc_out=Dataset(out_file,"r")
			new_data=country_data_object(outer_self=self,var_name=var_name,grid=nc_out.getncattr('original_grid'),**kwargs)
			new_data.raw_file=out_file
			new_data.add_data(raw=nc_out.variables[var_name][:,:,:],lat=nc_out.variables['lat'][:],lon=nc_out.variables['lon'][:],time=nc_out.variables['time'][:],year=nc_out.variables['year'][:],month=nc_out.variables['month'][:],day=nc_out.variables['day'][:])
			new_data.create_time_stamp()
			nc_out.close()

		if os.path.isfile(out_file.replace('.nc','_merged.nc')):
			nc_out=Dataset(out_file.replace('.nc','_merged.nc'),"r")
			new_data=country_data_object(outer_self=self,var_name=var_name,grid=nc_out.getncattr('original_grid'),**kwargs)
			new_data.raw_file=out_file.replace('.nc','_merged.nc')
			new_data.add_data(raw=nc_out.variables[var_name][:,:,:],lat=nc_out.variables['lat'][:],lon=nc_out.variables['lon'][:],time=nc_out.variables['time'][:],year=nc_out.variables['year'][:],month=nc_out.variables['month'][:],day=nc_out.variables['day'][:])
			new_data.create_time_stamp()
			nc_out.close()

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
			os.system('cdo -O sellonlatbox,'+str(min(lon))+','+str(max(lon))+','+str(min(lat))+','+str(max(lat))+' '+input_file+' '+out_file.replace('.nc','_tmp.nc'))

			new_data=country_data_object(outer_self=self,var_name=var_name,grid=grid,**kwargs)
			new_data.raw_file=out_file
			self.fill_gaps_in_time_axis(new_data,out_file.replace('.nc','_tmp.nc'),out_file)

	def hist_merge(self):
		'''
		merge historical and rcp-scenarios for each model
		'''
		# why are there files missing if I go through self._DATA only once???
		for data in self._DATA+self._DATA:
			if hasattr(data,'model'):
				data__=self.selection([data.model,data.var_name,data.data_type],show_selection=True)
				for hist in data__[:]:
					if hist.scenario.lower() in ['hist','historical']:
						delete_hist=False
						print '------- merging ',hist.model,hist.var_name,'-----------'
						print ' '.join([d_d.raw_file for d_d in data__])

						for ddd in data__[:]:
							if ddd.scenario.lower() not in ['hist','historical']:
								out_file = ddd.raw_file.replace('.nc','_merged.nc')
								tmp_file = ddd.raw_file.replace('.nc','_merged_tmp.nc')


								os.system('cdo -O mergetime '+ddd.raw_file+' '+hist.raw_file+' '+tmp_file)
								os.system('rm '+ddd.raw_file)

								ddd.raw_file=out_file
								self.fill_gaps_in_time_axis(ddd,tmp_file,out_file)

								delete_hist=True

						if delete_hist:
							self._DATA.remove(hist)
							os.system('rm '+hist.raw_file)

	def ensemble_statistic(self,stat='median',selection=None,write=True):
		'''
		identify ensembles and compute ensemble mean for each ensemble
		'''
		if selection is None:
			selection=self._DATA[:]

		# why do I need this????? see hist_merge
		for i in range(2):
			remaining=selection[:]
			for data in remaining:
				if hasattr(data,'model'):
					if data.model not in ['ensemble_median','ensemble_mean','ensemble_min','ensemble_max']:
						#print '__ ',data.name
						ensemble=self.find_ensemble([data.data_type,data.var_name,data.scenario])
						if ensemble[stat] is None:
							if stat=='median':command='cdo -O enspctl,50 '
							if stat=='min':command='cdo -O enspctl,0 '
							if stat=='max':command='cdo -O enspctl,100 '
							if stat=='25':command='cdo -O enspctl,25 '
							if stat=='75':command='cdo -O enspctl,75 '
							if stat=='mean':command='cdo -O ensmean '

							for member in ensemble['models'].values():
								command+=member.raw_file+' '

							out_file=self._working_directory_raw+'/'+self._iso+'_'+member.var_name+'_'+member.data_type+'_ensemble-'+stat+'_'+member.scenario+'.nc4'
							os.system(command+out_file)

							tags_=member.all_tags_dict.copy()
							tags_['model']='ensemble_'+stat
							tags_['raw_file']=out_file

							new_data=country_data_object(outer_self=self,**tags_)
							new_data.add_data(raw=Dataset(out_file).variables[member.original_var_name][:,:,:],lat=member.lat.copy(),lon=member.lon.copy())
							new_data.add_data(time=member.time.copy(),year=member.year.copy(),month=member.month.copy(),day=member.day.copy())
							new_data.create_time_stamp()

							os.system('ncatted -h -a tags_values,global,o,c,"'+'**'.join(new_data.all_tags_dict.values())+'" '+out_file)

							if write==False: os.system('rm '+out_file)

	def ensemble_mean(self,silent=True):
		'''
		identify ensembles and compute ensemble mean for each ensemble
		'''
		# why do I need this????? see hist_merge
		for i in range(2):
			remaining=self._DATA[:]
			for data in remaining:
				if hasattr(data,'model'):
					if data.model!='ensemble_mean':
						if silent==False:
							print '---------------- ensemble mean',data.var_name,'--------------------'
						ensemble=self.find_ensemble([data.data_type,data.var_name,data.scenario])
						if ensemble['mean'] is None:

							time_min,time_max=[],[]
							for member in ensemble['models'].values():
								remaining.remove(member)
								time_min.append(min(member.time_stamp))
								time_max.append(max(member.time_stamp))
							time_axis=self.create_complete_time_axis(member,yearmin=max(time_min)[0],yearmax=min(time_max)[0])
							#time_axis=np.around(np.arange(max(time_min),min(time_max),member.time_step),decimals=4)
							ensemble_mean=np.zeros([len(ensemble['models'].values()),len(time_axis),len(member.lat),len(member.lon)])*np.nan

							for member,i in zip(ensemble['models'].values(),range(len(ensemble['models'].values()))):
								for t in member.time_stamp:
									if t in time_axis:
										ensemble_mean[i,time_axis.index(t),:,:]=member.raw[member.time_stamp.index(t),:,:]

							ensemble_mean=np.nanmean(ensemble_mean,axis=0)

							tags_=member.all_tags_dict.copy()
							tags_['model']='ensemble_mean'
							out_file=self._working_directory_raw+'/'+self._iso+'_'+member.var_name+'_'+member.data_type+'_ensemble-mean_'+member.scenario+'.nc4'
							tags_['raw_file']=out_file

							new_data=country_data_object(outer_self=self,**tags_)
							new_data.add_data(raw=ensemble_mean,lat=member.lat,lon=member.lon)

							new_data.time_stamp=time_axis
							new_data.convert_time_stamp()

							# write ensemble_mean to file
							if silent==False:
								print member.name
								print member.raw_file
							nc_in=Dataset(member.raw_file,"r")
							os.system('rm '+out_file)
							nc_out=Dataset(out_file,"w")

							nc_out.createDimension('lat', len(new_data.lat))
							nc_out.createDimension('lon', len(new_data.lon))
							nc_out.createDimension('time', len(time_axis))

							varin=nc_in.variables['lat']
							outVar = nc_out.createVariable('lat', varin.datatype, varin.dimensions)
							outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
							outVar[:] = varin[:]

							varin=nc_in.variables['lon']
							outVar = nc_out.createVariable('lon', varin.datatype, varin.dimensions)
							outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
							outVar[:] = varin[:]

				 			outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=new_data.time
				 			outVar.setncatts({'calendar':'proleptic_gregorian','units':'days since 1950-1-1 00:00:00'})

				 			outVar = nc_out.createVariable('year', 'f', ('time',)) ; outVar[:]=new_data.year
				 			outVar = nc_out.createVariable('month', 'f', ('time',)) ; outVar[:]=new_data.month
 							outVar = nc_out.createVariable('day', 'f', ('time',)) ; outVar[:]=new_data.day

 							if silent==False:
 								print new_data.original_var_name
							varin=nc_in.variables[new_data.original_var_name]
							outVar = nc_out.createVariable(new_data.original_var_name, varin.datatype, varin.dimensions)
							outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
							outVar[:] = ensemble_mean[:,:,:]

				 			nc_out.setncattr('original_grid',new_data.grid)
				 			nc_out.setncattr('tags_keys','**'.join(new_data.all_tags_dict.keys()))
				 			nc_out.setncattr('tags_values','**'.join(new_data.all_tags_dict.values()))

							nc_out.close()
							nc_in.close()

	###########
	# analysis tools
	###########

	def area_average(self,mask_style='lat_weighted',selection=None,overwrite=False,regions=None):
		'''
		compute or load countrywide (and region-wide) averages
		mask_style: str: weighting used to compute countrywide averages
		selection: list: list of country_data_objects for which the are is computed
		overwrite: bool: if True, old file is deleted, new one created
		'''
		if selection is None:
			selection=self._DATA[:]

		for data in selection:
			out_file=self._working_directory+'/area_average/country_mean-'+data.name+'-'+mask_style+'.csv'

			if regions is None:
				regions=self._masks[data.grid][mask_style].keys()

			if mask_style not in data.area_average.keys():	data.area_average[mask_style]={}
			data.area_average[mask_style]['out_file']=out_file

			if os.path.isfile(out_file) and overwrite==False:
				country_mean_csv=pd.read_csv(out_file,sep=';')
				for key in country_mean_csv.keys():
					if key not in ['time','year','month']:
						data.area_average[mask_style][key]=np.array(country_mean_csv[key])

			else:
				# prepare table
				country_mean_csv = pd.DataFrame(index=range(len(data.time)))
				country_mean_csv['time']=data.time
				country_mean_csv['month']=data.month
				country_mean_csv['year']=data.year
				country_mean_csv['day']=data.day

			# load input data
			var_in=data.raw.copy()
			try:	# handle masked array
				masked=np.ma.getmask(var_in)
				var_in=np.ma.getdata(var_in)
				var_in[masked]=np.nan
			except: pass


			# get mask
			for name in regions:
				if name not in data.area_average[mask_style].keys() or overwrite:
					mask=self._masks[self._grid_dict[data.grid]][mask_style][name]

					country_area=np.where(mask>0)

					data.area_average[mask_style][name]=data.time.copy()*np.nan
					for i in range(len(data.time)):
						var_of_area=var_in[i,:,:][country_area]

						# NA handling: sum(mask*var)/sum(mask) for the area where var is not NA
						not_missing_in_var=np.where(np.isfinite(var_of_area))[0]	# np.where()[0] because of array([],)
						if len(not_missing_in_var)>0:
							data.area_average[mask_style][name][i]=sum(mask[country_area][not_missing_in_var]*var_of_area[not_missing_in_var])/sum(mask[country_area][not_missing_in_var])

						#print data.area_average[mask_style][name][i],sum(mask[country_area][not_missing_in_var]*var_of_area[not_missing_in_var]),sum(mask[country_area][not_missing_in_var])

					country_mean_csv[name]=data.area_average[mask_style][name]

			# delete index columns
			for key in country_mean_csv.keys():
				if key.split('.')[0]=='index': country_mean_csv.drop(key, axis=1, inplace=True)

			# save as csv
			country_mean_csv.to_csv(out_file,na_rep='NaN',sep=';',index_label='index',encoding='utf-8')

	def period_statistics(self,method='mean',threshold=None,below=False,selection=None,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]},ref_name='ref'):
		'''
		computes time averages for each grid-cell for a given period. possible to compute mean or frequency above or below threshold
		method: str: if 'mean', the period mean is computed. if frequency above threshold has to be calculated, method is a given name
		threshold: float | dict | array: threshold can be a value, an array with a different threshold for each grid-cell or a dict with a different threshold array for each time in year
		below: bool: if True, frequency below threshold instead of frequency above threshold
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		selection: list: list of country_data_objects for which the are is computed
		'''

		if selection is None:
			selection=self._DATA

		for data in selection:
			compute=True
			if 'ensemble_mean' in data.all_tags:
				compute=False

			# handle warming slices and fixed periods
			if depth(periods)==2:
				local_periods=periods
			if depth(periods)>2:
				try:	local_periods=periods[data.model][data.scenario]
				except:	compute=False

			if compute:
				if hasattr(data,'period')==False:	data.period={}
				if method not in data.period.keys(): data.period[method]={}

				if data.time_format=='monthly':		seasons=self._seasons
				if data.time_format=='10day':		seasons=self._seasons
				if data.time_format=='yearly':		seasons={'year':range(1,13)}

				if isinstance(threshold, (np.ndarray)) and threshold is not None:
					threshold_new=data.raw.copy()*np.nan
					for tt in range(len(data.time)):
						threshold_reshaped[tt,:,:]=threshold

				if isinstance(threshold, (float,int)) and threshold is not None:
					threshold_reshaped=np.ones(data.raw.shape)*threshold

				if isinstance(threshold, (dict)) and threshold is not None:
					threshold_reshaped=data.raw.copy()*np.nan
					for tt,i in zip(data.time_in_year,range(len(data.time_in_year))):
						threshold_reshaped[i,:,:]=threshold[tt]

				for period_name,period in zip(local_periods.keys(),local_periods.values()):
					for mons_in_sea,sea in zip(seasons.values(),seasons.keys()):
						if sea not in data.period[method].keys(): data.period[method][sea]={}
						# find relevant time steps
						relevant_years=np.where((data.year>=period[0]) & (data.year<period[1]))[0]
						relevenat_time_steps=data.get_relevant_time_steps_in_season(mons_in_sea,relevant_years)

						n_steps=len(relevenat_time_steps)
						if n_steps>0:
							# compute average
							if method=='mean':
								data.period[method][sea][period_name]=np.nanmean(data.raw[relevenat_time_steps,:,:],axis=0)
							if method=='max':
								data.period[method][sea][period_name]=np.nanmax(data.raw[relevenat_time_steps,:,:],axis=0)
							if method=='min':
								data.period[method][sea][period_name]=np.nanmin(data.raw[relevenat_time_steps,:,:],axis=0)

							# threshold exceeding
							if threshold is not None:
								data.period[method][sea][period_name]=np.zeros([data.raw.shape[1],data.raw.shape[2]])*np.nan
								for y in range(data.raw.shape[1]):
									for x in range(data.raw.shape[2]):
										diff=data.raw[relevenat_time_steps,y,x]-threshold_reshaped[relevenat_time_steps,y,x]
										if below==True:
											data.period[method][sea][period_name][y,x]=len(np.where(diff<0)[0])/float(n_steps)*100
										if below==False:
											data.period[method][sea][period_name][y,x]=len(np.where(diff>0)[0])/float(n_steps)*100

						else:
							#print 'years missing for',period_name,'in',data.name
							data.period[method][sea][period_name]=np.zeros([data.raw.shape[1],data.raw.shape[2]])*np.nan

				# period diff
				self.period_statistic_diff(data,method,sea,ref_name=ref_name)

	def period_statistic_diff(self,data,method,sea,ref_name='ref'):
		# period diff
		for sea in data.period[method].keys():
			for period_name in data.period[method][sea].keys():
				if period_name!=ref_name and ref_name in data.period[method][sea].keys() and period_name.split('_')[0]!='diff':
					data.period[method][sea]['diff_'+period_name+'-'+ref_name]=data.period[method][sea][period_name]-data.period[method][sea][ref_name]
					data.period[method][sea]['diff_relative_'+period_name+'-'+ref_name]=(data.period[method][sea][period_name]-data.period[method][sea][ref_name])/data.period[method][sea][ref_name]*100

	def period_model_agreement(self,ref_name='ref'):
		'''
		computes ensemble mean and model agreement for period means and frequencies
		'''
		remaining=self.selection(['ensemble_mean'],show_selection=False)

		for data in remaining:
			ensemble=self.find_ensemble([data.data_type,data.var_name,data.scenario])
			if hasattr(ensemble['mean'],'period')==False:	ensemble['mean'].period={}
			if hasattr(ensemble['mean'],'agreement')==False:	ensemble['mean'].agreement={}

			member=ensemble['models'].values()[0]

			if hasattr(member,'period'):
				if member.time_format=='monthly':		seasons=self._seasons
				if member.time_format=='10day':		seasons=self._seasons
				if member.time_format=='yearly':		seasons={'year':range(1,13)}

				for method in member.period.keys():
					ensemble['mean'].period[method]={}
					ensemble['mean'].agreement[method]={}
					for sea in seasons.keys():
						ensemble['mean'].period[method][sea]={}
						ensemble['mean'].agreement[method][sea]={}

						# ensemble mean
						for period in member.period[method][sea].keys():
							ensemble['mean'].period[method][sea][period]=member.period[method][sea][period].copy()*0
							for member in ensemble['models'].values():
								ensemble['mean'].period[method][sea][period]+=member.period[method][sea][period]
							ensemble['mean'].period[method][sea][period]/=float(len(ensemble['models'].values()))

						# model agreement
						for period in ensemble['mean'].period[method][sea].keys():
							if len(period.split('-'))>1 and period.split('-')[-1]==ref_name:
								agreement=ensemble['mean'].period[method][sea][period].copy()*0
								for member in ensemble['models'].values():
									agreement+=np.sign(member.period[method][sea][period])==np.sign(ensemble['mean'].period[method][sea][period])

								agreement[agreement<2./3.*len(ensemble['models'].values())]=0
								agreement[agreement>=2./3.*len(ensemble['models'].values())]=1
								ensemble['mean'].agreement[method][sea][period]=agreement

	def period_mean(self,selection=None,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]}):
		'''
		this function uses period_statistics and period_model_agreement to compute period mean
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		selection: list: list of country_data_objects for which the are is computed
		'''
		self.period_statistics(method='mean',selection=selection,periods=periods)
		self.period_model_agreement()

	def frequency_above_threshold(self,method,threshold=None,selection=None,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]}):
		'''
		this function uses period_statistics and period_model_agreement to compute the frequency above threshold for a period
		method: str: if 'mean', the period mean is computed. if frequency above threshold has to be calculated, method is a given name
		threshold: float | dict | array: threshold can be a value, an array with a different threshold for each grid-cell or a dict with a different threshold array for each time in year
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		selection: list: list of country_data_objects for which the are is computed
		'''
		self.period_statistics(method=method,threshold=threshold,below=False,selection=selection,periods=periods)
		self.period_model_agreement()

	def frequency_below_threshold(self,method,threshold=None,selection=None,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]}):
		'''
		see frequency_above_threshold with frequencies below threshold (below=True)
		'''
		self.period_statistics(method=method,threshold=threshold,below=True,selection=selection,periods=periods)
		self.period_model_agreement()

	def annual_cycle(self,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]},mask_style='lat_weighted',selection=None,regions=None,ref_name='ref'):
		'''
		computes the annual cycle in different periods and regions using the averages from area_average()
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		mask_style: str: weighting used to compute countrywide averages
		selection: list: list of country_data_objects for which the are is computed
		'''

		if selection is None:
			selection=self._DATA

		for data in selection:
			compute=True
			if 'ensemble_mean' in data.all_tags:
				compute=False

			if regions is None:
				regions=self._masks[data.grid][mask_style].keys()

			# handle warming slices and fixed periods
			if depth(periods)==2:
				local_periods=periods
			if depth(periods)>2:
				try:	local_periods=periods[data.model][data.scenario]
				except:	compute=False

			if compute:
				data.steps_in_year=sorted(set(data.time_in_year_num))
				if hasattr(data,'annual_cycle')==False:	data.annual_cycle={}
				if mask_style not in data.annual_cycle.keys():		data.annual_cycle[mask_style]={}
				#data.annual_cycle={key : val for key,val in zip(local_periods.keys(),local_periods.values())}

				for region in regions:
					if region not in data.annual_cycle[mask_style].keys():
						data.annual_cycle[mask_style][region]={}
					for period_name,period in zip(local_periods.keys(),local_periods.values()):
						annual_cycle=[]
						#print period_name,period,data.steps_in_year
						for tt in data.steps_in_year:
							relevant_time=np.where((data.year>=period[0]) & (data.year<period[1]) & (data.time_in_year_num==tt))[0]
							if len(relevant_time)>1:
								annual_cycle.append(np.nanmean(data.area_average[mask_style][region][relevant_time]))
						if len(annual_cycle)==len(data.steps_in_year):
							data.annual_cycle[mask_style][region][period_name]=np.array(annual_cycle)
						else:
							data.annual_cycle[mask_style][region][period_name]=np.zeros([len(data.steps_in_year)])*np.nan

						self.annual_cycle_diff(data=data,mask_style=mask_style,region=region,ref_name=ref_name)

	def annual_cycle_diff(self,data,mask_style,region,ref_name='ref'):
		for period_name in data.annual_cycle[mask_style][region].keys():
			if period_name!=ref_name and ref_name in data.annual_cycle[mask_style][region].keys():
				data.annual_cycle[mask_style][region]['diff_'+period_name+'-'+ref_name]=data.annual_cycle[mask_style][region][period_name]-data.annual_cycle[mask_style][region][ref_name]
				data.annual_cycle[mask_style][region]['diff_relative_'+period_name+'-'+ref_name]=(data.annual_cycle[mask_style][region][period_name]-data.annual_cycle[mask_style][region][ref_name])/data.annual_cycle[mask_style][region][ref_name]*100

	def annual_cycle_ensemble_mean(self,regions=None,selection=None):
		'''
		computes time averages for each grid-cell for a given period
		'''
		if selection is None:
			selection=self._DATA[:]

		remaining=selection[:]

		for data in remaining:
			compute=False
			if hasattr(data,'model'):
				if data.model!='ensemble_mean':
					compute=True

			if compute:
				if data.model in self.find_ensemble([data.data_type,data.var_name,data.scenario])['models'].keys():
					ensemble=self.find_ensemble([data.data_type,data.var_name,data.scenario])

					ensemble['mean'].steps_in_year=sorted(set(ensemble['mean'].time_in_year_num))

					for member in ensemble['models'].values():
						remaining.remove(member)

					# annual cycle ensemble mean
					if hasattr(member,'annual_cycle'):
						if hasattr(ensemble['mean'],'annual_cycle')==False:			ensemble['mean'].annual_cycle={}
						for mask_style in member.annual_cycle.keys():
							if mask_style not in ensemble['mean'].annual_cycle.keys():		ensemble['mean'].annual_cycle[mask_style]={}
							if regions is None:
								regions=self._masks[data.grid][mask_style].keys()
							for region in regions:
								ensemble['mean'].annual_cycle[mask_style][region]={}
								for period in member.annual_cycle[mask_style][region].keys():
									ensemble['mean'].annual_cycle[mask_style][region][period]=np.nanmean(np.vstack([member.annual_cycle[mask_style][region][period] for member in ensemble['models'].values()]),axis=0)

class country_data_object(object):
	'''
	objects of this class store information of one data-source.
	'''
	def __init__(SELF,outer_self,**kwargs):
		'''
		When loading a new data-input to the country_analysis-object meta-information is stored here
		outer_self: country_analysis object: functions of this class sometimes use variables and functions of the country_analysis class
		**kwargs: dict: meta information is given through kwargs
		The following keywords will be recognized when given in kwargs:
			mandatory:	raw_file, grid, time_format, var_name
			optional: 	given_var_name,  scenario, model, data_type
		'''

		SELF.outer_self=outer_self
		outer_self._DATA.append(SELF)
		SELF._iso=outer_self._iso

		if 'raw_file' in kwargs.keys():	SELF.raw_file=outer_self._working_directory_raw+'/'+kwargs['raw_file'].split('raw/')[-1]
		if 'grid' in kwargs.keys():	SELF.grid=kwargs['grid']
		if 'time_format' in kwargs.keys():	SELF.time_format=kwargs['time_format']

		SELF.original_var_name=kwargs['var_name']
		if 'given_var_name' in kwargs.keys():	SELF.var_name=kwargs['given_var_name']
		if 'given_var_name' not in kwargs.keys():	SELF.var_name=kwargs['var_name']

		if 'data_type' in kwargs.keys():	SELF.data_type=kwargs['data_type']

		if 'scenario' in kwargs.keys():
			kwargs['scenario']=kwargs['scenario'].replace('.','p').replace('45','4p5').lower()
			SELF.scenario=kwargs['scenario']

		if 'model' in kwargs.keys():	SELF.model=kwargs['model']

		SELF.area_average={}

		SELF.all_tags_dict=kwargs
		SELF.all_tags_dict.pop('raw_file', None)
		SELF.all_tags=kwargs.values()
		if 'given_var_name' in kwargs.keys():	SELF.all_tags.remove(kwargs['var_name'])


		SELF.name=SELF.var_name
		for key in sorted(kwargs.keys()):
			if key not in ['raw_file','grid','var_name','given_var_name']:
				SELF.name+='_'+kwargs[key]

	def add_data(SELF,**kwargs):
		'''
		add data to object
		Recognized data names: raw, time, year, month, day, lat, lon
		'''
		if 'raw' in kwargs.keys():
			tmp=kwargs['raw']
			tmp[np.isfinite(tmp)==False]=np.nan
			SELF.raw=tmp
		if 'time' in kwargs.keys():	SELF.time=kwargs['time']
		if 'year' in kwargs.keys():	SELF.year=kwargs['year']
		if 'month' in kwargs.keys():	SELF.month=kwargs['month']
		if 'day' in kwargs.keys():	SELF.day=kwargs['day']
		if 'lat' in kwargs.keys():	SELF.lat=kwargs['lat']
		if 'lon' in kwargs.keys():	SELF.lon=kwargs['lon']

	def create_time_stamp(SELF):
		'''
		creates a custom time_stamp and time_in_year which is required for the annual cycle
		'''
		if SELF.time_format=='yearly':
			SELF.day=SELF.year.copy()*0+15
			SELF.month=SELF.year.copy()*0+6
		elif SELF.time_format=='monthly':
			SELF.day=SELF.year.copy()*0+15
		elif SELF.time_format=='10day':
			SELF.day[SELF.day<10]=5
			SELF.day[np.where((SELF.day>10) & (SELF.day<20))]=15
			SELF.day[SELF.day>20]=25

		# uniform time-stamp
		SELF.time_stamp=[(SELF.year[i],SELF.month[i],SELF.day[i]) for i in range(len(SELF.year))]
		SELF.time_stamp_num=np.array([SELF.year[i]+((SELF.month[i]-1)*30 + SELF.day[i])/365. for i in range(len(SELF.year))])

		# time in year (required for annual_cycle)
		SELF.time_in_year=np.array([(SELF.month[i],SELF.day[i]) for i in range(len(SELF.year))])
		SELF.time_in_year_num=np.array([((SELF.month[i]-1)*30 + SELF.day[i])/365. for i in range(len(SELF.year))])
		SELF.plot_time=np.array([SELF.year[i]+((SELF.month[i]-1)*30 + SELF.day[i])/365. for i in range(len(SELF.year))])

	def convert_time_stamp(SELF):
		'''
		converts custom time_stamp from create_time_stamp to year, month, day and time
		'''
		SELF.year=np.array([int(t[0]) for t in SELF.time_stamp])
		SELF.month=np.array([int(t[1]) for t in SELF.time_stamp])
		SELF.day=np.array([int(t[2]) for t in SELF.time_stamp])

		SELF.time=np.array([(datetime.datetime(int(SELF.year[i]),int(SELF.month[i]),int(SELF.day[i])) - datetime.datetime(1950,1,1)).days for i in range(len(SELF.year))])

	def show_periods(SELF):
		periods=[]
		for method in SELF.period.keys():
			for sea in SELF.outer_self._seasons.keys():
				periods+=SELF.period[method][sea].keys()

		print '\n'.join(sorted(set(periods)))

	def cdo_command(SELF,new_var_name,cdo_command='yearmax'):
		'''
		computes annual max using 'cdo yearmax'. this function could be used as a template for other functions of this kind
		new_var_name: str: given new var name
		'''
		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean':
				return 'you might want to apply the cdo command to the individual models and compute the ensemble mean afterwards'

		kwargs=SELF.all_tags_dict
		kwargs['var_name']=SELF.original_var_name
		kwargs['given_var_name']=new_var_name
		kwargs['time_format']='yearly'

		out_file=SELF.raw_file.replace(SELF.var_name,new_var_name)
		os.system('cdo '+cdo_command+' '+SELF.raw_file+' '+out_file.replace('.nc','_tmp.nc'))

		new_data=country_data_object(outer_self=SELF.outer_self,**kwargs)
		new_data.raw_file=out_file
		new_data.original_var_name=SELF.original_var_name
		new_data.add_data(lat=SELF.lat,lon=SELF.lon)
		SELF.outer_self.fill_gaps_in_time_axis(new_data,out_file.replace('.nc','_tmp.nc'),out_file)

	def seas_max(SELF,new_var_name,season):
		'''
		computes annual max using 'cdo yearmax'. this function could be used as a template for other functions of this kind
		new_var_name: str: given new var name
		'''
		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean':
				return 'not a good idea to apply yearmax on ensemble mean'

		kwargs=SELF.all_tags_dict
		kwargs['var_name']=SELF.original_var_name
		kwargs['given_var_name']=new_var_name
		kwargs['time_format']='yearly'

		out_file=SELF.raw_file.replace(SELF.var_name,new_var_name)
		os.system('cdo -seasmax '+SELF.raw_file+' '+out_file.replace('.nc','_tmp1.nc'))
		os.system('cdo -selseas,'+season+' '+out_file.replace('.nc','_tmp1.nc')+' '+out_file.replace('.nc','_tmp.nc'))
		os.system('rm '+out_file.replace('.nc','_tmp1.nc'))

		new_data=country_data_object(outer_self=SELF.outer_self,**kwargs)
		new_data.raw_file=out_file
		new_data.original_var_name=SELF.original_var_name
		new_data.add_data(lat=SELF.lat,lon=SELF.lon)
		SELF.outer_self.fill_gaps_in_time_axis(new_data,out_file.replace('.nc','_tmp.nc'),out_file)

	def year_max(SELF,new_var_name):
		'''
		computes annual max using 'cdo yearmax'. this function could be used as a template for other functions of this kind
		new_var_name: str: given new var name
		'''
		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean':
				return 'not a good idea to apply yearmax on ensemble mean'

		year_new=np.array(sorted(set(SELF.year)))
		dat=np.zeros([len(year_new),len(SELF.lat),len(SELF.lon)])*np.nan
		for yr,i in zip(year_new,range(len(year_new))):
			dat[i,:,:]=np.nanmax(SELF.raw[np.where(SELF.year==yr)[0],:,:],axis=0)

		kwargs=SELF.all_tags_dict
		kwargs['var_name']=SELF.original_var_name
		kwargs['given_var_name']=new_var_name
		kwargs['time_format']='yearly'

		out_file=SELF.raw_file.replace(SELF.var_name,new_var_name)
		os.system('cdo yearmax '+SELF.raw_file+' '+out_file.replace('.nc','_tmp.nc'))

		new_data=country_data_object(outer_self=SELF.outer_self,**kwargs)
		new_data.raw_file=out_file
		new_data.original_var_name=SELF.original_var_name
		new_data.add_data(lat=SELF.lat,lon=SELF.lon)
		SELF.outer_self.fill_gaps_in_time_axis(new_data,out_file.replace('.nc','_tmp.nc'),out_file)

	def details(SELF,detailed=True):
		'''
		show details of country_data_object
		detailed: bool: if True, detailed details are shown
		'''
		if detailed:
			text=SELF.name.replace('_',' ')
			text+='\t\t\tcoverage: '+str(int(min(SELF.year)))+'-'+str(int(max(SELF.year)))
			# if hasattr(SELF,'period'):
			# 	for key in SELF.period.keys():
			# 		text+='\nperiod '+key+': '+','.join(sorted(SELF.period[key][].keys()))
			return(text)
		else:
			return(SELF.name)

	def get_relevant_time_steps_in_season(SELF,months_in_season,relevant_years=None):
		'''
		get relevant time indices
		months_in_season: list: months (as 1:12)
		relevant years: list: years as int
		'''
		if relevant_years is None:
			relevant_years=range(len(SELF.year))
		relevenat_time_steps=[]
		for yr in relevant_years:
			if int(SELF.month[yr]) in months_in_season:
				relevenat_time_steps.append(yr)
		return(relevenat_time_steps)

	def plot_map(SELF,to_plot,limits=None,ax=None,out_file=None,title='',polygons=None,grey_area=None,color_bar=True,color_label='',color_palette=plt.cm.plasma,color_range=None,highlight_region=None,show_all_adm_polygons=False,show_region_names=False,show_merged_region_names=False):
		'''
		this function creates a map
		to_plot: np.ndarray: values to plot
		lat: array: latitudes
		lon: array: longitudes
		limits: list: [xmin,xmax,ymin,ymax]
		ax: sub_plot_ax: axes on which to plot. If None new figure is created
		out_file: file path: path where the plot is saved
		title: str: plot title
		polygons: list: list of polygon objects to be overlayed
		grey_area: np.ndarray: same dimensions as to_plot. 0 will be colored in gray, 1 will be left transparent
		color_bar: bool: if True, a color bar is added to the plot
		color_label: str: label of color bar
		color_palette: plt.cm.object: color palette of the plot
		color_range: list: [min-value,max-value]
		'''
		if ax is None: show=True
		if ax is not None: show=False

		lat=SELF.lat.copy()
		lon=SELF.lon.copy()

		if to_plot=='empty':
			to_plot=np.zeros([len(lat),len(lon)])*np.nan
			color_range=[0,1]

		if ax is None:
			asp=(len(lon)/float(len(lat)))**0.5
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(3*asp+1,3/asp))


		# handle limits
		if limits is None:
			half_lon_step=abs(np.diff(lon.copy(),1)[0]/2)
			half_lat_step=abs(np.diff(lat.copy(),1)[0]/2)
			limits=[np.min(lon)-half_lon_step,np.max(lon)+half_lon_step,np.min(lat)-half_lat_step,np.max(lat)+half_lat_step]


		# handle 0 to 360 lon
		if max(lon)>180:
			problem_start=np.where(lon>180)[0][0]
			new_order=np.array(range(problem_start,len(lon))+range(0,problem_start))
			to_plot=to_plot[:,new_order]
			lon=lon[new_order]
			lon[lon>180]-=360

		if limits[0]>180:limits[0]-=360
		if limits[1]>180:limits[1]-=360


		m = Basemap(ax=ax,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
		m.drawmapboundary(fill_color='1.')

		# get color_range
		if color_range is None:
			color_range=[np.min(to_plot[np.isfinite(to_plot)]),np.max(to_plot[np.isfinite(to_plot)])]

		x,y=lon.copy(),lat.copy()
		x-=np.diff(x,1)[0]/2.
		y-=np.diff(y,1)[0]/2.
		x=np.append(x,[x[-1]+np.diff(x,1)[0]])
		y=np.append(y,[y[-1]+np.diff(y,1)[0]])
		x,y=np.meshgrid(x,y)
		im = m.pcolormesh(x,y,np.ma.masked_invalid(to_plot),cmap=color_palette,vmin=color_range[0],vmax=color_range[1])


		# mask some grid-cells
		if grey_area is not None:
			to_plot=grey_area.copy()
			to_plot[to_plot==0]=0.5
			to_plot[to_plot==1]=np.nan
			to_plot=np.ma.masked_invalid(to_plot)
			im2 = m.pcolormesh(x,y,to_plot,cmap=plt.cm.Greys,vmin=0,vmax=1)

		# show coastlines and borders
		m.drawcoastlines()
		m.drawstates()
		m.drawcountries()

		# add polygons
		if polygons is None and hasattr(SELF.outer_self,'_adm_polygons'):
			polygons=SELF.outer_self._adm_polygons
			if show_all_adm_polygons:
				for name,polygon in polygons.items():
					try:
						x,y=polygon.exterior.xy
						m.plot(x,y,color='black',linewidth=0.5)
						if show_region_names:
							if (show_merged_region_names or len(name.split('+'))==1) and name!=SELF._iso:
								ctr=polygons[name].centroid.xy
								plt.text(ctr[0][0],ctr[1][0],SELF.outer_self._regions[name],horizontalalignment='center',verticalalignment='center',fontsize=8)

					except Exception,e:
						try:
							areas=[]
							for shape in polygon:
								areas.append(shape.area)
								x,y=shape.exterior.xy
								m.plot(x,y,color='black',linewidth=0.5)
							if show_region_names:
								if (show_merged_region_names or len(name.split('+'))==1) and name!=SELF._iso:
									ctr=polygons[name][areas.index(max(areas))].centroid.xy
									plt.text(ctr[0][0],ctr[1][0],SELF.outer_self._regions[name],horizontalalignment='center',verticalalignment='center',fontsize=8)
						except:
							# no idea what goes wrong here, its an issue for some CPV shapes...
							pass

			# highlight one region
			if highlight_region is not None:
				if highlight_region!=SELF._iso:
					try:
						x,y=polygons[highlight_region].exterior.xy
						m.plot(x,y,color='black',linewidth=1.5)
						m.plot(x,y,color='yellow',linewidth=1.5,linestyle='--')
					except:
						for shape in polygons[highlight_region]:
							x,y=shape.exterior.xy
							m.plot(x,y,color='black',linewidth=1.5)
							m.plot(x,y,color='yellow',linewidth=1.5,linestyle='--')


		# add colorbar
		if color_bar==True:
			cb = m.colorbar(im,'right', size="5%", pad="2%")
			tick_locator = mpl.ticker.MaxNLocator(nbins=5)
			cb.locator = tick_locator
			cb.update_ticks()
			cb.set_label(color_label, rotation=90)

		if title != '':
			ax.set_title(title.replace('_',' '))

		if out_file is None and show==True:plt.show()
		if out_file is not None:plt.savefig(out_file)

		return(im,color_range)

	def display_map(SELF,period=None,method='mean',season='year',show_agreement=True,limits=None,ax=None,out_file=None,title=None,polygons=None,color_bar=True,color_label=None,color_palette=None,color_range=None,time=None,highlight_region=None,show_all_adm_polygons=True,show_region_names=False,show_merged_region_names=False):
		'''
		plot maps of data.
		period: str: if  the averag over a period is to be plotted specify the period name
		method: str: 'mean' for period mean. 'custom_name' for frequency above/below threshold (see period_statistcs())
		season: str: name of the season
		show_agreement: bool: if True, if ensemble_mean is plotted, model-disagreement is masked in gray
		see plot_map() for remaining kwargs
		'''

		grey_area=None

		if period is None:
			if time is None:
				time=int(len(SELF.time)/2)
				print 'no time specified. '+str(int(SELF.month[time]))+'/'+str(int(SELF.year[time]))+' selected'
				to_plot=SELF.raw[time,:,:].copy()
				if title is None:title='_'.join([SELF.name]+[str(int(SELF.month[time])),'/',str(int(SELF.year[time]))])
		else:
			to_plot=SELF.period[method][season][period].copy()
			if title is None:title=SELF.name+'_'+method+'_'+period+'_'+season

			# mask
			mask=SELF.outer_self._masks[SELF.outer_self._grid_dict[SELF.grid]]['lat_weighted'][SELF.outer_self._iso].copy()
			to_plot[np.isfinite(mask)==False]=np.nan

			if hasattr(SELF,'agreement') and show_agreement:
				if period in SELF.agreement[method][season].keys():
					grey_area=SELF.agreement[method][season][period].copy()
					grey_area[np.isfinite(mask)==False]=1


		if len(np.where(np.isfinite(to_plot)==False)[0])==len(to_plot.flatten()):
			print 'nothing to plot here'
			return(None,None)

		if color_label is None:color_label=SELF.var_name.replace('_',' ')

		if color_range is None:
			if np.sign(np.nanpercentile(to_plot,[10]))==np.sign(np.nanpercentile(to_plot,[90])):
				color_range=np.nanpercentile(to_plot,[10,90])
			else:
				abs_boundary=max(abs(np.nanpercentile(to_plot,[10,90])))
				color_range=[-abs_boundary,abs_boundary]

		if color_palette is None:
			if SELF.var_name in ['pr','RX1','year_RX5']:
				if np.mean(color_range)==0:						color_palette=plt.cm.RdBu
				elif np.mean(color_range)<0:					color_palette=plt.cm.Reds_r
				elif np.mean(color_range)>0:					color_palette=plt.cm.Blues
			elif SELF.var_name in ['tas','TXx']:
				if np.mean(color_range)==0:						color_palette=plt.cm.RdBu_r
				elif np.mean(color_range)<0:					color_palette=plt.cm.Blues_r
				elif np.mean(color_range)>0:					color_palette=plt.cm.Reds
			else:
				if np.mean(color_range)==0:						color_palette=plt.cm.RdBu
				elif np.mean(color_range)<0:					color_palette=plt.cm.Reds
				elif np.mean(color_range)>0:					color_palette=plt.cm.Blues_r


		im,color_range=SELF.plot_map(to_plot,color_bar=color_bar,color_label=color_label,color_palette=color_palette,color_range=color_range,grey_area=grey_area,limits=limits,ax=ax,out_file=out_file,title=title,polygons=polygons,highlight_region=highlight_region,show_all_adm_polygons=show_all_adm_polygons,show_region_names=show_region_names,show_merged_region_names=show_merged_region_names)
		return(im,color_range)

	def display_mask(SELF,mask_style=None,region=None):
		'''
		show a map of the used mask
		mask_style: str: mask_style to be shown (see create_mask_country or create_mask_admin)
		region: str: region to be shown
		'''

		if region is None:	region=SELF._iso

		outer_self=SELF.outer_self
		if mask_style is None:
			all_masks=glob.glob(outer_self._working_directory+'/masks/'+outer_self._iso+'_'+SELF.grid+'*')
			mask_styles=[mask_file.split('_')[-2]+'_weighted' for mask_file in all_masks]
			mask_style=mask_styles[0]
			print 'mask_style: '+mask_style+'\nother available mask_styles: '+', '.join(mask_styles[1:-1])

		try:
			to_plot=outer_self._masks[outer_self._grid_dict[SELF.grid]][mask_style][region]
			lat=outer_self._masks[SELF.grid]['lat_mask'].copy()
			lon=outer_self._masks[SELF.grid]['lon_mask'].copy()
			half_lon_step=abs(np.diff(lon.copy(),1)[0]/2)
			half_lat_step=abs(np.diff(lat.copy(),1)[0]/2)
			cou_mask=outer_self._masks[SELF.grid][mask_style][SELF._iso]
			relevant_lats=lat[np.where(np.isfinite(cou_mask))[0]]
			relevant_lons=lon[np.where(np.isfinite(cou_mask))[1]]
			limits=[np.min(relevant_lons)-half_lon_step,np.max(relevant_lons)+half_lon_step,np.min(relevant_lats)-half_lat_step,np.max(relevant_lats)+half_lat_step]
			SELF.plot_map(to_plot,title=SELF.grid+' '+mask_style,color_label='importance for area average',color_palette=plt.cm.plasma,limits=limits)
		except:
			print 'no mask has been created for '+mask_style+' and '+region

	def plot_transients(SELF,mask_style='lat_weighted',season='year',region=None,running_mean_years=1,ax=None,out_file=None,title=None,ylabel=None,label='',color='blue',y_range=None,x_range=[1960,2090],ref_period=None,shading_range=None,shading_opacity=0.2,plot_median=False,show_all_models=False,offset=0.0):
		'''
		plot transient of countrywide average
		mask_style: str: weighting used to compute countrywide averages
		running_mean: int: years to be averaged in moving average
		ax: subplot: subplot on which the map will be plotted
		out_file: str: location where the plot is saved
		title: str: title of the plot
		ylabel: str: labe to put on y-axis
		show: logical: show the subplot?
		'''

		if ax is not None:
			show=False

		if ax is None:
			show=True
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

		if region is None:
			region=SELF.outer_self._iso

		if SELF.time_format=='monthly':
			running_mean=running_mean_years*len(SELF.outer_self._seasons[season])
			relevenat_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season])

		if SELF.time_format=='10day':
			running_mean=running_mean_years*len(SELF.outer_self._seasons[season])*3
			relevenat_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season])

		if SELF.time_format=='yearly':
			running_mean=running_mean_years
			relevenat_time_steps=range(len(SELF.year))
			if season!='year':
				return(0)

		if ref_period is None:
			y=running_mean_func(SELF.area_average[mask_style][region][relevenat_time_steps], running_mean)+offset
		if ref_period is not None:
			ref_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season],np.where((SELF.year>=ref_period[0]) & (SELF.year<ref_period[1]))[0])
			ref_mean=np.nanmean(SELF.area_average[mask_style][region][ref_time_steps])
			y=(running_mean_func(SELF.area_average[mask_style][region][relevenat_time_steps], running_mean)-ref_mean)/ref_mean*100+offset

		if plot_median==False:
			ax.plot(SELF.plot_time[relevenat_time_steps],y,linestyle='-',label=label,color=color)

		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean':
				ensemble=SELF.outer_self.find_ensemble([SELF.data_type,SELF.var_name,SELF.scenario])

				time_axis=ensemble['mean'].time_stamp_num

				ensemble_range=np.zeros([len(ensemble['models'].values()),len(time_axis[relevenat_time_steps])])*np.nan

				for member,i in zip(ensemble['models'].values(),range(len(ensemble['models'].values()))):
					if ref_period is None:
						member_runmean=running_mean_func(member.area_average[mask_style][region][relevenat_time_steps], running_mean)+offset
					if ref_period is not None:
						ref_mean=np.nanmean(member.area_average[mask_style][region][ref_time_steps])
						member_runmean=(running_mean_func(member.area_average[mask_style][region][relevenat_time_steps], running_mean)-ref_mean)/ref_mean*100+offset

					for t in time_axis:
						ensemble_range[i,np.where(time_axis[relevenat_time_steps]==t)]=member_runmean[np.where(member.time_stamp_num[relevenat_time_steps]==t)]

					if show_all_models:	ax.plot(SELF.plot_time[relevenat_time_steps],ensemble_range[i,:],linestyle='--',label=member.model,color=color)

				if shading_range is not None:
					ax.fill_between(SELF.plot_time[relevenat_time_steps],np.nanpercentile(ensemble_range,shading_range[0],axis=0),np.nanpercentile(ensemble_range,shading_range[1],axis=0),alpha=shading_opacity,color=color)
				if plot_median:
					ax.plot(SELF.plot_time[relevenat_time_steps],np.nanpercentile(ensemble_range,50,axis=0),color=color,linestyle='-',label=label)


		if ylabel is None:ylabel=SELF.var_name.replace('_',' ')
		ax.set_ylabel(ylabel)
		if title is None:title=SELF.name.replace('_',' ')
		ax.set_title(title)

		if y_range is not None:
			ax.set_ylim(y_range)

		if x_range is not None:
			ax.set_xlim(x_range)

		if show==True:ax.legend(loc='best')
		if out_file is None and show==True:plt.show()
		if out_file is not None:plt.savefig(out_file)

		return(1)

	def plot_ensemble_transients(SELF,mask_style='lat_weighted',season='year',region=None,running_mean_years=1,ax=None,out_file=None,title=None,ylabel=None,label='',color='blue',y_range=None,x_range=[1960,2100],ref_period=None,shading_range=None,shading_opacity=0.2,plot_median=False):
		'''
		plot transient of countrywide average
		mask_style: str: weighting used to compute countrywide averages
		running_mean: int: years to be averaged in moving average
		ax: subplot: subplot on which the map will be plotted
		out_file: str: location where the plot is saved
		title: str: title of the plot
		ylabel: str: labe to put on y-axis
		show: logical: show the subplot?
		'''

		'''
		restructure the compeletely messy prol_transients
		'''

		if ax is not None:
			show=False

		if ax is None:
			show=True
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

		if region is None:
			region=SELF.outer_self._iso

		if SELF.time_format=='monthly':
			running_mean=running_mean_years*len(SELF.outer_self._seasons[season])
			relevenat_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season])

		if SELF.time_format=='10day':
			running_mean=running_mean_years*len(SELF.outer_self._seasons[season])*3
			relevenat_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season])

		if SELF.time_format=='yearly':
			running_mean=running_mean_years
			relevenat_time_steps=range(len(SELF.year))
			if season!='year':
				return(0)

		if ref_period is None:
			y=running_mean_func(SELF.area_average[mask_style][region][relevenat_time_steps], running_mean)
		if ref_period is not None:
			ref_time_steps=SELF.get_relevant_time_steps_in_season(SELF.outer_self._seasons,np.where((SELF.year>=ref_period[0]) & (SELF.year<ref_period[1]))[0])
			ref_mean=np.nanmean(SELF.area_average[mask_style][region][ref_time_steps])
			y=(running_mean_func(SELF.area_average[mask_style][region][relevenat_time_steps], running_mean)-ref_mean)/ref_mean*100

		if plot_median==False:
			ax.plot(SELF.plot_time[relevenat_time_steps],y,linestyle='-',label=label,color=color)

		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean' and shading_range is not None:
				ensemble=SELF.outer_self.find_ensemble([SELF.data_type,SELF.var_name,SELF.scenario])

				time_axis=ensemble['mean'].time_stamp_num

				ensemble_range=np.zeros([len(ensemble['models'].values()),len(time_axis)])*np.nan

				for member,i in zip(ensemble['models'].values(),range(len(ensemble['models'].values()))):
					if SELF.time_format=='yearly':	relevenat_time_steps=range(len(member.year))
					if SELF.time_format!='yearly':	relevenat_time_steps=member.get_relevant_time_steps_in_season(SELF.outer_self._seasons[season])

					if ref_period is None:
						member_runmean=running_mean_func(member.area_average[mask_style][region][relevenat_time_steps], running_mean)
					if ref_period is not None:
						ref_mean=np.nanmean(member.area_average[mask_style][region][ref_time_steps])
						member_runmean=(running_mean_func(member.area_average[mask_style][region][relevenat_time_steps], running_mean)-ref_mean)/ref_mean*100
						print member_runmean

					for t in member.time_stamp_num:
						ensemble_range[i,np.where(time_axis==t)]=member_runmean[np.where(member.time_stamp_num==t)]

				ax.fill_between(time_axis,np.nanpercentile(ensemble_range,shading_range[0],axis=0),np.nanpercentile(ensemble_range,shading_range[1],axis=0),alpha=shading_opacity,color=color)
				if plot_median:
					ax.plot(time_axis,np.nanpercentile(ensemble_range,50,axis=0),color=color,linestyle='-',label=label)


		if ylabel is None:ylabel=SELF.var_name.replace('_',' ')
		ax.set_ylabel(ylabel)
		if title is None:title=SELF.name.replace('_',' ')
		ax.set_title(title)

		if y_range is not None:
			ax.set_ylim(y_range)

		if x_range is not None:
			ax.set_xlim(x_range)

		if show==True:ax.legend(loc='best')
		if out_file is None and show==True:plt.show()
		if out_file is not None:plt.savefig(out_file)

		return(1)

	def plot_annual_cycle(SELF,mask_style='lat_weighted',region=None,period=None,ax=None,out_file=None,title=None,ylabel=None,label='',color='blue',xlabel=True,shading_range=None,shading_opacity=0.2,plot_median=False):
		'''
		plot transient of countrywide average
		meta_data: list of strs: meta information required to acces data
		mask_style: str: weighting used to compute countrywide averages
		running_mean: int: years to be averaged in moving average
		ax: subplot: subplot on which the map will be plotted
		out_file: str: location where the plot is saved
		title: str: title of the plot
		ylabel: str: labe to put on y-axis
		show: logical: show the subplot?
		'''

		if ax is not None:
			show=False

		if ax is None:
			show=True
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

		if period is None:
			period=[min(SELF.year),max(SELF.year)]

		if region is None:
			region=SELF._iso

		if plot_median==False:
			ax.plot(SELF.steps_in_year,SELF.annual_cycle[mask_style][region][period],linestyle='-',label=label,color=color)

		if hasattr(SELF,'model'):
			if SELF.model=='ensemble_mean' and shading_range is not None:
				ensemble=SELF.outer_self.find_ensemble([SELF.data_type,SELF.var_name,SELF.scenario])
				ensemble_annual_cycle=np.vstack([member.annual_cycle[mask_style][region][period] for member in ensemble['models'].values()])

				ax.fill_between(SELF.steps_in_year,np.nanpercentile(ensemble_annual_cycle,shading_range[0],axis=0),np.nanpercentile(ensemble_annual_cycle,shading_range[1],axis=0),alpha=shading_opacity,color=color)
				if plot_median:
					ax.plot(SELF.steps_in_year,np.nanpercentile(ensemble_annual_cycle,50,axis=0),color=color,linestyle='-',label=label)

		if xlabel==True:
			ax.set_xticks(np.arange(1/24.,25./24.,1/12.))
			ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
		if xlabel==False:
			ax.set_xticks([0.5])
			ax.set_xticklabels([''])

		if ylabel is None:ylabel=SELF.var_name.replace('_',' ')
		ax.set_ylabel(ylabel)
		if title is None:title=SELF.name.replace('_',' ')+' '+region
		ax.set_title(title)

		if show==True:ax.legend(loc='best')
		if out_file is None and show==True:plt.show()
		if out_file is not None:plt.savefig(out_file)

	def historical_index_validation(SELF,extreme_type='flood',mask_style='lat_weighted',region=None,ax=None,out_file=None,title=None,ylabel=None,label='',color='blue',y_range=None):

		if ax is not None:
			show=False

		if ax is None:
			show=True
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

		if region is None:
			region=SELF.outer_self._iso

		ax.fill_between(SELF.time_stamp,SELF.time_stamp.copy()*0,SELF.area_average[mask_style][region],linestyle='-',color='black')

		lower,upper=ax.get_ylim()

		for event in SELF.outer_self._extreme_events[extreme_type]:
			if region in event['region'] or event['region']=='n.a.':
				print event
				yr,mn,dy=event['start'].split('-')
				start_time=int(yr)+(30*int(mn)+int(dy))/365.
				yr,mn,dy=event['end'].split('-')
				end_time=int(yr)+(30*int(mn)+int(dy))/365.
				print start_time,end_time
				ax.fill_between([start_time,end_time],[lower,lower],[upper,upper],linestyle='-',color='blue',alpha=0.5)

		if show==True:ax.legend(loc='best')
		if out_file is None and show==True:plt.show()
		if out_file is not None:plt.savefig(out_file,dpi=300)






	# def seasonal_statistic(self,method='mean',selection=None,overwrite=False,seasons=None):
	# 	if selection is None:
	# 		selection=self._DATA

	# 	if seasons is None:
	# 		seasons=self._seasons

	# 	for data in selection:
	# 		if hasattr(data,'seasonal')==False:
	# 			data.seasonal={'year_axis':np.array(sorted(set(data.year)))}
	# 		data.seasonal[method]={}
	# 		for months,season in zip(seasons.values(),seasons.keys()):
	# 			data.seasonal[method][season]=np.zeros([len(set(data.year)),len(data.lat),len(data.lon)])*np.nan
	# 			for yr,t in zip(sorted(set(data.year)),range(len(set(data.year)))):
	# 				relevenat_time_steps=data.get_relevant_time_steps_in_season(months,np.where(data.year==yr)[0])
	# 				if method=='max':
	# 					data.seasonal['max'][season][t,:,:]=np.nanmax(data.raw[relevenat_time_steps,:,:],axis=0)
	# 				if method=='min':
	# 					data.seasonal['min'][season]['min'][t,:,:]=np.nanmin(data.raw[relevenat_time_steps,:,:],axis=0)
	# 				if method=='mean':
	# 					#print data.raw[relevenat_time_steps,:,:],np.nanmean(data.raw[relevenat_time_steps,:,:],axis=0)
	# 					data.seasonal['mean'][season][t,:,:]=np.nanmean(data.raw[relevenat_time_steps,:,:],axis=0)

	# def area_average_seasonal_statistic(self,method='mean',mask_style='lat_weighted',selection=None,overwrite=False,regions=None,seasons=None):

	# 	if selection is None:
	# 		selection=self._DATA

	# 	if seasons is None:
	# 		seasons=self._seasons

	# 	if regions is None:
	# 		regions=self._masks[data.grid][mask_style].keys()

	# 	for data in selection:
	# 		for season,season_name in zip(seasons.values(),seasons.keys()):
	# 			for name in regions:
	# 				data.average[mask_style][name]['seasonal'][season]['max']=data.year.copy()*np.nan
	# 				for yr,t in zip(data.year,range(len(data.year))):
	# 					relevenat_time_steps=data.get_relevant_time_steps_in_season(season,[yr])
	# 					if method=='max':
	# 						data.average[mask_style][name]['seasonal'][season]['max'][t]=np.nanmax(data.average[mask_style][name]['year'][relevenat_time_steps])
	# 					if method=='min':
	# 						data.average[mask_style][name]['seasonal'][season]['min'][t]=np.nanmin(data.average[mask_style][name]['year'][relevenat_time_steps])
	# 					if method=='mean':
	# 						data.average[mask_style][name]['seasonal'][season]['mean'][t]=np.nanmean(data.average[mask_style][name]['year'][relevenat_time_steps])
