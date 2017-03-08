'''
Class to cliamte data on the national (& sub national) scale
Peter Pfleiderer
peter.pfleiderer@climateanalytics.org
'''

import sys,glob,os,itertools,datetime,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import random as random
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pylab as plt 
import matplotlib as mpl
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator
import seaborn as sns
from matplotlib.colors import ListedColormap

class country_analysis(object):

	def __init__(self,iso,working_directory):
		self._iso=iso
		self._working_directory=working_directory+iso

		self._masks={}
		os.system('mkdir '+self._working_directory+'/masks')

		if os.path.isdir(working_directory)==False:os.system('mkdir '+working_directory)
		if os.path.isdir(self._working_directory+'/plots')==False:os.system('mkdir '+self._working_directory+'/plots')


	def create_mask_country(self,input_file,var_name,shape_file,mask_style='lat_weighted',pop_mask_file=''):
		'''
		create country mask
		input_file: str: location of example input data (required for the identification of the grid)
		var_name: str: variable name of input file
		shape_file: str: location of the shape_file used to identify country borders
		mask_style: str: name under which the mask will be stored (important for further analysis)
		pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
		'''

		# get information about grid of input data
		nc_in=Dataset(input_file,'r')
		try:
			lat = nc_in.variables['lat'][:]								
			lon = nc_in.variables['lon'][:].squeeze()
		except:
			lat = nc_in.variables['latitude'][:]								
			lon = nc_in.variables['longitude'][:].squeeze()			
		# formerly shift_lon
		if max(lon)>200:	lon_shift=-180.0
		else:				lon_shift=0.0	
		lon+=lon_shift

		nx = len(lon)	;	ny = len(lat)
		grid=str(ny)+'x'+str(nx)

		if grid not in self._masks.keys():self._masks[grid]={}
		if mask_style not in self._masks[grid].keys():self._masks[grid][mask_style]={}

		mask_file=self._working_directory+'/masks/'+self._iso+'_'+grid+'_'+mask_style+'.nc4'

		if os.path.isfile(mask_file):
			# load existing mask
			nc_mask=Dataset(mask_file,'r')
			self._masks[grid][mask_style][self._iso] = nc_mask.variables[self._iso][:,:]  
			self._masks[grid]['lat_mask'] = nc_mask.variables['lat'][:]  
			self._masks[grid]['lon_mask'] = nc_mask.variables['lon'][:]  

		else:
			# loop over the grid to get grid polygons
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
			xmin, xmax, ymin, ymax = country_polygons.bounds
			ext = [(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin),(xmin,ymin)]
			ext_poly = Polygon(ext)

			# load population mask		
			if pop_mask_file=='':	
				pop_mask = np.ones((len(lat),len(lon)))
			else:
				# regrid population mask 
				mygrid=open(self._working_directory+'/masks/'+grid+'.txt','w')
				mygrid.write('gridtype=lonlat\nxsize='+str(len(lon))+'\nysize='+str(len(lat))+'\nxfirst='+str(lon[0])+'\nxinc='+str(np.mean(np.diff(lon,1)))+'\nyfirst='+str(lat[0])+'\nyinc='+str(np.mean(np.diff(lat,1))))
				mygrid.close()
				os.system('cdo remapbil,'+self._working_directory+'/masks/'+grid+'.txt '+pop_mask_file+' '+self._working_directory+'/masks/'+grid+'_'+mask_style+'.nc')	
				nc_pop_mask = Dataset(self._working_directory+'/masks/'+grid+'_'+mask_style+'.nc')
				pop_mask = np.array(nc_pop_mask.variables['mask'][:,:]).squeeze()
				pop_mask = np.roll(pop_mask,shift,axis=1)

			# compute overlap
			overlap = np.zeros((ny,nx))
			for i in range(nx):
				for j in range(ny):
					# check gridcell is relevant
					if grid_polygons[i,j].intersects(ext_poly):
						# get fraction of grid-cell covered by polygon
						intersect = grid_polygons[i,j].intersection(country_polygons).area/grid_polygons[i,j].area*country_polygons.area
						if pop_mask_file!='':
							# population weighting
							overlap[j,i] = intersect*pop_mask[j,i]
						if mask_style=='lat_weighted':
							# multiply overlap with latitude weighting
							overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

			# renormalize overlap to get sum(mask)=1
			overlap_zwi=overlap.copy()
			overlap_sum=sum(overlap_zwi.flatten())
			if overlap_sum!=0:
				output=np.zeros(overlap.shape)
				output=overlap/overlap_sum
				# mask zeros
				output[output==0]=np.nan
				output=np.ma.masked_invalid(output)
				# shift back to original longitudes
				self._masks[grid][mask_style][self._iso]=np.roll(output,shift,axis=1)

			# save mask
			print mask_file
			nc_mask=Dataset(mask_file,'w')
			nc_mask.createDimension('lat', len(lat))
			nc_mask.createDimension('lon', len(lon))
 			outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
 			outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')
 			outVar = nc_mask.createVariable(self._iso, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][self._iso][:,:]

 			# close nc_files
			nc_in.close()
			nc_mask.close()

	def create_mask_admin(self,input_file,var_name,shape_file,mask_style='lat_weighted',pop_mask_file=''):
		'''
		create country mask
		input_file: str: location of example input data (required for the identification of the grid)
		var_name: str: variable name of input file
		shape_file: str: location of the shape_file used to identify country borders
		mask_style: str: name under which the mask will be stored (important for further analysis)
		pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
		'''

		# get information about grid of input data
		nc_in=Dataset(input_file,'r')
		try:
			lat = nc_in.variables['lat'][:]								
			lon = nc_in.variables['lon'][:].squeeze()
		except:
			lat = nc_in.variables['latitude'][:]								
			lon = nc_in.variables['longitude'][:].squeeze()			
		# formerly shift_lon
		if max(lon)>200:	lon_shift=-180.0
		else:				lon_shift=0.0	
		lon+=lon_shift

		nx = len(lon)	;	ny = len(lat)
		grid=str(ny)+'x'+str(nx)

		if grid not in self._masks.keys():self._masks[grid]={}
		if mask_style not in self._masks[grid].keys():self._masks[grid][mask_style]={}

		mask_file=self._working_directory+'/masks/'+self._iso+'_admin_'+grid+'_'+mask_style+'.nc4'

		if os.path.isfile(mask_file):
			# load existing mask
			nc_mask=Dataset(mask_file,'r')
			self._masks[grid]['lat_mask'] = nc_mask.variables['lat'][:]  
			self._masks[grid]['lon_mask'] = nc_mask.variables['lon'][:] 
			# get all variables (regions)
			for name in nc_mask.variables.keys():
				if name not in ['lat','lon']:
					self._masks[grid][mask_style][name] = nc_mask.variables[name][:,:]  
 

		else:
			# loop over the grid to get grid polygons
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
			
			# here loop over shape objects in shapefile

			# load shape file
			m = Basemap()
			m.readshapefile(shape_file, 'admin', drawbounds=False)

			# collect all shapes of region
			region_polygons={}
			for shape, region in zip(m.admin, m.admin_info):
				region = {k.lower():v for k,v in region.items()}	
				name = region['name_1']
				if name in region_polygons.keys():
					region_polygons[name] = \
					region_polygons[name].symmetric_difference(Polygon(shape))
				else:
					region_polygons[name] = Polygon(shape)


			# get boundaries for faster computation
			xmins, xmaxs, ymins, ymaxs = [], [], [], []
			for name in region_polygons.keys():
				bounds=region_polygons[name].bounds
				xmins.append(bounds[0])
				xmaxs.append(bounds[1])
				ymins.append(bounds[2])
				ymaxs.append(bounds[3])
			xmin, xmax, ymin, ymax = min(xmins), max(xmaxs), min(ymins), max(ymaxs)
			ext = [(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin),(xmin,ymin)]
			ext_poly = Polygon(ext)

			# prepare outputfile
			os.system('rm '+mask_file)
			nc_mask=Dataset(mask_file,'w')
			nc_mask.createDimension('lat', len(lat))
			nc_mask.createDimension('lon', len(lon))
 			outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
 			outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')

			for name in region_polygons.keys():
				region_shape = region_polygons[name]
				overlap = np.zeros((ny,nx))
				for i in range(nx):
					for j in range(ny):
						# check gridcell is relevant
						if grid_polygons[i,j].intersects(ext_poly):
							# get fraction of grid-cell covered by polygon
							intersect = grid_polygons[i,j].intersection(region_shape).area/grid_polygons[i,j].area*region_shape.area
							if pop_mask_file!='':
								# population weighting
								overlap[j,i] = intersect*pop_mask[j,i]
							if mask_style=='lat_weighted':
								# multiply overlap with latitude weighting
								overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

				# renormalize overlap to get sum(mask)=1
				overlap_zwi=overlap.copy()
				overlap_sum=sum(overlap_zwi.flatten())
				if overlap_sum!=0:
					output=np.zeros(overlap.shape)
					output=overlap/overlap_sum
					# mask zeros
					output[output==0]=np.nan
					output=np.ma.masked_invalid(output)
					# shift back to original longitudes
					self._masks[grid][mask_style][name]=np.roll(output,shift,axis=1)
					outVar = nc_mask.createVariable(name, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][name][:,:]

			# close nc_files
			nc_in.close()
			nc_mask.close()

	def country_zoom(self,in_file,var_name,meta_data=['CMIP5','rcp2p6','hadgem2-es'],mask_style='lat_weighted',time_units=None,time_calendar=None):
		'''
		zoom input_file to area relevant for the country
		in_file: type str: file to be processed
		out_file: type str: filepath where output is going to be stored
		var: type str: variable name in netcdf file
		iso: type str: iso3c of country
		mask_path: type str: path to where the masks are stored
		'''

		if '_data' not in dir(self): 
			self._data={}
			self._meta=[]
		if os.path.isdir(self._working_directory+'/raw')==False:os.system('mkdir '+self._working_directory+'/raw')

		out_file=self._working_directory+'/raw/'+in_file.split('/')[-1].replace('.nc','_'+self._iso+'.nc')

		if os.path.isfile(out_file):
			nc_out=Dataset(out_file,"r")

			tmp = self._data
			for meta_info in meta_data:
				if meta_info not in tmp:tmp[meta_info]={}
				tmp=tmp[meta_info]

			tmp['data']=nc_out.variables[var_name][:,:,:]
			tmp['lat']=nc_out.variables['lat'][:]
			tmp['lon']=nc_out.variables['lon'][:]
			tmp['time']=nc_out.variables['time'][:]
			tmp['year']=nc_out.variables['year'][:]
			tmp['month']=nc_out.variables['month'][:]
			tmp['grid']=nc_out.variables[var_name].getncattr('grid')

			self._meta.append(meta_data)

			nc_out.close()


		else:
			# open file to get information
			print in_file
			nc_in=Dataset(in_file,"r")
			lat_in=nc_in.variables['lat'][:]
			lon_in=nc_in.variables['lon'][:]
			grid=str(len(lat_in))+'x'+str(len(lon_in))

			country_mask=self._masks[grid][mask_style][self._iso]
			lat_mask=self._masks[grid]['lat_mask']
			lon_mask=self._masks[grid]['lon_mask']
			
			# check whether lon and lat are similarly defined in mask and in_file
			if (lat_in in lat_mask) == False:
				lat_mask=lat_mask[::-1]
				country_mask=country_mask[::-1,:]
			if (lon_in in lon_mask) == False:
				lon_mask=lon_mask[::-1]
				country_mask=country_mask[:,::-1]

			# find relevant area (as rectangle)
			lon_mean=np.mean(country_mask,0)
			lons=sorted(np.where(lon_mean!=0)[0])

			lat_mean=np.mean(country_mask,1)
			lats=sorted(np.where(lat_mean!=0)[0])

			nx,ny=len(lons),len(lats)

			var_in=nc_in.variables[var_name][:,list(lats),list(lons)]


			try:	# handle masked array
				masked=np.ma.getmask(var_in)
				var_in=np.ma.getdata(var_in)
				var_in[masked]=np.nan
			except: pass
			# creat a 1-NA mask
			red_mask = country_mask[np.ix_(list(lats),list(lons))]
			red_mask[red_mask>0]=1
			red_mask[red_mask==0]=np.nan
			country_data=var_in*red_mask

			# handle time information
			time=nc_in.variables['time'][:]
			datevar = []
			# if specified units and calendar
			if time_units!=None and time_calendar!=None:
				datevar.append(num2date(time,units = time_units,calendar= time_calendar))
			# if no specification
			if time_units==None and time_calendar==None:
				time_unit=nc_in.variables['time'].units
				try:	
					cal_temps = nc_in.variables['time'].calendar
					datevar.append(num2date(time,units = time_unit,calendar = cal_temps))
				except:
					datevar.append(num2date(time,units = time_unit))
			# create index variable
			year=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])
			month=np.array([int(str(date).split("-")[1])	for date in datevar[0][:]])

			# write zoomed file
			nc_out=Dataset(out_file,"w")
			nc_out.createDimension('time', len(time))
			nc_out.createDimension('lat', ny)
			nc_out.createDimension('lon', nx)
			# lat lon 
 			outVar = nc_out.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat_mask[list(lats)]	;	outVar.setncattr('units','deg south')
 			outVar = nc_out.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon_mask[list(lons)]	;	outVar.setncattr('units','deg east')
 			# time
 			outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=time	;	outVar.setncatts({k: nc_in.variables['time'].getncattr(k) for k in nc_in.variables['time'].ncattrs()})
 			outVar = nc_out.createVariable('year', 'f', ('time',)) ; outVar[:]=year
 			outVar = nc_out.createVariable('month', 'f', ('time',)) ; outVar[:]=month
 			# data
 			outVar = nc_out.createVariable(var_name, 'f', ('time','lat','lon',),fill_value=np.nan)
 			for k in nc_in.variables[var_name].ncattrs():
 				if k!='_FillValue':outVar.setncatts({k:nc_in.variables[var_name].getncattr(k)})
 			outVar.setncattr('grid',grid)
			outVar[:] = country_data[:,:,:]
			# close the output file
			nc_out.close()
			nc_in.close()
			print out_file

			# store in dictionary
			tmp = self._data
			for meta_info in meta_data:
				if meta_info not in tmp: tmp[meta_info]={}
				tmp=tmp[meta_info]

			tmp['data']=country_data
			tmp['lat']=lat_mask[list(lats)]
			tmp['lon']=lon_mask[list(lons)]
			tmp['year']=year
			tmp['month']=month
			tmp['time']=time
			tmp['grid']=grid

			self._meta.append(meta_data)


	def average(self,mask_style='lat_weighted',meta_data=[]):
		'''
		compute countrywide averages for all loaded datasets
		mask_style: str: weighting used to compute countrywide averages
		meta_data: list of strs: meta_data restrictions. data files for which the meta_data is as specified, the average is computed 
		'''

		# preprare dictionary structure as for self._data
		if '_transient' not in dir(self): self._transient={}
		if os.path.isdir(self._working_directory+'/country_average')==False:os.system('mkdir '+self._working_directory+'/country_average')

		for meta_list in self._meta:
			tmp_in = self._data
			tmp_out = self._transient
			compute=True

			for i in range(len(meta_list)):
				meta_info=meta_list[i]
				if i>=len(meta_data):
					if meta_info not in tmp_out: tmp_out[meta_info]={}
					tmp_in=tmp_in[meta_info]
					tmp_out=tmp_out[meta_info]

				# in case of meta_data restrictions
				if i<len(meta_data):
					if meta_info==meta_data[i] or meta_data[i]==None:
						if meta_info not in tmp_out: tmp_out[meta_info]={}
						tmp_in=tmp_in[meta_info]
						tmp_out=tmp_out[meta_info]
						continue
					else: 
						compute=False
						break
					break

			if compute==True:
				# check if file has been saved
				out_file=self._working_directory+'/country_average/country_mean_'+'_'.join(meta_list)+'_'+mask_style+'.csv'
				if os.path.isfile(out_file):
					tmp_out['time']=np.array(pd.read_csv(out_file,sep=';')['time'])
					tmp_out['month']=np.array(pd.read_csv(out_file,sep=';')['month'])
					tmp_out['year']=np.array(pd.read_csv(out_file,sep=';')['year'])
					tmp_out[mask_style]=np.array(pd.read_csv(out_file,sep=';')['_'.join(meta_list)])	

				# if not compute
				else:
					# load input data
					var_in=tmp_in['data'][:,:,:]			
					try:	# handle masked array
						masked=np.ma.getmask(var_in)
						var_in=np.ma.getdata(var_in)
						var_in[masked]=np.nan
					except: pass

					# get mask
					mask=self._masks[tmp_in['grid']][mask_style]
					lat_mask=self._masks[tmp_in['grid']]['lat_mask']
					lon_mask=self._masks[tmp_in['grid']]['lon_mask']

					# find relevant area (as rectangle) and check whether lon and lat are correct (on same grid differences in lat decreasing or increasing could arise)
					lat_mean=np.mean(mask,1)
					lats=np.where(lat_mean!=0)
					if lat_mask[lats][0]!=tmp_in['lat'][0]:
						var_in=var_in[:,:,::-1]
						if lat_mask[lats][0]!=tmp_in['lat'][-1]:
							print 'problem with lat' ; return('error')

					lon_mean=np.mean(mask,0)
					lons=np.where(lon_mean!=0)
					if lon_mask[lons][0]!=tmp_in['lon'][0]:
						var_in=var_in[:,:,::-1]
						if lon_mask[lons][0]!=tmp_in['lon'][-1]:
							print 'problem with lon' ; return('error')

					# zoom mask to relevant area
					mask=mask[np.ix_(list(lats[0]),list(lons[0]))]
					country_area=np.where(mask>0)

					tmp_out['time']=tmp_in['time']
					tmp_out['month']=tmp_in['month']
					tmp_out['year']=tmp_in['year']
					tmp_out[mask_style]=tmp_in['time'].copy()*np.nan
					for i in range(len(tmp_in['time'])):
						var_of_area=var_in[i,:,:][country_area]
						# NA handling: sum(mask*var)/sum(mask) for the area where var is not NA
						not_missing_in_var=np.where(np.isfinite(var_of_area))[0]	# np.where()[0] because of array([],)
						if len(not_missing_in_var)>0:
							tmp_out[mask_style][i]=sum(mask[country_area][not_missing_in_var]*var_of_area[not_missing_in_var])/sum(mask[country_area][not_missing_in_var])
			
					# save as csv 
					country_mean_csv = pd.DataFrame(index=range(len(tmp_in['time'])))
					country_mean_csv['time']=tmp_in['time']
					country_mean_csv['month']=tmp_in['month']
					country_mean_csv['year']=tmp_in['year']
					country_mean_csv['_'.join(meta_list)]=tmp_out[mask_style]
					country_mean_csv.to_csv(out_file,na_rep='NaN',sep=';',index_label='index')




	def period_averages(self,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]},meta_data=[]):
		'''
		computes time averages for each grid-cell for a given period
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		meta_data: list of strs: meta_data restrictions. data files for which the meta_data is as specified, the average is computed 
		'''

		# preprare dictionary structure as for self._data
		for meta_list in self._meta:
			tmp = self._data
			compute=True

			for i in range(len(meta_list)):
				meta_info=meta_list[i]
				if i>=len(meta_data):
					tmp=tmp[meta_info]

				# in case of meta_data restrictions
				if i<len(meta_data):
					if meta_info==meta_data[i] or meta_data[i]==None:
						tmp=tmp[meta_info]
						continue
					else: 
						compute=False
						break
					break

			if compute==True:
				tmp['period']={}
				tmp['period_meta']=periods

				for period_name in periods.keys():	
					tmp['period'][period_name]=np.mean(np.ma.masked_invalid(tmp['data']),axis=0)

				for period_name in periods.keys():
					if period_name!='ref' and 'ref' in periods.keys():
						tmp['period'][period_name+'-'+'ref']=tmp['period'][period_name]-tmp['period']['ref']


	def plot_transient(self,meta_data,mask_style='lat_weighted',running_mean=1,ax=None,out_file=None,title=None,ylabel=None,show=True):
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

		tmp = self._transient
		for meta_info in meta_data:
			tmp=tmp[meta_info]

		if ax==None:
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

		if int(np.mean(np.diff(tmp['year'],1)))!=1:print 'not yearly data! please consider this for the running mean'

		ax.plot(tmp['time'],pd.rolling_mean(tmp[mask_style],running_mean),linestyle='-',label=mask_style)

		ax.set_xticks(tmp['time'][range(0,len(tmp['time']),240)]) 
		ax.set_xticklabels(tmp['year'][range(0,len(tmp['time']),240)])

		if ylabel==None:ylabel=meta_data[0]
		ax.set_ylabel(ylabel)
		if title==None:title=' '.join(meta_data)
		ax.set_title(title)
		ax.legend(loc='best')
		
		if out_file==None and show==True:plt.show()
		if out_file!=None:plt.savefig(out_file)




	def plot_map(self,meta_data,source='_data',period=None,time=None,color_bar=True,color_label=None,color_palette=plt.cm.plasma,color_range=None,grey_area=None,limits=None,ax=None,out_file=None,title=None,show=True):
		'''
		plot maps of data. 
		meta_data: list of strs: meta information required to acces data
		source: str: default='_data'. if masks are to be plotted, specify source='_masks'
		period: str: if  the averag over a period is to be plotted specify the period name
		time: int: index in time axis of data (to be plotted)
		color_bar: logical: if True, color-scale is plotted besides the map
		color_label: str: label of the color-scale
		color_palette: plt.cm. object: colors used
		color_range: [float,float]: minimal and maximal value on color-scale
		limits: [lon_min,lon_max,lat_min,lat_max]: extend of the map
		ax: subplot: subplot on which the map will be plotted
		out_file: str: location where the plot is saved
		title: str: title of the plot
		show: logical: show the subplot?
		'''
		if source=='_masks':
			tmp = self._masks[meta_data[0]]
			to_plot=tmp[meta_data[1]][meta_data[2]]
			lat=tmp['lat_mask']
			lon=tmp['lon_mask']				

			if color_label==None:color_label='importance of grid-cell\nfor countrywide average'
			if title==None:title=' '.join(meta_data)

		if source=='_data':
			tmp = self._data

			for meta_info in meta_data:
				tmp=tmp[meta_info]

			if period==None:
				if time==None:
					time=int(len(tmp['time'])/2)
					print 'no time specified. '+str(int(tmp['month'][time]))+'/'+str(int(tmp['year'][time]))+' selected'
					to_plot=tmp['data'][time,:,:]
					if title==None:title=' '.join(meta_data+[str(int(tmp['month'][time])),'/',str(int(tmp['year'][time]))])
			else:
				to_plot=tmp['period'][period]
				if title==None:title=' '.join(meta_data+[period])

			lat=tmp['lat']
			lon=tmp['lon']
			if color_label==None:color_label=meta_data[0]


		if ax==None:
			fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))		

		# handle 0 to 360 lon
		if max(lon)>180:
			problem_start=np.where(lon>180)[0][0]
			new_order=np.array(range(problem_start,len(lon))+range(0,problem_start))
			to_plot=to_plot[:,new_order]
			lon=lon[new_order]
			lon[lon>180]-=360

		# handle limits
		half_lon_step=abs(np.diff(lon,1)[0]/2)
		half_lat_step=abs(np.diff(lat,1)[0]/2)
		if limits==None:
			limits=[np.min(lon)-half_lon_step,np.max(lon)+half_lon_step,np.min(lat)-half_lat_step,np.max(lat)+half_lat_step]
		if limits!=None:
			lon_select=np.where((lon>=limits[0])	&	(lon<=limits[1]))[0]
			lat_select=np.where((lat>=limits[2])	&	(lat<=limits[3]))[0]
			to_plot=to_plot[lat_select,:][:,lon_select]
		# correct limits if necessary
		extent=[np.min(lon)-half_lon_step,np.max(lon)+half_lon_step,np.min(lat)-half_lat_step,np.max(lat)+half_lat_step]

		if limits[0]<extent[0]:limits[0]=extent[0]; comment='limits changed to extend of data! new limits =',limits
		if limits[1]>extent[1]:limits[1]=extent[1]; comment='limits changed to extend of data! new limits =',limits
		if limits[2]<extent[2]:limits[2]=extent[2]; comment='limits changed to extend of data! new limits =',limits
		if limits[3]>extent[3]:limits[3]=extent[3]; comment='limits changed to extend of data! new limits =',limits
		if 'comment' in dir(): print comment

		m = Basemap(ax=ax,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
		m.drawmapboundary(fill_color='1.')

		# imshow does not support decreasing lat or lon
		to_plot=np.ma.masked_invalid(to_plot)
		if lat[0]>lat[1]:to_plot=to_plot[::-1,:]
		if lon[0]>lon[1]:to_plot=to_plot[:,::-1]



		# get color_range
		if color_range==None:
			color_range=[np.min(to_plot[np.isfinite(to_plot)]),np.max(to_plot[np.isfinite(to_plot)])]

		im = m.imshow(to_plot,cmap=color_palette,vmin=color_range[0],vmax=color_range[1],interpolation='none',extent=extent)

		# mask some grid-cells
		if grey_area!=None:
			to_plot=np.ma.masked_invalid(grey_area.copy())
			if lat[0]>lat[1]:to_plot=to_plot[::-1,:]
			if lon[0]>lon[1]:to_plot=to_plot[:,::-1]
			im2 = m.imshow(to_plot,cmap=plt.cm.Greys,vmin=0,vmax=1,interpolation='none',extent=extent)

		# show coastlines and borders
		m.drawcoastlines()
		m.drawstates()
		m.drawcountries()

		# add colorbar
		if color_bar==True:
			cb = m.colorbar(im,'right', size="5%", pad="2%")
			tick_locator = ticker.MaxNLocator(nbins=5)
			cb.locator = tick_locator
			cb.update_ticks()
			cb.set_label(color_label, rotation=90)

		ax.set_title(title)
		ax.legend(loc='best')
		
		if out_file==None and show==True:plt.show()
		if out_file!=None:plt.savefig(out_file)
		return(im)









