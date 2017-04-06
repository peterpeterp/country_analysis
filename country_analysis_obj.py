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

from matplotlib import rc
rc('text', usetex=True)

def plot_map(to_plot,lat,lon,color_bar=True,color_label='',color_palette=plt.cm.plasma,color_range=None,grey_area=None,limits=None,ax=None,out_file=None,title=''):
	# this actualy plots the map

		if ax==None: show=True
		if ax!=None: show=False

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
		if limits==None:
			half_lon_step=abs(np.diff(lon,1)[0]/2)
			half_lat_step=abs(np.diff(lat,1)[0]/2)
			relevant_lats=lat[np.where(np.isfinite(to_plot))[0]]
			relevant_lons=lon[np.where(np.isfinite(to_plot))[1]]
			limits=[np.min(relevant_lons)-half_lon_step,np.max(relevant_lons)+half_lon_step,np.min(relevant_lats)-half_lat_step,np.max(relevant_lats)+half_lat_step]

		m = Basemap(ax=ax,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
		m.drawmapboundary(fill_color='1.')

		# get color_range
		if color_range==None:
			color_range=[np.min(to_plot[np.isfinite(to_plot)]),np.max(to_plot[np.isfinite(to_plot)])]

		x,y=lon.copy(),lat.copy()
		x-=np.diff(x,1)[0]/2.
		y-=np.diff(y,1)[0]/2.
		x=np.append(x,[x[-1]+np.diff(x,1)[0]])
		y=np.append(y,[y[-1]+np.diff(y,1)[0]])
		x,y=np.meshgrid(x,y)
		im = m.pcolormesh(x,y,to_plot,cmap=color_palette,vmin=color_range[0],vmax=color_range[1])


		# mask some grid-cells
		if grey_area!=None:
			to_plot=np.ma.masked_invalid(grey_area.copy())
			if y[0]>y[1]:to_plot=to_plot[::-1,:]
			if x[0]>x[1]:to_plot=to_plot[:,::-1]
			im2 = m.pcolormesh(x,y,to_plot,cmap=plt.cm.Greys,vmin=0,vmax=1)

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

class country_analysis(object):

	def __init__(self,iso,working_directory):
		self._iso=iso
		self._working_directory=working_directory+iso

		self._masks={}
		self._regions=[iso]

		self._data={}
		self._meta=[]

		self._DATA=[]
		self._META={}

		self._file_info={'masks':{},'raw':{},'country_average':{}}

		if os.path.isdir(self._working_directory)==False:os.system('mkdir '+self._working_directory)
		if os.path.isdir(self._working_directory+'/masks')==False:os.system('mkdir '+self._working_directory+'/masks')
		if os.path.isdir(self._working_directory+'/plots')==False:os.system('mkdir '+self._working_directory+'/plots')
		if os.path.isdir(self._working_directory+'/raw')==False:os.system('mkdir '+self._working_directory+'/raw')
		if os.path.isdir(self._working_directory+'/country_average')==False:os.system('mkdir '+self._working_directory+'/country_average')

	class data(object):
		def __init__(SELF,tags,outer_self,raw_file,original_var_name):
			SELF.index=len(outer_self._DATA)
			SELF.raw_file=outer_self._working_directory+'/raw/'+raw_file.split('raw/')[-1]
			SELF.original_var_name=original_var_name
			SELF.tags=tags
			if 'var' in tags.keys():	SELF.var=tags['var']
			if 'type' in tags.keys():	SELF.type=tags['type']
			if 'scenario' in tags.keys():	SELF.scenario=tags['scenario']
			if 'model' in tags.keys():	SELF.model=tags['model']

			SELF.average={}

			SELF.all_tags=[]

			for key in tags.keys():
				SELF.all_tags.append(tags[key])
				if key not in outer_self._META.keys(): outer_self._META[key]={}
				if tags[key] not in outer_self._META[key].keys(): outer_self._META[key][tags[key]]=[]
				outer_self._META[key][tags[key]].append(SELF.index)

			SELF.name='_'.join(SELF.all_tags)



		def get_time_stamp(SELF):
			SELF.time_stamp=np.array([int(SELF.year[i])+int(SELF.month[i])/100. for i in range(len(SELF.year))])

		def display_map(SELF,period=None,time=None,color_bar=True,color_label=None,color_palette=None,color_range=None,grey_area=None,limits=None,ax=None,out_file=None,title=None):
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

			# this prepares plotting

			if period==None:
				if time==None:
					time=int(len(SELF.time)/2)
					print 'no time specified. '+str(int(SELF.month[time]))+'/'+str(int(SELF.year[time]))+' selected'
					to_plot=SELF.raw[time,:,:].copy()
					if title==None:title='_'.join([SELF.name]+[str(int(SELF.month[time])),'/',str(int(SELF.year[time]))])
			else:
				to_plot=SELF.period[period].copy()
				if title==None:title=SELF.name+'_'+period

			lat=SELF.lat.copy()
			lon=SELF.lon.copy()
			if color_label==None:color_label=SELF.var

			if color_palette==None:
				if SELF.var=='pr':		color_palette=plt.cm.YlGnBu
				elif SELF.var=='tas':	color_palette=plt.cm.YlOrBr
				elif SELF.var==	'SPEI':	color_palette=plt.cm.plasma
				else:					color_palette=plt.cm.plasma

			im=plot_map(to_plot,lat,lon,color_bar=color_bar,color_label=color_label,color_palette=color_palette,color_range=color_range,grey_area=grey_area,limits=limits,ax=ax,out_file=out_file,title=title)
			return(im)

		def plot_transient(SELF,mask_style=None,region=None,running_mean=1,ax=None,out_file=None,title=None,ylabel=None,label=''):
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

			if ax!=None:
				show=False

			if ax==None:
				show=True
				fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

			if mask_style!=None:
				mask_styles=[mask_style]
			if mask_style==None:
				mask_styles=SELF.average.keys()

			if region!=None:
				regions=[region]
			if region==None:
				regions=SELF.average[mask_styles[0]].keys()

			for mask_style in mask_styles:
				for region in regions:
					ax.plot(SELF.time,pd.rolling_mean(SELF.average[mask_style][region],running_mean),linestyle='-',label=region+' '+mask_style.replace('_',' '))

			ax.set_xticks(SELF.time[range(0,len(SELF.time),240)]) 
			ax.set_xticklabels(SELF.year[range(0,len(SELF.time),240)])

			if ylabel==None:ylabel=SELF.var.replace('_',' ')
			ax.set_ylabel(ylabel)
			if title==None:title=SELF.name.replace('_',' ')
			ax.set_title(title)
			
			if show==True:ax.legend(loc='best')
			if out_file==None and show==True:plt.show()
			if out_file!=None:plt.savefig(out_file)

			if int(np.mean(np.diff(SELF.year,1)))!=1:
				return 'not yearly data! please consider this for the running mean'

		def plot_annual_cycle(SELF,mask_style='lat_weighted',region=None,period=None,ax=None,out_file=None,title=None,ylabel=None,label=''):
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

			if ax!=None:
				show=False

			if ax==None:
				show=True
				fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,4))

			if period==None:
				period=[min(SELF.year),max(SELF.year)]

			if region!=None:
				regions=[region]
			if region==None:
				regions=SELF.average[mask_style].keys()

			for region in regions:
				annual_cycle=[]
				for mon in range(1,13):
					relevant_time=np.where((SELF.year>=period[0]) & (SELF.year<period[1]) & (SELF.month==mon))[0]
					annual_cycle.append(np.nanmean(SELF.average[mask_style][region][relevant_time]))


				ax.plot(range(0,12),annual_cycle,linestyle='-',label=label)

			ax.set_xticks(range(0,12)) 
			ax.set_xticklabels(['JAN','FEB','MAR','APR','MAI','JUN','JUL','AUG','SEP','OCT','NOV','DEC'])

			if ylabel==None:ylabel=SELF.var.replace('_',' ')
			ax.set_ylabel(ylabel)
			if title==None:title=SELF.name.replace('_',' ')
			ax.set_title(title)
			
			if show==True:ax.legend(loc='best')
			if out_file==None and show==True:plt.show()
			if out_file!=None:plt.savefig(out_file)



	def display_mask(self,grid=None,mask_style=None):
		if grid==None:
			print 'Please select a grid:'
			for grid in self._masks.keys():
				print grid
			return None
		if mask_style==None:
			print 'Please select a mask-style:'
			for key in self._masks[grid].keys():
				if key not in ['lat_mask','lon_mask']:
					print key
			return None

		else:
			toplo=self._masks[grid][mask_style][self._iso]
			lat=self._masks[grid]['lat_mask'].copy()
			lon=self._masks[grid]['lon_mask'].copy()

			plot_map(toplo,lat,lon,title=grid+' '+mask_style)			

	def display(self,selection=None):
		if selection==None:
			selection=self._DATA
		for data in selection:
			print data.index,data.name,min(data.year),max(data.year)

	def reorder_data():
		for data in self._DATA:
			print 66

	def unit_conversions(self):
		for data in self._DATA:
			if data.var=='tas':
				if np.nanmean(data.raw)>100:
					data.raw-=273.15
			if data.var=='pr':
				if np.nanmax(data.raw)<10:
					data.raw*=86400

	def selection(self,filters):
		selection=[]
		for data in self._DATA:
			selected=True
			for key in filters:
				if key not in data.all_tags:
					selected=False
			if selected:
				selection.append(data)

		self.display(selection)
		return selection

	def ensemble(self,filters):
		ensemble={}
		for data in self._DATA:
			selected=True
			for key in filters:
				if key not in data.all_tags:
					selected=False
			if selected:
				ensemble[data.model]=data
				print data.model+': '+data.name

		return ensemble

	def prepare_for_download(self):
		sys.path.append(self._working_directory)
		sys.path.append('../')

		file_info={'masks':{},'raw':{},'country_average':{}}

		for grid in self._masks.keys():
			for mask_style in self._masks[grid].keys():
				if mask_style not in ['lat_mask','lon_mask']:
					file_info['masks'][self._working_directory+'/masks/'+self._iso+'_'+grid+'_'+mask_style+'.nc4']=[grid,mask_style]

		for data in self._DATA:
			file_info['raw'][data.raw_file]={'tags':data.tags,'var_name':data.original_var_name}
			for mask_style in data.average.keys():
				file_info['country_average'][data.average[mask_style]['out_file']]={'mask_style':mask_style,'name':data.name}

		#print file_info

		output = open(self._working_directory+'/'+self._iso+'_file_info.pkl', 'wb')
		pickle.dump(file_info, output)	;	output.close()

		os.chdir(self._working_directory)
		os.chdir('../')
		os.system('tar -zcf '+self._working_directory+'.tar.gz '+self._iso)


	def load_from_tar(self,path):
		os.system('tar -zxf '+self._working_directory+'.tar.gz -C '+self._working_directory.replace(self._iso,''))

		pkl_file = open(self._working_directory+'/'+self._iso+'_file_info.pkl', 'rb')
		self._file_info = pickle.load(pkl_file)	;	pkl_file.close()  

		for file in self._file_info['masks'].keys():
			meta=self._file_info['masks'][file]
			grid,mask_style = meta[0],meta[1]
			file_new=self._working_directory+'/masks'+file.split('masks')[-1]
			if grid not in self._masks.keys():
				self._masks[grid]={}
			if mask_style not in self._masks[grid].keys():
				self._masks[grid][mask_style]={}
			self.load_masks(file_new,grid,mask_style)


		for file in self._file_info['raw'].keys():
			tags=self._file_info['raw'][file]['tags']
			var_name=self._file_info['raw'][file]['var_name']

			file_new=self._working_directory+'/raw'+file.split('raw')[-1]

			#print file_new
			nc_out=Dataset(file_new,"r")

			new_data=self.data(tags=tags,outer_self=self,raw_file=file,original_var_name=var_name)

			new_data.raw=nc_out.variables[var_name][:,:,:]
			new_data.grid=nc_out.variables[var_name].getncattr('grid')
			new_data.lat=nc_out.variables['lat'][:]
			new_data.lon=nc_out.variables['lon'][:]
			new_data.time=nc_out.variables['time'][:]
			new_data.year=nc_out.variables['year'][:]
			new_data.month=nc_out.variables['month'][:]
			new_data.get_time_stamp()


			self._DATA.append(new_data)



		for file in self._file_info['country_average'].keys():
			mask_style=self._file_info['country_average'][file]['mask_style']
			name=self._file_info['country_average'][file]['name']

			file_new=self._working_directory+'/country_average'+file.split('country_average')[-1]

			for data in self._DATA:
				if data.name==name:
					table=pd.read_csv(file_new,sep=';')
					for key in table.keys():
						if key not in ['time','year','month','index']:
							if mask_style not in data.average.keys():	data.average[mask_style]={}
							data.average[mask_style][key]=np.array(table[key])




	def identify_grid(self,input_file,lat_name,lon_name):
		# get information about grid of input data
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
		grid=str(ny)+'x'+str(nx)
		nc_in.close()

		return lon,lat,grid,lon_shift	

	def load_masks(self,mask_file,grid,mask_style):
		# load existing mask
		nc_mask=Dataset(mask_file,'r')
		self._masks[grid]['lat_mask'] = nc_mask.variables['lat'][:]  
		self._masks[grid]['lon_mask'] = nc_mask.variables['lon'][:] 

		# get all variables (regions)
		for name in nc_mask.variables.keys():
			if name not in ['lat','lon']:
				self._masks[grid][mask_style][name] = nc_mask.variables[name][:,:]  
				if name not in self._regions: self._regions.append(name)	



	def get_grid_polygons(self,grid,lon,lat,lon_shift):
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

	def regrid_pop_mask(self,grid,lon,lat,shift,pop_mask_file,mask_style):
		# regrid population mask 
		mygrid=open(self._working_directory+'/masks/'+grid+'.txt','w')
		mygrid.write('gridtype=lonlat\nxsize='+str(len(lon))+'\nysize='+str(len(lat))+'\nxfirst='+str(lon[0])+'\nxinc='+str(np.mean(np.diff(lon,1)))+'\nyfirst='+str(lat[0])+'\nyinc='+str(np.mean(np.diff(lat,1))))
		mygrid.close()
		os.system('cdo remapbil,'+self._working_directory+'/masks/'+grid+'.txt '+pop_mask_file+' '+self._working_directory+'/masks/'+mask_style+'_'+grid+'.nc')	
		nc_pop_mask = Dataset(self._working_directory+'/masks/'+mask_style+'_'+grid+'.nc')
		pop_mask = np.array(nc_pop_mask.variables['mask'][:,:]).squeeze()
		pop_mask = np.roll(pop_mask,shift,axis=1)

		return pop_mask

	def grid_polygon_overlap(self,grid,lon,lat,grid_polygons,country_polygons,shift,mask_style,ext_poly,name,pop_mask=None):
		nx = len(lon)	;	ny = len(lat)

		overlap = np.zeros((ny,nx))
		for i in range(nx):
			for j in range(ny):
				# check gridcell is relevant
				if grid_polygons[i,j].intersects(ext_poly):
					# get fraction of grid-cell covered by polygon
					intersect = grid_polygons[i,j].intersection(country_polygons).area/grid_polygons[i,j].area*country_polygons.area
					if pop_mask!=None:
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

	def create_mask_country(self,input_file,var_name,shape_file,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon'):
		'''
		create country mask
		input_file: str: location of example input data (required for the identification of the grid)
		var_name: str: variable name of input file
		shape_file: str: location of the shape_file used to identify country borders
		mask_style: str: name under which the mask will be stored (important for further analysis)
		pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
		'''

		lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

		if grid not in self._masks.keys():
			self._masks[grid]={}
		if mask_style not in self._masks[grid].keys():
			self._masks[grid][mask_style]={}

		mask_file=self._working_directory+'/masks/'+self._iso+'_'+grid+'_'+mask_style+'.nc4'
		self._file_info['masks'][mask_file]=[grid,mask_style]

		if os.path.isfile(mask_file) and overwrite==False:
 			self.load_masks(mask_file,grid,mask_style)

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
			xmin, xmax, ymin, ymax = country_polygons.bounds
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
			nc_mask.close()

		

	def create_mask_admin(self,input_file,var_name,shape_file,mask_style='lat_weighted',pop_mask_file='',overwrite=False,lat_name='lat',lon_name='lon'):
		'''
		create country mask
		input_file: str: location of example input data (required for the identification of the grid)
		var_name: str: variable name of input file
		shape_file: str: location of the shape_file used to identify country borders
		mask_style: str: name under which the mask will be stored (important for further analysis)
		pop_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
		'''

		lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

		if grid not in self._masks.keys():
			self._masks[grid]={}
		if mask_style not in self._masks[grid].keys():
			self._masks[grid][mask_style]={}

		mask_file=self._working_directory+'/masks/'+self._iso+'_admin_'+grid+'_'+mask_style+'.nc4'
		self._file_info['masks'][mask_file]=[grid,mask_style]

		if os.path.isfile(mask_file) and overwrite==False:
 			self.load_masks(mask_file,grid,mask_style)

		else:
			grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)
			
			# load shape file
			m = Basemap()
			m.readshapefile(shape_file, 'admin', drawbounds=False)

			# collect all shapes of region
			region_polygons={}
			for shape, region in zip(m.admin, m.admin_info):
				region = {k.lower():v for k,v in region.items()}	
				name = region['name_1']
				self._regions.append(name)
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

			for name in region_polygons.keys():
				region_shape = region_polygons[name]
				self.grid_polygon_overlap(grid,lon, lat, grid_polygons, region_shape, shift, mask_style, ext_poly, name, pop_mask)
				outVar = nc_mask.createVariable(name, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][name][:,:]
			
			nc_mask.close()

	def understand_time_format(self,nc_in,time_units=None,time_calendar=None):
		time=nc_in.variables['time'][:]
		datevar = []
		# if specified units and calendar
		if time_units!=None and time_calendar!=None:
			datevar.append(num2date(time,units = time_units,calendar= time_calendar))
		# if no specification
		time[time<0]=0
		if time_units==None and time_calendar==None:
			time_unit=nc_in.variables['time'].units
			if True:	
				cal_temps = nc_in.variables['time'].calendar
				datevar.append(num2date(time,units = time_unit,calendar = cal_temps))
			# except:
			# 	datevar.append(num2date(time,units = time_unit))
		# create index variable
		year=np.array([int(str(date).split("-")[0])	for date in datevar[0][:]])
		month=np.array([int(str(date).split("-")[1])	for date in datevar[0][:]])

		return(time,year,month)

	def country_zoom(self,input_file,var_name,tags=None,mask_style='lat_weighted',time_units=None,time_calendar=None,lat_name='lat',lon_name='lon',overwrite=False,additional_tag=''):
		'''
		zoom input_file to area relevant for the country
		input_file: type str: file to be processed
		out_file: type str: filepath where output is going to be stored
		var: type str: variable name in netcdf file
		iso: type str: iso3c of country
		mask_path: type str: path to where the masks are stored
		'''

		out_file=self._working_directory+'/raw/'+input_file.split('/')[-1].replace('.nc',additional_tag+'_'+self._iso+'.nc')

		if os.path.isfile(out_file) and overwrite==False:
			nc_out=Dataset(out_file,"r")

			new_data=self.data(tags=tags,outer_self=self,raw_file=out_file,original_var_name=var_name)

			print var_name
			print tags
			print nc_out.variables.keys()
			new_data.raw=nc_out.variables[var_name][:,:,:]
			new_data.grid=nc_out.variables[var_name].getncattr('grid')
			new_data.lat=nc_out.variables['lat'][:]
			new_data.lon=nc_out.variables['lon'][:]
			new_data.time=nc_out.variables['time'][:]
			new_data.year=nc_out.variables['year'][:]
			new_data.month=nc_out.variables['month'][:]
			new_data.get_time_stamp()

			self._DATA.append(new_data)

			nc_out.close()


		else:
			# open file to get information
			print input_file
			lon_in,lat_in,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)
			nc_in=Dataset(input_file,"r")

			country_mask=self._masks[grid][mask_style][self._iso]
			lat_mask=self._masks[grid]['lat_mask']
			lon_mask=self._masks[grid]['lon_mask']
			
			# check whether lon and lat are similarly defined in mask and input_file
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
			time,year,month=self.understand_time_format(nc_in)


			# write zoomed file
			nc_out=Dataset(out_file,"w")
			nc_out.createDimension('time', len(time))
			nc_out.createDimension('bnds', 2)
			nc_out.createDimension('lat', ny)
			nc_out.createDimension('lon', nx)
			# lat lon 
 			outVar = nc_out.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat_mask[list(lats)]	;	outVar.setncattr('units','deg south')
 			outVar = nc_out.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon_mask[list(lons)]	;	outVar.setncattr('units','deg east')
 			# time
 			outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=time	;	outVar.setncatts({k: nc_in.variables['time'].getncattr(k) for k in nc_in.variables['time'].ncattrs()})

 			# still not sure how to handle these time bounds!
 			try:
 				outVar = nc_out.createVariable('time_bnds', 'f', ('time','bnds')) ; outVar[:]=nc_in.variables['time_bnds'][:,:]	;	outVar.setncatts({k: nc_in.variables['time_bnds'].getncattr(k) for k in nc_in.variables['time_bnds'].ncattrs()})
 			except:
 				print 'issue with time_bnds '+out_file

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

			new_data=self.data(tags=tags,outer_self=self,raw_file=out_file,original_var_name=var_name)
			
			new_data.raw=country_data
			new_data.grid=grid
			new_data.lat=lat_mask[list(lats)]
			new_data.lon=lon_mask[list(lons)]
			new_data.time=time
			new_data.year=year
			new_data.month=month
			new_data.get_time_stamp()

			self._DATA.append(new_data)


	def hist_merge(self):

		for data in self._DATA:
			if hasattr(data,'model'):
				data__=[]
				for dd in self._DATA:
					if hasattr(dd,'model'):
						if dd.model==data.model and dd.var==data.var and dd.type==data.type:	data__.append(dd)

				for hist in data__:
					if hist.scenario.lower() in ['hist','historical']:
						delete_hist=False
						print '------- merging -----------'
						print hist.model,hist.var,hist.scenario
						print data__

						for ddd in data__:	
							if ddd.scenario.lower() not in ['hist','historical']:
								os.system('cdo -O mergetime '+hist.raw_file+' '+ddd.raw_file+' '+ddd.raw_file.replace('.nc','_tmp.nc'))

								nc_in=Dataset(ddd.raw_file.replace('.nc','_tmp.nc'))
								time,year,month=self.understand_time_format(nc_in)

								os.system('rm '+ddd.raw_file.replace('.nc','_merged.nc'))
								nc_out=Dataset(ddd.raw_file.replace('.nc','_merged.nc'),"w")
								nc_out.createDimension('time', len(nc_in.variables['time'][:]))
								nc_out.createDimension('lat', len(ddd.lat))
								nc_out.createDimension('lon', len(ddd.lon))
								# lat lon 
					 			outVar = nc_out.createVariable('lat', 'f', ('lat',)) ; outVar[:]=ddd.lat	;	outVar.setncattr('units','deg south')
					 			outVar = nc_out.createVariable('lon', 'f', ('lon',)) ; outVar[:]=ddd.lon	;	outVar.setncattr('units','deg east')
					 			# time
					 			outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=time	;	outVar.setncatts({k: nc_in.variables['time'].getncattr(k) for k in nc_in.variables['time'].ncattrs()})
					 			outVar = nc_out.createVariable('year', 'f', ('time',)) ; outVar[:]=year
					 			outVar = nc_out.createVariable('month', 'f', ('time',)) ; outVar[:]=month
					 			outVar = nc_out.createVariable(ddd.original_var_name, 'f', ('time','lat','lon',),fill_value=np.nan)
					 			for k in nc_in.variables[ddd.original_var_name].ncattrs():
					 				if k!='_FillValue':outVar.setncatts({k:nc_in.variables[ddd.original_var_name].getncattr(k)})
					 			outVar.setncattr('grid',ddd.grid)
								outVar[:] = nc_in.variables[ddd.original_var_name][:,:,:]

								nc_in.close()
								os.system('rm '+ddd.raw_file.replace('.nc','_tmp.nc'))
							

								print ddd.model,ddd.var,ddd.scenario
								

								print ddd.original_var_name
								ddd.raw=nc_out.variables[ddd.original_var_name][:,:,:]
								ddd.time=nc_out.variables['time'][:]
								ddd.year=nc_out.variables['year'][:]
								ddd.month=nc_out.variables['month'][:]
								ddd.get_time_stamp()
								ddd.raw_file=ddd.raw_file.replace('.nc','_merged.nc')

								nc_out.close()

								delete_hist=True
				
						if delete_hist:	
							self._DATA.remove(hist)


	def average(self,mask_style='lat_weighted',filters={},overwrite=False):
		'''
		compute countrywide averages for all loaded datasets
		mask_style: str: weighting used to compute countrywide averages
		meta_data: list of strs: meta_data restrictions. data files for which the meta_data is as specified, the average is computed 
		'''

		for data in self._DATA:
			compute=True
			for key in filters.keys():
				if filters[key] not in data.all_tags:
					compute=False

			if compute:
				# check if file has been saved
				out_file=self._working_directory+'/country_average/country_mean_'+'_'.join(data.all_tags)+'_'+mask_style+'.csv'
				if mask_style not in data.average.keys():	data.average[mask_style]={}
				data.average[mask_style]['out_file']=out_file

				if os.path.isfile(out_file) and overwrite==False:
					table=pd.read_csv(out_file,sep=';')
					for key in table.keys():
						if key not in ['time','year','month']:
							data.average[mask_style][key]=np.array(table[key])

				else:
					# prepare table
					country_mean_csv = pd.DataFrame(index=range(len(data.time)))
					country_mean_csv['time']=data.time
					country_mean_csv['month']=data.month
					country_mean_csv['year']=data.year

					# load input data
					var_in=data.raw.copy()			
					try:	# handle masked array
						masked=np.ma.getmask(var_in)
						var_in=np.ma.getdata(var_in)
						var_in[masked]=np.nan
					except: pass

					# find relevant area (as rectangle) and check whether lon and lat are correct (on same grid differences in lat decreasing or increasing could arise)
					mask=self._masks[data.grid][mask_style][self._iso]
					lat_mask=self._masks[data.grid]['lat_mask']
					lon_mask=self._masks[data.grid]['lon_mask']

					lat_mean=np.mean(mask,1)
					lats=np.where(lat_mean!=0)
					if lat_mask[lats][0]!=data.lat[0]:
						var_in=var_in[:,:,::-1]
						if lat_mask[lats][0]!=data.lat[-1]:
							print 'problem with lat' ; return('error')

					lon_mean=np.mean(mask,0)
					lons=np.where(lon_mean!=0)
					if lon_mask[lons][0]!=data.lon[0]:
						var_in=var_in[:,:,::-1]
						if lon_mask[lons][0]!=data.lon[-1]:
							print 'problem with lon' ; return('error')

					# get mask
					for name in self._masks[data.grid][mask_style].keys():
						mask=self._masks[data.grid][mask_style][name]

						# zoom mask to relevant area
						mask=mask[np.ix_(list(lats[0]),list(lons[0]))]
						country_area=np.where(mask>0)

						data.average[mask_style][name]=data.time.copy()*np.nan
						for i in range(len(data.time)):
							var_of_area=var_in[i,:,:][country_area]
							# NA handling: sum(mask*var)/sum(mask) for the area where var is not NA
							not_missing_in_var=np.where(np.isfinite(var_of_area))[0]	# np.where()[0] because of array([],)
							if len(not_missing_in_var)>0:
								data.average[mask_style][name][i]=sum(mask[country_area][not_missing_in_var]*var_of_area[not_missing_in_var])/sum(mask[country_area][not_missing_in_var])
				
						country_mean_csv[name.encode('utf-8')]=data.average[mask_style][name]

					# save as csv 
					country_mean_csv.to_csv(out_file,na_rep='NaN',sep=';',index_label='index')



	def period_averages(self,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]},filters={}):
		'''
		computes time averages for each grid-cell for a given period
		periods: dict={'period_name':[start_year,end_year], ...}: start and end years of period
		meta_data: list of strs: meta_data restrictions. data files for which the meta_data is as specified, the average is computed 
		'''

		for data in self._DATA:
			compute=True
			for key in filters.keys():
				if filters[key] not in data.all_tags:
					compute=False

			if compute:
				#print data.name

				data.period={}
				data.period_meta=periods

				for period_name in periods.keys():	
					if len(np.where((data.year>=periods[period_name][0]) & (data.year<periods[period_name][1]))[0])>0:
						years_in_period=np.where((data.year>=periods[period_name][0]) & (data.year<periods[period_name][1]))

						data.period[period_name]=np.mean(np.ma.masked_invalid(data.raw[years_in_period,:,:][0,:,:,:]),axis=0)

				for period_name in periods.keys():
					if period_name!='ref' and 'ref' in periods.keys() and period_name in data.period.keys():
						data.period[period_name+'-'+'ref']=data.period[period_name]-data.period['ref']



	def ensemble_mean(self):

		remaining=self._DATA[:]

		for data in remaining:
			if hasattr(data,'model'):
				if data.model != 'ensemble_mean':
					print '----------------'
					ensemble=[]
					time_min,time_max=[],[]
					for dd in self._DATA:

						if hasattr(dd,'model'):
							if data.model != 'ensemble_mean':
								if dd.type==data.type and dd.var==data.var and dd.scenario==data.scenario:	
									ensemble.append(dd)
									remaining.remove(dd)
									time_min.append(min(dd.time_stamp))
									time_max.append(max(dd.time_stamp))


					time_min=max(time_min)
					time_max=min(time_max)

					tags_=ensemble[0].tags.copy()
					tags_['model']='ensemble_mean'

					new_data=self.data(tags=tags_,outer_self=self,raw_file=self._working_directory+'/raw/'+ensemble[0].var+'_'+ensemble[0].type+'_ensemble-mean_'+ensemble[0].scenario+'.nc',original_var_name=ensemble[0].original_var_name)


					command='cdo -O ensmean '
					for member,i in zip(ensemble,range(len(ensemble))):
						os.system('cdo -selyear,1951/2098 '+member.raw_file+' '+self._working_directory+'/raw/'+str(i)+'tmp.nc4')
						command+=self._working_directory+'/raw/'+str(i)+'tmp.nc4 '

					os.system(command+' '+self._working_directory+'/raw/ensmean_tmp.nc4')

					nc_ens=Dataset(self._working_directory+'/raw/ensmean_tmp.nc4',"r")
					new_data.raw=nc_ens.variables[ensemble[0].original_var_name][:,:,:]
					time,year,month=self.understand_time_format(nc_ens)

					new_data.grid=member.grid
					new_data.lat=member.lat
					new_data.lon=member.lon
					new_data.time=time
					new_data.year=year
					new_data.month=month

					self._DATA.append(new_data)

					nc_in=Dataset(ensemble[0].raw_file,"r")
					os.system('rm '+new_data.raw_file)
					nc_out=Dataset(new_data.raw_file,"w")
					nc_out.createDimension('time', len(time))
					nc_out.createDimension('lat', len(new_data.lat))
					nc_out.createDimension('lon', len(new_data.lon))
					# lat lon 
		 			outVar = nc_out.createVariable('lat', 'f', ('lat',)) ; outVar[:]=new_data.lat	;	outVar.setncattr('units','deg south')
		 			outVar = nc_out.createVariable('lon', 'f', ('lon',)) ; outVar[:]=new_data.lon	;	outVar.setncattr('units','deg east')
		 			# time
		 			outVar = nc_out.createVariable('time', 'f', ('time',)) ; outVar[:]=time	;	outVar.setncatts({k: nc_in.variables['time'].getncattr(k) for k in nc_in.variables['time'].ncattrs()})
		 			outVar = nc_out.createVariable('year', 'f', ('time',)) ; outVar[:]=year
		 			outVar = nc_out.createVariable('month', 'f', ('time',)) ; outVar[:]=month
		 			outVar = nc_out.createVariable(new_data.original_var_name, 'f', ('time','lat','lon',),fill_value=np.nan)
		 			for k in nc_in.variables[new_data.original_var_name].ncattrs():
		 				if k!='_FillValue':outVar.setncatts({k:nc_in.variables[new_data.original_var_name].getncattr(k)})
		 			outVar.setncattr('grid',new_data.grid)
					outVar[:] = new_data.raw

					nc_out.close()			
					nc_in.close()	

					os.system('rm '+self._working_directory+'/raw/*tmp*')

















