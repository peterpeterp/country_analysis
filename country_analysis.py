
#############
# 
#############

import sys,glob,os,itertools,datetime,pickle
import numpy as np
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import random as random
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pylab as plt 


class country_analysis(object):

	def __init__(self,iso,working_directory):
		self._iso=iso
		self._working_directory=working_directory

		if os.path.isdir(working_directory+iso)==False:os.system('mkdir '+working_directory+iso)
		if os.path.isdir(self._working_directory+self._iso+'/plots')==False:os.system('mkdir '+self._working_directory+self._iso+'/plots')


	def create_mask(self,filename,var_name,shape_file,shift_lon=0.0,mask_style='lat_weighted',pop_mask_file=''):

		if '_masks' not in dir(self): self._masks={}
		if os.path.isdir(self._working_directory+self._iso+'/masks')==False:os.system('mkdir '+self._working_directory+self._iso+'/masks')

		# get information about grid of input data
		nc_in=Dataset(filename,'r')
		input_data=nc_in.variables[var_name][1,:,:]
		lat = nc_in.variables['lat'][:]								
		lon = nc_in.variables['lon'][:].squeeze()+shift_lon 	
		nx = len(lon)	;	ny = len(lat)
		grid=str(ny)+'x'+str(nx)

		if grid not in self._masks.keys():self._masks[grid]={}
		if mask_style not in self._masks[grid].keys():self._masks[grid][mask_style]={}


		mask_file=self._working_directory+self._iso+'/masks/'+self._iso+'_'+grid+'_'+mask_style+'.nc4'

		if os.path.isfile(mask_file):
			# load existing mask
			nc_mask=Dataset(mask_file,'r')
			self._masks[grid][mask_style] = nc_mask.variables[self._iso][:,:]  
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
			lon=lon-shift_lon
			shift = len(lon)-np.where(lon==lon[0]-shift_lon)[0][0]

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
				mygrid=open(self._working_directory+self._iso+'/masks/'+grid+'.txt','w')
				mygrid.write('gridtype=lonlat\nxsize='+str(len(lon))+'\nysize='+str(len(lat))+'\nxfirst='+str(lon[0])+'\nxinc='+str(np.mean(np.diff(lon,1)))+'\nyfirst='+str(lat[0])+'\nyinc='+str(np.mean(np.diff(lat,1))))
				mygrid.close()
				os.system('cdo remapbil,'+self._working_directory+self._iso+'/masks/'+grid+'.txt '+pop_mask_file+' '+self._working_directory+self._iso+'/masks/'+grid+'_'+mask_style+'.nc')	
				nc_pop_mask = Dataset(self._working_directory+self._iso+'/masks/'+grid+'_'+mask_style+'.nc')
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
				self._masks[grid][mask_style]=np.roll(output,shift,axis=1)

			# save mask
			print mask_file
			nc_mask=Dataset(mask_file,'w')
			nc_mask.createDimension('lat', len(lat))
			nc_mask.createDimension('lon', len(lon))
 			outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
 			outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')
 			outVar = nc_mask.createVariable(self._iso, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[grid][mask_style][:,:]

 			# close nc_files
			nc_in.close()
			nc_mask.close()

	def country_zoom(self,in_file,var_name,meta_data=['CMIP5','rcp2p6','hadgem2-es'],mask_style='lat_weighted',time_units=None,time_calendar=None):
		'''
		compute weighted country average for each timestep
		in_file: type str: file to be processed
		out_file: type str: filepath where output is going to be stored
		var: type str: variable name in netcdf file
		iso: type str: iso3c of country
		mask_path: type str: path to where the masks are stored
		'''

		if '_data' not in dir(self): 
			self._data={}
			self._meta=[]
		if os.path.isdir(self._working_directory+self._iso+'/raw')==False:os.system('mkdir '+self._working_directory+self._iso+'/raw')


		out_file=self._working_directory+self._iso+'/raw/'+in_file.split('/')[-1].replace('.nc','_'+self._iso+'.nc')
		print out_file

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

			country_mask=self._masks[grid][mask_style]
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


	def average(self,mask_style):
		'''
		mask_style: str: type of mask weighting 
		computes country average for all loaded datasets
		'''

		# preprare dictionary structure as for self._data
		if '_transcient' not in dir(self): self._transcient={}
		if os.path.isdir(self._working_directory+self._iso+'/country_average')==False:os.system('mkdir '+self._working_directory+self._iso+'/country_average')

		for meta_list in self._meta:
			tmp_in = self._data
			tmp_out = self._transcient
			for meta_info in meta_list:
				if meta_info not in tmp_out: tmp_out[meta_info]={}
				tmp_in=tmp_in[meta_info]
				tmp_out=tmp_out[meta_info]

			# check if file has been saved
			out_file=self._working_directory+self._iso+'/country_average/country_mean_'+'_'.join(meta_list)+'_'+mask_style+'.csv'
			if os.path.isfile(out_file):
				tmp_out['time']=pd.read_csv(out_file,sep=';')['time']
				tmp_out['month']=pd.read_csv(out_file,sep=';')['month']
				tmp_out['year']=pd.read_csv(out_file,sep=';')['year']
				tmp_out[mask_style]=pd.read_csv(out_file,sep=';')['_'.join(meta_list)]		

			# if not compute
			else:
				# load input data
				var_in=tmp_in['data']			

				# get mask
				mask=self._masks[tmp_in['grid']][mask_style]
				lat_mask=self._masks[tmp_in['grid']]['lat_mask']
				lon_mask=self._masks[tmp_in['grid']]['lon_mask']

				# find relevant area (as rectangle) and check whether lon and lat are correct (on same grid differences in lat decreasing or increasing could arise)
				lat_mean=np.mean(mask,1)
				lats=np.where(lat_mean!=0)
				if lat_mask[lats][0]!=tmp_in['lat'][0]:
					print lat_mask[lats][0],tmp_in['lat'][0]
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


	def plot_transcient(self,meta_data,mask_style):
		'''

		'''

		tmp = self._transcient
		for meta_info in meta_data:
			tmp=tmp[meta_info]

		fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(5,5))

		#plt.plot(tmp['time'],tmp[mask_style])
		plt.plot(tmp['time'],tmp[mask_style],linestyle=':',color='blue')
		plt.plot(tmp['time'],pd.rolling_mean(tmp[mask_style],20),linestyle='-',color='red')

		plt.savefig(self._working_directory+self._iso+'/plots/'+'_'.join([self._iso]+meta_data+[mask_style])+'.png')



	def period_averages(self,periods={'ref':[1986,2006],'2030s':[2025,2045],'2040s':[2035,2055]}):
		'''

		'''

		# preprare dictionary structure as for self._data
		for meta_list in self._meta:
			tmp = self._data
			for meta_info in meta_list:
				tmp=tmp[meta_info]

			tmp['period']={}
			tmp['period_meta']=periods

			for period_name in periods.keys():
				tmp['period'][period_name]=np.mean(tmp['data'][np.where((tmp['year']>=periods[period_name][0]) & (tmp['year']<periods[period_name][1]))[0],:,:],axis=0)

			for period_name in periods.keys():
				if period_name!='ref' and 'ref' in periods.keys():
					tmp['period'][period_name+'-'+'ref']=tmp['period'][period_name]-tmp['period']['ref']



















