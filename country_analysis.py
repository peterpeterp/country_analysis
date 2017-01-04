
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

	def __init__(self,iso,var_name,working_directory):
		self._iso=iso
		self._var_name=var_name
		self._working_directory=working_directory


	def create_mask(self,filename,var_name,shape_file,shift_lon=0.0,mask_style='lat_weighted'):

		if '_masks' not in dir(self): self._masks={}

		mask_file=self._working_directory+'tmp/masks/'+self._iso+'_'+mask_style+'.nc4'

		if os.path.isfile(mask_file):
			# load existing mask
			nc_mask=Dataset(mask_file,'r')
			self._masks[mask_style] = nc_mask.variables[self._iso][:,:]  
			self._masks['lat_mask'] = nc_mask.variables['lat'][:]  
			self._masks['lon_mask'] = nc_mask.variables['lon'][:]  

		else:
			# get information about grid of input data
			nc_in=Dataset(filename,'r')
			input_data=nc_in.variables[var_name][1,:,:]
			lat = nc_in.variables['lat'][:]								; self._masks['lat_mask']=lat
			lon = nc_in.variables['lon'][:].squeeze()+shift_lon 		; self._masks['lon_mask']=lon
			nx = len(lon)
			ny = len(lat)

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

			overlap = np.zeros((ny,nx))
			for i in range(nx):
				for j in range(ny):
					# check gridcell is relevant
					if grid_polygons[i,j].intersects(ext_poly):
						# get fraction of grid-cell covered by polygon
						intersect = grid_polygons[i,j].intersection(country_polygons).area/grid_polygons[i,j].area*country_polygons.area
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
				self._masks[mask_style]=np.roll(output,shift,axis=1)

			# save mask
			nc_mask=Dataset(mask_file,'w')
			nc_mask.createDimension('lat', len(lat))
			nc_mask.createDimension('lon', len(lon))
 			outVar = nc_mask.createVariable('lat', 'f', ('lat',)) ; outVar[:]=lat[:]	;	outVar.setncattr('units','deg south')
 			outVar = nc_mask.createVariable('lon', 'f', ('lon',)) ; outVar[:]=lon[:]	;	outVar.setncattr('units','deg east')
 			outVar = nc_mask.createVariable(self._iso, 'f', ('lat','lon',),fill_value='NaN') ; outVar[:]=self._masks[mask_style][:,:]

 			# close nc_files
			nc_in.close()
			nc_mask.close()

	def country_zoom(self,in_file,var_name,mask_style='lat_weighted',meta_data=['CMIP5','rcp2p6','hadgem2-es']):
		'''
		compute weighted country average for each timestep
		in_file: type str: file to be processed
		out_file: type str: filepath where output is going to be stored
		var: type str: variable name in netcdf file
		iso: type str: iso3c of country
		mask_path: type str: path to where the masks are stored
		'''

		if '_data' not in dir(self): self._data={}

		# open file to get information
		print in_file
		nc_in=Dataset(in_file,"r")
		lon_in=nc_in.variables['lon'][:]
		lat_in=nc_in.variables['lat'][:]

		# find relevant area (as rectangle)
		lon_mean=np.mean(self._masks[mask_style],0)
		lons=np.where(lon_mean!=0)[0]

		lat_mean=np.mean(self._masks[mask_style],1)
		lats=np.where(lat_mean!=0)[0]

		nx,ny=len(lons),len(lats)

		# copy netcdf and write zoomed file
		out_file=self._working_directory+'tmp/raw/'+in_file.split('/')[-1].replace('.nc','_'+self._iso+'.nc')
		print out_file
		os.system("rm "+out_file)
		nc_out=Dataset(out_file,"w")
		for dname, the_dim in nc_in.dimensions.iteritems():
			if dname=='lon':nc_out.createDimension(dname, nx)
			elif dname=='lat':nc_out.createDimension(dname, ny)
			else:nc_out.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

		# Copy variables
		for v_name, varin in nc_in.variables.iteritems():
			print v_name
			outVar = nc_out.createVariable(v_name, varin.datatype, varin.dimensions)
						    
			# Copy variable attributes
			outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
						    
			if v_name=='lat':	outVar[:] = self._masks['lat_mask'][list(lats)]
			elif v_name=='lon':	outVar[:] = self._masks['lon_mask'][list(lons)]

			elif v_name==var_name:	
				# check whether lon and lat are similarly defined in mask and in_file
				if (lat_in in self._masks['lat_mask']) == False:	lats= -np.array(lats)+len(lat_in)
				if (lon_in in self._masks['lon_mask']) == False:	lons= -np.array(lons)+len(lon_in)
				lats,lons=sorted(lats),sorted(lons)

				var_in=nc_in.variables[var_name][:,list(lats),list(lons)]
				try:	# handle masked array
					masked=np.ma.getmask(var_in)
					var_in=np.ma.getdata(var_in)
					var_in[masked]=np.nan
				except: pass
				# creat a 1-NA mask
				red_mask = self._masks[mask_style][np.ix_(list(lats),list(lons))]
				red_mask[red_mask>0]=1
				red_mask[red_mask==0]=np.nan
				country_data=var_in*red_mask
				outVar[:] = country_data[:,:,:]

			else:	outVar[:] = varin[:]

		# close the output file
		nc_out.close()
		nc_in.close()

		tmp = self._data
		for meta_info in meta_data:
			if meta_info not in tmp: tmp[meta_info]={}
			tmp=tmp[meta_info]

		tmp[self._var_name]=country_data
		tmp['lat']=self._masks['lat_mask'][list(lats)]
		tmp['lon']=self._masks['lon_mask'][list(lons)]



