import numpy as np
from mpl_toolkits.basemap import Basemap, cm
import mpl_toolkits.basemap
import matplotlib.pylab as plt
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator

def plot_map(ax,lon,lat,Z,color_palette=plt.cm.plasma,color_range=None,color_label=None,grey_area=None,limits=None):
	# handle 0 to 360 lon
	if max(lon)>180:
		problem_start=np.where(lon>180)[0][0]
		new_order=np.array(range(problem_start,len(lon))+range(0,problem_start))
		Z=Z[:,new_order]
		lon=lon[new_order]
		lon[lon>180]-=360

	# handle limits
	if limits==None:
		limits=[np.min(lon),np.max(lon),np.min(lat),np.max(lat)]
	if limits!=None:
		lon_select=np.where((lon>=limits[0])	&	(lon<=limits[1]))[0]
		lat_select=np.where((lat>=limits[2])	&	(lat<=limits[3]))[0]
		Z=Z[lat_select,:][:,lon_select]

	m = Basemap(ax=ax,llcrnrlon=limits[0],urcrnrlon=limits[1],llcrnrlat=limits[2],urcrnrlat=limits[3],resolution="l",projection='cyl')
	m.drawmapboundary(fill_color='1.')

	# imshow does not support decreasing lat or lon
	Z=np.ma.masked_invalid(Z)
	if lat[0]>lat[1]:Z=Z[::-1,:]
	if lon[0]>lon[1]:Z=Z[:,::-1]



	# get volor_range
	if color_range==None:
		color_range=[np.min(Z[np.isfinite(Z)]),np.max(Z[np.isfinite(Z)])]

	im = m.imshow(Z,cmap=color_palette,vmin=color_range[0],vmax=color_range[1],interpolation='none',extent=[np.min(lon),np.max(lon),np.min(lat),np.max(lat)])

	# mask some grid-cells
	if grey_area!=None:
		Z=np.ma.masked_invalid(grey_area.copy())
		if lat[0]>lat[1]:Z=Z[::-1,:]
		if lon[0]>lon[1]:Z=Z[:,::-1]
		im2 = m.imshow(Z,cmap=plt.cm.Greys,vmin=0,vmax=1,interpolation='none',extent=[np.min(lon),np.max(lon),np.min(lat),np.max(lat)])

	# show coastlines and borders
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()

	# add colorbar
	if color_label!=None:
		cb = m.colorbar(im,'right', size="5%", pad="2%")
		tick_locator = ticker.MaxNLocator(nbins=5)
		cb.locator = tick_locator
		cb.update_ticks()
		cb.set_label(color_label, rotation=90)

	return(ax,im)






