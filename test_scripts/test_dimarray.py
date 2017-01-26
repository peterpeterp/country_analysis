import dimarray as da
import numpy as np

lat=range(10)
lon=range(20)
time=range(5)
data=np.zeros([5,10,20])

models=['mod1','mod2','Na']
scenarios=['rcp26','rcp85','Na']
datasets=['ncep','cru','Na']

test1=da.DimArray(axes=[np.asarray(models),np.asarray(scenarios),np.asarray(datasets),time,lat,lon],dims=['model','scenario','dataset','time','lat','lon'])
test1['mod1']['rcp26']['Na']=data

test1['Na']['Na']['cru']=data









models2=['mod1','mod2','mod3']

test2=da.DimArray(axes=[np.asarray(models2),np.asarray(models2),np.asarray(models2),time,lat,lon],dims=['model','scenario','dataset','time','lat','lon'])
for model in models:
	test2[model]=test1[model]
test2['mod3']=data
