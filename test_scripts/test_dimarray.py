import dimarray as da
import numpy as np

lat=range(10)
lon=range(20)
time=range(5)
models=['mod1','mod2']
data=np.zeros([5,10,20])

test1=da.DimArray(axes=[np.asarray(models),time,lat,lon],dims=['model','time','lat','lon'])
test1['mod1']=data

models2=['mod1','mod2','mod3']

test2=da.DimArray(axes=[np.asarray(models2),time,lat,lon],dims=['model','time','lat','lon'])
for model in models:
	test2[model]=test1[model]
test2['mod3']=data
