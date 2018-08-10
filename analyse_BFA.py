# -*- coding: utf-8 -*-

import sys,glob,os,pickle
import numpy as np
from netCDF4 import Dataset,num2date
import pandas as pd

iso = 'BFA'
print(iso)
os.chdir('/Users/peterpfleiderer/Projects/country_analysis/')
import country_analysis as country_analysis; reload(country_analysis)

# data will be stored in working_directory
COU=country_analysis.country_analysis(iso,working_directory='data/'+iso+'/',additional_tag='')


COU.load_data()

members_tas=COU.classify_ensemble(['tas'])
COU.ensemble_statistic(members_tas,'tas_ensmedian')
COU.save_data()

COU.period_statistic([1986,2006],'ref_mean')
COU.period_statistic([2025,2045],'2030s_mean')
COU.period_statistic([2035,2055],'2040s_mean')

COU.period_diff('ref_mean','2040s_mean','diff_ref_2040s')





#
