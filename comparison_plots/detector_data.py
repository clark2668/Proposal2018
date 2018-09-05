# -*- coding: utf-8 -*-
import numpy as np #import numpy
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp

def get_limit(limit_name):
	"""
	get_limit

	Parameters
	----------
	limit_name: string
		the limit curve you want to plot
	Returns
	-------
	energy_bins: ndarray
		the energy bins the limit is defined over in the units of log10(eV)
	limit_values: ndarray
		the limit in the energy bins in units of 1/cm^2/s/sr
	"""

	if(limit_name=='ara2_2016'):
		data = np.genfromtxt("data/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5])
		return energy, limit