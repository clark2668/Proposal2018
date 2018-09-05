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


def get_aeff(limit_name):
	"""
	get_limit

	Parameters
	----------
	limit_name: string
		the aeff curve you want
	Returns
	-------
	energy_bins: ndarray
		the energy bins the aeff is defined over in the units of log(eV)
	limit_values: ndarray
		the aeff values in the energy bins in units of cm^2 * str
	"""

	if(limit_name=='ara2_2016_single'):
		data = np.genfromtxt("data/ara2_aeff.csv",delimiter=',',skip_header=0,names=['energy','aeff'])
		#this digitized plot is in units of m^2, so multiply by 1e4 and 4pi
		aeff=data['aeff']*1.e4*4.*np.pi
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5])
		return energy, aeff

	if(limit_name=='arianna_sp'):
		data = np.genfromtxt("data/arianna_sp_aff.csv",delimiter=',',skip_header=1,names=['logenergy','aeff'])
		aeff=data['aeff']
		energy = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.])
		return energy, aeff