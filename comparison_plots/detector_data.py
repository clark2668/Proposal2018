# -*- coding: utf-8 -*-
import numpy as np #import numpy
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp
import constants as const

def get_limit(resource_name):
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
	if(resource_name=='ara_testbed'):
		#this is the ARA testbed limit (445 days x 1 Station) at the analysis level
		#direct digitization of Fig 13 in https://arxiv.org/abs/1404.5285
		data = np.genfromtxt("data/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']
		energy = np.array([17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.])

		return energy, limit

	if(resource_name=='ara_testbed_1year'):
		#scaling of the ara testbed limit
		#to get to 1 year, 1 station value
		#and rescale to get 1 year: 445./365.
		energy, limit = get_limit('ara_testbed')
		limit= limit*(445./365.)

		return energy, limit

	if(resource_name=='ara2_2016'):
		#this is the ARA2 station limit (7.5 months x Two Stations) at the analysis level
		#direct digitization of Fig 37 in https://arxiv.org/abs/1507.08991
		data = np.genfromtxt("data/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5])

		return energy, limit

	if(resource_name=='ara_200m_1year'):
		#scaling of the ara 2 station limit
		#to get to 1 year, 1 station value
		#we want to rip out the two station fact (2.)
		#and rescale to get 1 year: 225./365.
		energy, limit = get_limit('ara2_2016')
		limit= limit*2.*(225./365.)

		return energy, limit

	if(resource_name=='arianna_hra3_singlestation_1year'):
		#this is the ARIANNA HRA3 limit (170 days x three Stations) at the analysis level
		#scaling of Fig 22 in https://arxiv.org/abs/1410.7352
		data = np.genfromtxt("data/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']/(data['energy']/1e9)*3.*(170./365.)
		#to get to 0.58 year 1 station value
		# we want to rip out the three station factor (3.)
		#and rescale to get 1 year: 170./365./0.58
		energy = np.array([17.,17.5,18.,18.5,19.,19.5,20.,20.5])

		return energy, limit


def get_aeff(resource_name):
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


	##first, all the aeffs which come from inverting limit curves

	if(resource_name=='ara_testbed_fromlimit'):
		#this comes from inverting the analysis level limit curve of Fig 37 in https://arxiv.org/abs/1507.08991
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_testbed_1year')
		single_station_aeff = 2.3/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	if(resource_name=='ara_200m_1year_fromlimit'):
		#this comes from inverting the analysis level limit curve of Fig 37 in https://arxiv.org/abs/1507.08991	
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_200m_1year')
		single_station_aeff = 2.44/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	if(resource_name=='arianna_hra3_single_fromlimit'):
		#this comes from inverting the analysis level limit curve of Fig 22 in https://arxiv.org/abs/1410.7352
		#remove the statistical factor of 2.3 for 90% CL (ARIANNA uses 2.3)
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year
		
		energy_logev, limit = get_limit('arianna_hra3_singlestation_1year')
		single_station_aeff = 2.3/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	##then, all the aeffs which come from direct digitiziations

	if(resource_name=='arianna_sp'):
		#this come straight from the google drive sensitivty curve by Chris Persichelli
		data = np.genfromtxt("data/arianna_sp_aff.csv",delimiter=',',skip_header=1,names=['logenergy','aeff'])
		aeff=data['aeff']
		energy = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.])
		return energy, aeff

	if(resource_name=='ara_200m_1year_fromfigure'):
		#this come straight from a digitization of Fig 14 in https://arxiv.org/abs/1507.08991
		data = np.genfromtxt("data/ara2_aeff.csv",delimiter=',',skip_header=0,names=['energy','aeff'])

		aeff=data['aeff']*1.e4*4.*np.pi
		#this digitized plot is in units of m^2, so multiply by 1e4 and 4pi
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5])

		return energy, aeff