# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import interp1d, splrep, splev
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

	if(resource_name=='ara2_2016_trigger'):
		#this is the ARA2 station limit (7.5 months x Two Stations) at the analysis level
		#direct digitization of Fig 37 in https://arxiv.org/abs/1507.08991
		data = np.genfromtxt("data/ara2_limit_trigger.csv",delimiter=',',skip_header=1,names=['energy','limit'])
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

	if(resource_name=='ara_200m_1year_trigger'):
		#scaling of the ara 2 station limit at trigger level
		#to get to 1 year, 1 station value
		#we want to rip out the two station fact (2.)
		#and rescale to get 1 year: 225./365.
		energy, limit = get_limit('ara2_2016_trigger')
		limit= limit*2.*(225./365.)

		return energy, limit

	if(resource_name=='ara_100m_1year'):
		#going to invert the 100m Aeff which we got by geometric averaging the 200m and Testbed
		#this one is a touch hokey because we don't actually have this measurement
		energy, aeff = get_aeff('ara_100m')
		limit= 2.44/np.log(10.)/aeff/(const.SecPerYear)

		return energy, limit

	if(resource_name=='ara_phased_1year'):
		#take the ARA 200m station, slide it over by a factor of two
		#then re-interpolate it to the original half-energy decade bins
		#the interpolation is most robust in log-log space where the curve is well behaved
		energy, limit = get_limit('ara_200m_1year')
		energy = np.log10(np.power(10.,energy)/2.) #move the limit over by a factor of two in linear space
		interpolator = splrep(energy, np.log10(limit),k=2)
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20])
		limit_interp = np.power(10.,splev(energy, interpolator))
		return energy, limit_interp

	if(resource_name=='ara_phased_100m_1year'):
		#take the ARA 200m station, rescale it to 100m, and slide it over by a factor of two
		#then re-interpolate it to the original half-energy decade bins
		#the interpolation is most robust in log-log space where the curve is well behaved
		energy, limit = get_limit('ara_100m_1year')
		energy = np.log10(np.power(10.,energy)/2.) #move the limit over by a factor of two in linear space
		interpolator = splrep(energy, np.log10(limit),k=2)
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20])
		limit_interp = np.power(10.,splev(energy, interpolator))
		return energy, limit_interp

	if(resource_name=='ara_phased_100m_1year_trigger'):
		#take the ARA 200m station, rescale it to 100m, and slide it over by a factor of two
		#then re-interpolate it to the original half-energy decade bins
		#the interpolation is most robust in log-log space where the curve is well behaved
		energy, limit = get_limit('ara_100m_1year_trigger')
		energy = np.log10(np.power(10.,energy)/2.) #move the limit over by a factor of two in linear space
		interpolator = splrep(energy, np.log10(limit),k=2)
		energy = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20])
		limit_interp = np.power(10.,splev(energy, interpolator))
		return energy, limit_interp

	if(resource_name=='arianna_hra3_singlestation_1year'):
		#this is the ARIANNA HRA3 limit (170 days x three Stations) at the analysis level
		#scaling of Fig 22 in https://arxiv.org/abs/1410.7352
		#to get to 0.58 year 1 station value
		# we want to rip out the three station factor (3.)
		#and rescale to get 1 year: 170./365./0.58
		data = np.genfromtxt("data/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']/(data['energy']/1e9)*3.*(170./365.)
		energy = np.array([17.,17.5,18.,18.5,19.,19.5,20.,20.5])

		return energy, limit

	if(resource_name=='icecube_2018'):
		#this is the IceCube Limit
		#Fig 6 of https://arxiv.org/abs/1807.01820
		#we need to kill the E^2 bit
		data = np.genfromtxt("data/icecube_2018.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']/(data['energy'])
		#energy = np.log10(data['energy']*1e9)
		energy = np.array([16,16.5,17.,17.5,18.,18.5,19.,19.5,20.])
		#energy = data['energy']
		#limit = np.power(10.,data['limit'])

		return energy, limit

	if(resource_name=='anita_2018'):
		#this is the ANITA 3-flight limit
		#Fig 6 of https://arxiv.org/abs/1803.02719
		data = np.genfromtxt("data/anita_threeflight_2018.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']
		energy = np.array([18.,18.5,19.,19.5,20.,20.5,21])

		return energy, limit

	if(resource_name=='arianna_2017icrc_1year'):
		#this is the 1296, 5 year ARIANNA ICRC
		#Fig 6 of https://pos.sissa.it/301/977/pdf
		data = np.genfromtxt("data/arianna_2017icrc.csv",delimiter=',',skip_header=1,names=['energy','limit'])
		limit=data['limit']/(data['energy'])*1296.*(1.5*143./365.)*5.
		energy = np.array([15.5, 16, 16.5, 17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.])

		return energy, limit

	if(resource_name=='arianna_2017icrc_1year_trigger'):
		#undo 90% analysis efficiency
		energy, limit = get_limit(arianna_2017icrc_1year)
		limit = limit * 0.9 #remove the analysis efficiency of 90% (make it *smaller*)

		return energy, limit

def get_exposure(resource_name):
	"""
	get_limit

	Parameters
	----------
	limit_name: string
		the exposure curve you want
	Returns
	-------
	energy_bins: ndarray
		the energy bins the aeff is defined over in the units of log10(eV)
	expsoure_values: ndarray
		the exposure values in the energy bins in units of cm^2 * sr * s
	"""
	if(resource_name=='best_existing'):

		icecube_logeV, icecube_limit = get_limit('icecube_2018')
		anita_logeV, anita_limit = get_limit('anita_2018')

		icecube_exposure = 2.44/np.log(10)/0.5/icecube_limit
		anita_exposure = 2.44/np.log(10)/0.5/anita_limit
	
		existing_exposure = icecube_exposure[:8]
		existing_energy = icecube_logeV[:8]

		#this is the exposure we need to get our hands on
		existing_exposure = np.append(existing_exposure,anita_exposure[4:5])
		existing_energy = np.append(existing_energy,20)

		return existing_energy, existing_exposure

	if(resource_name=='target'):

		energy, exposure = get_exposure('best_existing')
		target_exposure = 10. * exposure

		return energy, target_exposure

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
		the energy bins the aeff is defined over in the units of log10(eV)
	limit_values: ndarray
		the aeff values in the energy bins in units of cm^2 * sr
	"""


	##first, all the aeffs which come from inverting limit curves


	if(resource_name=='arianna_icrc2017_fromlimit'):
		#this comes from inverting the analysis level limit curve of Fig 6 in https://pos.sissa.it/301/977/pdf
		#remove the statistical factor of 2.3 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('arianna_2017icrc_1year')
		single_station_aeff = 2.3/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	if(resource_name=='arianna_icrc2017_fromlimit_trigger'):
		#this comes from inverting the analysis level limit curve of Fig 6 in https://pos.sissa.it/301/977/pdf
		#remove the statistical factor of 2.3 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year
		#make this a trigger level curve

		energy_logev, limit = get_limit('arianna_2017icrc_1year_trigger')
		single_station_aeff = 2.3/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

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

	if(resource_name=='ara_200m_1year_fromlimit_trigger'):
		#this comes from inverting the trigger level limit curve of Fig 37 in https://arxiv.org/abs/1507.08991	
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_200m_1year_trigger')
		single_station_aeff = 2.44/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	if(resource_name=='ara_100m'):
		ara_logeV, ara_200m_aeff = get_aeff('ara_200m_1year_fromlimit')
		testbed_logeV, testbed_aeff =get_aeff('ara_testbed_fromlimit')
		testbed_logeV=testbed_logeV[:-1]
		testbed_aeff=testbed_aeff[:-1]
	
		ara_logeV=ara_logeV[2:]	
		ara_200m_aeff=ara_200m_aeff[2:]
	
		final = np.sqrt(testbed_aeff*ara_200m_aeff)

		interpolator = splrep(testbed_logeV, np.log10(final),k=2)
		test_energy = np.arange(16,21,.5)
		aeff_interp = np.power(10.,splev(test_energy, interpolator))
		aeff_interp[0]=aeff_interp[0]*.2
		aeff_interp[1]=aeff_interp[1]*.8
		return test_energy, aeff_interp

	if(resource_name=='ara_100m_trigger'):
		ara_logeV, ara_200m_aeff = get_aeff('ara_200m_1year_fromlimit_trigger')
		testbed_logeV, testbed_aeff =get_aeff('ara_testbed_fromlimit')
		testbed_logeV=testbed_logeV[:-1]
		testbed_aeff=testbed_aeff[:-1]
	
		ara_logeV=ara_logeV[2:]	
		ara_200m_aeff=ara_200m_aeff[2:]
	
		final = np.sqrt(testbed_aeff*ara_200m_aeff)

		interpolator = splrep(testbed_logeV, np.log10(final),k=2)
		test_energy = np.arange(16,21,.5)
		aeff_interp = np.power(10.,splev(test_energy, interpolator))
		aeff_interp[0]=aeff_interp[0]*.2
		aeff_interp[1]=aeff_interp[1]*.8
		return test_energy, aeff_interp


	if(resource_name=='ara_phased_100m'):
		#this comes from inverting the prediction for the 100m ARA phased enhaned array station
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_phased_100m_1year')
		single_station_aeff = 2.44/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	if(resource_name=='ara_phased_100m_trigger'):
		#this comes from inverting the prediction for the 100m ARA phased enhaned array station
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_phased_100m_1year_trigger')
		single_station_aeff = 2.44/np.log(10)/0.5/limit/(const.SecPerYear)
		
		return energy_logev, single_station_aeff

	# if(resource_name=='ara_100m'):
	# 	#we are going to rescale the 200m value to land roughly in between the ARA deep and ARA testbed

	# 	energy_logev, aeff = get_aeff('ara_200m_1year_fromlimit')
	# 	single_station_aeff = aeff/5.
		
	# 	return energy_logev, single_station_aeff

	if(resource_name=='ara_phased_fromlimit'):
		#this comes from inverting the analysis level limit curve of the phased array estimate
		#remove the statistical factor of 2.44 for 90% CL
		#remove ln10 for conversion to logarithmic bins
		#remove half-decade wide energy bins (0.5)
		#remove seconds per year

		energy_logev, limit = get_limit('ara_phased_1year')
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