# -*- coding: utf-8 -*-
import numpy as np #import numpy
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp
import constants as const

def get_counts(resource_name):
	"""
	get_counts

	Parameters
	----------
	logev: ndarray
		energies in log10(eV)
	curve: ndarray
		the values we want to integrate

	Returns
	-------
	counts: ndarray
		Counts per energy bin
	"""

	'''
	To be done
	'''