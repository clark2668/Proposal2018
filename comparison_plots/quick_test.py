# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp
import detector_data as detector
import constants as const


def main():

	ara_logeV, ara_200m_aeff = detector.get_aeff('ara_200m_1year_fromlimit')
	phased_logeV, phased_aeff = detector.get_aeff('ara_phased_fromlimit')

	fig = plt.figure(figsize=(2.*11,8.5))
	ax_limit = fig.add_subplot(1,2,2)

	ax_aeff.plot(np.power(10.,ara_logeV),ara_200m_aeff,'-<', linewidth=2.0,color='blue',label=r'1 ARA 200m Station @ SP')
	beautify_aeff(ax_aeff)
	fig.savefig("try_super_station.png",edgecolor='none',bbox_inches="tight") #save the figure


def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e2 #the lower x limit
	yup = 1.e10 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('[A$\Omega]_{eff}$  [cm$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='lower left',title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	#this_ax.set_yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])

main()
