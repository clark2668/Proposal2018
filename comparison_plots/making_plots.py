# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp
import detector_data as detector


def main():



	ara2_logeV, ara2_limit = detector.get_limit('ara2_2016')
	ara_logeV, ara_limit = detector.get_limit('ara2_2016_singlestation_1year')
	arianna_logeV, arianna_limit = detector.get_limit('arianna_hra3_singlestation_1year')
	testbed_logeV, testbed_limit = detector.get_limit('ara_testbed_1year')

	ara_logeV, ara_aeff = detector.get_aeff('ara2_2016_single_fromlimit')
	arianna_logeV, arianna_aeff = detector.get_aeff('arianna_hra3_single_fromlimit')
	testbed_logeV, testbed_aeff = detector.get_aeff('ara_testbed_fromlimit')

	fig = plt.figure(figsize=(2.*11,8.5))
	ax_limit = fig.add_subplot(1,2,2)
	ax_aeff = fig.add_subplot(1,2,1)	


	ax_limit.plot(np.power(10.,ara_logeV),ara_limit,'-s', linewidth=2.0,color='blue',label=r'1 ARA 200m Station @ SP, 1 Year')
	ax_limit.plot(np.power(10.,arianna_logeV),arianna_limit/0.58,'-o', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB, 0.58 Year')
	ax_limit.plot(np.power(10.,testbed_logeV),testbed_limit,'-^', linewidth=2.0,color='red',label=r'1 ARA Testbed Surface Station @ SP, 1 Year')
	beautify_limit(ax_limit, 0)

	ax_aeff.plot(np.power(10.,ara_logeV),ara_aeff,'-s', linewidth=2.0,color='blue',label=r'1 ARA 200m Station @ SP')
	ax_aeff.plot(np.power(10.,arianna_logeV),arianna_aeff,'-o', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB')
	ax_aeff.plot(np.power(10.,testbed_logeV),testbed_aeff,'-^', linewidth=2.0,color='red',label=r'1 ARA Testbed Surface Station @ SP')
	beautify_aeff(ax_aeff)
	
	fig.savefig("current_limits_and_aeff.png",edgecolor='none',bbox_inches="tight") #save the figure


def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e3 #the lower x limit
	yup = 1.e10 #the uppper x limit
	this_ax.set_xlabel('EnergZy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('[A$\Omega]_{eff}$  [cm$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='lower right',title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	#this_ax.set_yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])



#pass it the axes and the number of theory curves you're including
#always pass the theory curves first
def beautify_limit(this_ax, num_theory):
	sizer=20
	xlow =  1.e15 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-20 #the lower y limit
	yup = 1e-10 #the upper y limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('E F(E) [$cm^{-2} s^{-1} sr^{-1}$]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_legend = this_ax.legend(loc='lower right',title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	# handles, labels = this_ax.get_legend_handles_labels()
	# legend1 = this_ax.legend(handles[0:num_theory], labels[0:num_theory], loc='lower left',title='Flux Models')
	# this_ax.add_artist(legend1)
	# legend2 = this_ax.legend(handles[num_theory:], labels[num_theory:], loc='upper right',title='Analysis Level')
	# setp(legend1.get_texts(), fontsize=17)
	# setp(legend1.get_title(), fontsize=17)
	# setp(legend2.get_texts(), fontsize=17)
	# setp(legend2.get_title(), fontsize=17)
	this_ax.grid()

#pass it the axes and the number of theory curves you're including
#always pass the theory curves first
def beautify_counts(this_ax):
	sizer=20
	xlow =  1.e15 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-20 #the lower y limit
	yup = 1e-10 #the upper y limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('Events',size=sizer)
	#this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	#this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='lower right',title='Analysis Level Event Counts')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)

main()
