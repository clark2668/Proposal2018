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



	'''
	Load theory information
	'''
	ahlers_data = np.genfromtxt("data/ahlers_2012.csv",delimiter=',',skip_header=1,names=['energy','flux'])
	ahlers_energy_logev = ahlers_data['energy']
	ahlers_energy = np.power(10.,ahlers_energy_logev)
	ahlers_limit_log = ahlers_data['flux']
	ahlers_limit = np.power(10.,ahlers_limit_log)
	ahlers_interpolator = splrep(ahlers_energy_logev, ahlers_limit_log,k=4)

	icecube_energy_logev=np.array([15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5])
	icecube_energy = np.power(10.,icecube_energy_logev)

	#we can first have the IceCube thru-mu fluxes
	#we can have the nomina, and upper and lower 1 sigma possibilities
	def icecube_thrumu_function(E):
		return 3.03 * ((E/1.e14)**-2.19) * 1e-27
	def icecube_thrumu_upper_function(E):
		return 3.81 * ((E/1.e14)**-2.09) * 1e-27
	def icecube_thrumu_lower_function(E):
		return 2.34 * ((E/1.e14)**-2.29) * 1e-27
	

	def icecube_combined_function(E):
		return 6.7 * ((E/1.e14)**-2.50) * 1e-27
	def icecube_combined_upper_function(E):
		return 7.8 * ((E/1.e14)**-2.41) * 1e-27
	def icecube_combined_lower_function(E):
		return 5.5 * ((E/1.e14)**-2.59) * 1e-27
	#to compute to a curve on EF(E) we have to multiply by energy
	icecube_thrumu_efe = icecube_thrumu_function(icecube_energy) * icecube_energy
	icecube_thrumu_upper_efe = icecube_thrumu_upper_function(icecube_energy) * icecube_energy
	icecube_thrumu_lower_efe = icecube_thrumu_lower_function(icecube_energy) * icecube_energy

	icecube_combined_efe = icecube_combined_function(icecube_energy) * icecube_energy
	icecube_combined_upper_efe = icecube_combined_upper_function(icecube_energy) * icecube_energy
	icecube_combined_lower_efe = icecube_combined_lower_function(icecube_energy) * icecube_energy
	

	'''
	Load detector information
	'''
	ara_logeV, ara_100m_aeff = detector.get_aeff('ara_phased_100m')
	arianna_logeV, arianna_aeff = detector.get_aeff('arianna_icrc2017_fromlimit')

	ara_logeV, ara_100m_limit = detector.get_limit('ara_phased_100m_1year')
	arianna_logeV, arianna_limit = detector.get_limit('arianna_2017icrc_1year')

	#trim the ARIANNA sensitivity to match the ara array in energy binning
	arianna_logeV = arianna_logeV[:-2]
	arianna_aeff = arianna_aeff[:-2]
	arianna_limit = arianna_limit[:-2]
	#arianna_logeV = arianna_logeV[-1:]
	#arianna_aeff = arianna_aeff[-1:]
	#arianna_limit = arianna_limit[-1:]

	app_ara = ara_100m_aeff * const.SecPerYear
	app_arianna = arianna_aeff * const.SecPerYear * 0.58

	exisiting_logeV, existing_best_exp = detector.get_exposure('best_existing')
	target_logeV, target_exp = detector.get_exposure('target')

	icecube_logeV, icecube_limit = detector.get_limit('icecube_2018')
	anita_logeV, anita_limit = detector.get_limit('anita_2018')

	fig = plt.figure(figsize=(2.*11,8.5))
	ax_limit = fig.add_subplot(1,2,2)
	ax_aeff = fig.add_subplot(1,2,1)

	ax_limit.plot(ahlers_energy,ahlers_limit,'-', linewidth=3.0,color='gray',label="Ahlers 2012")
	ax_limit.fill_between(icecube_energy,icecube_thrumu_lower_efe,icecube_thrumu_upper_efe,facecolor='orange',alpha=0.3)
	ax_limit.fill_between(icecube_energy,icecube_combined_lower_efe,icecube_combined_upper_efe,facecolor='magenta',alpha=0.3)	
	ax_limit.plot(icecube_energy,icecube_thrumu_efe,'-.', linewidth=3.0,color='orange',label=r'IceCube Thru-Mu (E$^{-2.19}$)')
	ax_limit.plot(icecube_energy,icecube_combined_efe,':', linewidth=3.0,color='magenta',label='IceCube Combined (E$^{-2.50}$)')	
	ax_limit.plot(np.power(10.,ara_logeV),ara_100m_limit,'-^', linewidth=2.0,color='blue',label=r'1 ARA Phased 100m Station @ SP, 1 Year')
	ax_limit.plot(np.power(10.,arianna_logeV),arianna_limit/0.58,'-v', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB, 0.58 Year')
	ax_limit.plot(np.power(10.,icecube_logeV),icecube_limit/0.58,'->', linewidth=2.0,color='black',label=r'IceCube 2018')
	ax_limit.plot(np.power(10.,anita_logeV),anita_limit/0.58,'-<', linewidth=2.0,color='brown',label=r'ANITA 2018')	
	beautify_limit_withtheory(ax_limit,3)

	ax_aeff.plot(np.power(10.,exisiting_logeV),existing_best_exp,'--', linewidth=2.0,color='grey',label=r'Existing Best Exposure')
	ax_aeff.plot(np.power(10.,target_logeV),target_exp,'--', linewidth=2.0,color='red',label=r'NGRA Target Exposure')	
	ax_aeff.plot(np.power(10.,ara_logeV),app_ara,'-^', linewidth=2.0,color='blue',label=r'1 ARA Phased 100m Station @ SP, 1 Year')
	ax_aeff.plot(np.power(10.,arianna_logeV),app_arianna,'-v', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB, 0.58 Year')
	beautify_aeff(ax_aeff)
	fig.savefig("exposures_ara100phased_arianna.png",edgecolor='none',bbox_inches="tight") #save the figure



def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e11 #the lower x limit
	yup = 1.e21 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('Exposure=[TA$\Omega]_{eff}$  [s cm$^2$ sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='upper left')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)


#pass it the axes and the number of theory curves you're including
#always pass the theory curves first
def beautify_limit_withtheory(this_ax, num_theory):
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
	handles, labels = this_ax.get_legend_handles_labels()
	legend1 = this_ax.legend(handles[0:num_theory], labels[0:num_theory], loc='lower left',title='Flux Models')
	this_ax.add_artist(legend1)
	legend2 = this_ax.legend(handles[num_theory:], labels[num_theory:], loc='upper right',title='Analysis Level Single Event Sensitivity')
	setp(legend1.get_texts(), fontsize=17)
	setp(legend1.get_title(), fontsize=17)
	setp(legend2.get_texts(), fontsize=17)
	setp(legend2.get_title(), fontsize=17)
	this_ax.grid()

main()
