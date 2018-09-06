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

	ara_logeV, ara_200m_limit = detector.get_limit('ara_200m_1year')
	arianna_logeV, arianna_limit = detector.get_limit('arianna_hra3_singlestation_1year')
	testbed_logeV, testbed_limit = detector.get_limit('ara_testbed_1year')
	ara_100m_logeV, ara_100m_limit = detector.get_limit('ara_100m_1year')
	phased_logeV, phased_limit = detector.get_limit('ara_phased_1year')

	ara_logeV, ara_200m_aeff = detector.get_aeff('ara_200m_1year_fromlimit')
	arianna_logeV, arianna_aeff = detector.get_aeff('arianna_hra3_single_fromlimit')
	testbed_logeV, testbed_aeff = detector.get_aeff('ara_testbed_fromlimit')
	ara_100m_logeV, ara_100m_aeff = detector.get_aeff('ara_100m_1year_fromlimit')
	phased_logeV, phased_aeff = detector.get_aeff('ara_phased_fromlimit')

	fig = plt.figure(figsize=(2.*11,8.5))
	ax_limit = fig.add_subplot(1,2,2)
	ax_aeff = fig.add_subplot(1,2,1)

	ax_limit.plot(np.power(10.,phased_logeV),phased_limit,'-v', linewidth=2.0,color='magenta',label=r'1 ARA Phased 200m Station @ SP, 1 Year')
	ax_limit.plot(np.power(10.,ara_logeV),ara_200m_limit,'-<', linewidth=2.0,color='blue',label=r'1 ARA 200m Station @ SP, 1 Year')
	ax_limit.plot(np.power(10.,ara_100m_logeV),ara_100m_limit,'->', linewidth=2.0,color='black',label=r'1 ARA 100m Station @ SP, 1 Year (approx)')
	ax_limit.plot(np.power(10.,arianna_logeV),arianna_limit/0.58,'-o', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB, 0.58 Year')
	ax_limit.plot(np.power(10.,testbed_logeV),testbed_limit,'-^', linewidth=2.0,color='red',label=r'1 ARA Testbed Surface Station @ SP, 1 Year')
	beautify_limit(ax_limit)

	ax_aeff.plot(np.power(10.,phased_logeV),phased_aeff,'-<', linewidth=2.0,color='magenta',label=r'1 ARA Phased 200m Station @ SP, 1 Year')
	ax_aeff.plot(np.power(10.,ara_logeV),ara_200m_aeff,'-<', linewidth=2.0,color='blue',label=r'1 ARA 200m Station @ SP')
	ax_aeff.plot(np.power(10.,ara_100m_logeV),ara_100m_aeff,'->', linewidth=2.0,color='black',label=r'1 ARA 100m Station @ SP (approx)')
	ax_aeff.plot(np.power(10.,arianna_logeV),arianna_aeff,'-o', linewidth=2.0,color='green',label=r'1 ARIANNA Surface Station @ MB')
	ax_aeff.plot(np.power(10.,testbed_logeV),testbed_aeff,'-^', linewidth=2.0,color='red',label=r'1 ARA Testbed Surface Station @ SP')
	beautify_aeff(ax_aeff)
	fig.savefig("current_single_station_limits_and_aeff.png",edgecolor='none',bbox_inches="tight") #save the figure



	#okay, now we can glue together all our copies of the various detectors and scale them for livetime

	#need to set up joing limit energy bins
	total_energy_logeV = np.array([16.,16.5,17.,17.5,18.,18.5,19,19.5,20.,20.5,21.])

	#need to insert two zeros in both ARA testbed and ARIANNA station
	#for 16 and 16.5 energy bins (no data there)
	testbed_aeff = np.insert(ara_200m_aeff,0,0.,axis=0)
	testbed_aeff = np.insert(ara_200m_aeff,0,0.,axis=0)
	arianna_limit = np.insert(arianna_limit,0,0.,axis=0)
	arianna_limit = np.insert(arianna_limit,0,0.,axis=0)

	#need to insert two zeros for 20.5 and 21 in phased array station
	phased_aeff = np.append(phased_aeff,0.)
	phased_aeff = np.append(phased_aeff,0.)

	#need to insert a 21 bin in ARA200m, ARA100m, ARIANNA
	ara_200m_aeff = np.append(ara_200m_aeff,0.)
	ara_100m_aeff = np.append(ara_100m_aeff,0.)
	arianna_aeff = np.append(arianna_aeff,0.)

	#number of stations of each type of station
	num_testbed = 1.
	num_100m = 1.
	num_200m = 3.
	num_phased = 1.

	#livetime in seconds for each type of station
	live_testbed = 632. * const.SecPerDay
	live_100m = 731. * const.SecPerDay
	live_200m = 2991. * const.SecPerDay
	live_phased = 186. * const.SecPerDay

	#we multiply the number of stations by the livetime
	testbed_app_total = num_testbed * live_testbed * testbed_aeff
	ara_100m_app_total = num_100m * live_100m * ara_100m_aeff
	ara_200m_app_total = num_200m * live_200m * ara_200m_aeff
	phased_app_total = num_phased * live_phased * phased_aeff

	#add!
	#units here are seconds * cm^2 * sr
	total_aeff = (num_testbed*testbed_aeff) + (num_100m*ara_100m_aeff) + (num_200m*ara_200m_aeff) + (num_phased * phased_aeff)
	total_app = testbed_app_total + ara_100m_app_total + ara_200m_app_total + phased_app_total
	total_limit = 2.44/np.log(10)/0.5/total_app

	ahlers_data = np.genfromtxt("data/ahlers_2012.csv",delimiter=',',skip_header=1,names=['energy','flux'])
	ahlers_energy_logev = ahlers_data['energy']
	ahlers_energy = np.power(10.,ahlers_energy_logev)
	ahlers_limit_log = ahlers_data['flux']
	ahlers_limit = np.power(10.,ahlers_limit_log)

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


	fig2 = plt.figure(figsize=(2.*11,8.5))
	ax2_aeff = fig2.add_subplot(1,2,1)
	ax2_lim = fig2.add_subplot(1,2,2)
	
	ax2_lim.plot(ahlers_energy,ahlers_limit,'-', linewidth=3.0,color='gray',label="Ahlers 2012")
	ax2_lim.fill_between(icecube_energy,icecube_thrumu_lower_efe,icecube_thrumu_upper_efe,facecolor='orange',alpha=0.3)
	ax2_lim.fill_between(icecube_energy,icecube_combined_lower_efe,icecube_combined_upper_efe,facecolor='magenta',alpha=0.3)	
	ax2_lim.plot(icecube_energy,icecube_thrumu_efe,'-.', linewidth=3.0,color='orange',label=r'IceCube Thru-Mu (E$^{-2.19}$)')
	ax2_lim.plot(icecube_energy,icecube_combined_efe,':', linewidth=3.0,color='magenta',label='IceCube Combined (E$^{-2.50}$)')

	ax2_aeff.plot(np.power(10.,total_energy_logeV),total_aeff,'-o',linewidth=2.0,color='blue',label=r'Total Deployed ARA Sensitivity')
	ax2_lim.plot(np.power(10.,total_energy_logeV),total_limit,'-o',linewidth=2.0,color='blue',label=r'Total Deployed ARA Sensitivity')
	beautify_aeff(ax2_aeff)
	beautify_limit_withtheory(ax2_lim,3)
	ax2_aeff.set_ylim([1e2,1e12])
	ax2_lim.set_ylim([1e-22,1e-10])
	fig2.savefig("current_total_limits_and_aeff.png",edgecolor='none',bbox_inches="tight") #save the figure




def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e2 #the lower x limit
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

#pass it the axes and the number of theory curves you're including
#always pass the theory curves first
def beautify_limit(this_ax):
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
	this_legend = this_ax.legend(loc='upper right',title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
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
