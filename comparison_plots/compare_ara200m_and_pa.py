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


	ara_logeV, ara_200m_aeff = detector.get_aeff('ara_200m_1year_fromlimit')
	phased_logeV, phased_aeff = detector.get_aeff('ara_phased_fromlimit')

	#trim the ARA sensitivity to match the phased array in energy binning
	ara_logeV = ara_logeV[:-1]
	ara_200m_aeff = ara_200m_aeff[:-1]

	app_ara = ara_200m_aeff * const.SecPerYear
	app_pa = phased_aeff * const.SecPerYear

	ara_app_interpolator = splrep(ara_logeV, np.log10(app_ara),k=4)
	pa_app_interpolator = splrep(phased_logeV, np.log10(app_pa),k=4)

	counts_ara_icecubethrumu=[]
	counts_ara_icecubecombined=[]
	counts_ara_ahlers=[]
	energy_bins=[]
	bins = np.arange(16,19.5,0.5)
	for bin in bins:
		temp_logev = np.arange(bin,bin+0.5,0.01)
		temp_energy = np.power(10.,temp_logev)
		#the ara sensitivity
		temp_ara_aeff = np.power(10.,splev(temp_logev, ara_app_interpolator))
		#the flux estimates
		temp_icecube_thrumu = icecube_thrumu_function(temp_energy)
		temp_icecube_combined = icecube_combined_function(temp_energy)
		temp_ahlers = np.power(10.,splev(temp_logev, ahlers_interpolator))/temp_energy
		#the counts
		temp_counts_thrumu = np.trapz(temp_icecube_thrumu*temp_ara_aeff,temp_energy)
		temp_counts_combined = np.trapz(temp_icecube_combined*temp_ara_aeff,temp_energy)
		temp_counts_ahlers = np.trapz(temp_ahlers*temp_ara_aeff,temp_energy)
		#push back the counts
		counts_ara_icecubethrumu.append(temp_counts_thrumu)
		counts_ara_icecubecombined.append(temp_counts_combined)
		counts_ara_ahlers.append(temp_counts_ahlers)
		#save the energy bin
		energy_bins.append(np.power(10.,bin))
	
	counts_ara_icecubethrumu = np.array(counts_ara_icecubethrumu)
	counts_ara_icecubecombined = np.array(counts_ara_icecubecombined)
	counts_ara_ahlers = np.array(counts_ara_ahlers)
	energy_bins = np.array(energy_bins)

	#what's the phased array contribution
	counts_pa_icecubethrumu=[]
	counts_pa_icecubecombined=[]
	counts_pa_ahlers=[]
	pa_energy_bins=[]
	bins = np.arange(16,19.5,0.5)
	for bin in bins:
		temp_logev = np.arange(bin,bin+0.5,0.01)
		temp_energy = np.power(10.,temp_logev)
		temp_pa_aeff = np.power(10.,splev(temp_logev, pa_app_interpolator))

		#the flux estimates
		temp_icecube_thrumu = icecube_thrumu_function(temp_energy)
		temp_icecube_combined = icecube_combined_function(temp_energy)
		temp_ahlers = np.power(10.,splev(temp_logev, ahlers_interpolator))/temp_energy

		temp_counts_thrumu_pa = np.trapz(temp_icecube_thrumu*temp_pa_aeff,temp_energy)
		temp_counts_combined_pa = np.trapz(temp_icecube_combined*temp_pa_aeff,temp_energy)
		temp_counts_ahlers_pa = np.trapz(temp_ahlers*temp_pa_aeff,temp_energy)

		counts_pa_icecubethrumu.append(temp_counts_thrumu_pa)
		counts_pa_icecubecombined.append(temp_counts_combined_pa)
		counts_pa_ahlers.append(temp_counts_ahlers_pa)

		pa_energy_bins.append(np.power(10.,bin))
	counts_pa_icecubethrumu = np.array(counts_pa_icecubethrumu)
	counts_pa_icecubecombined = np.array(counts_pa_icecubecombined)
	counts_pa_ahlers = np.array(counts_pa_ahlers)
	pa_energy_bins = np.array(pa_energy_bins)

	ratio_ahlers = counts_pa_ahlers/counts_ara_ahlers
	ratio_icecubethrumu = counts_pa_icecubethrumu/counts_ara_icecubethrumu
	ratio_icecubecombined = counts_pa_icecubecombined/counts_ara_icecubecombined

	fig_counts = plt.figure(figsize=(3.*11,8.5))
	ax_count_ara = fig_counts.add_subplot(1,3,1)
	ax_count_pa = fig_counts.add_subplot(1,3,2)
	ax_count_ratio = fig_counts.add_subplot(1,3,3)	
	ax_count_ara.set_title("ARA 200m",fontsize=24)
	ax_count_pa.set_title("ARA 200m w/ PA",fontsize=24)
	ax_count_ratio.set_title("Ratio ARA200/PA",fontsize=24)	
	
	n, bins, patches= ax_count_ara.hist(energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_ara_ahlers,
										label=r'Ahlers 2012: %.2f'%counts_ara_ahlers.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='red',
										linewidth=4)
	n, bins, patches= ax_count_ara.hist(energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_ara_icecubethrumu,
										label=r'IceCube Thru-Mu E$^{-2.19}$: %.2f'%counts_ara_icecubethrumu.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='blue',
										linewidth=4)
	n, bins, patches= ax_count_ara.hist(energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_ara_icecubecombined,
										label=r'IceCube Combined E$^{-2.5}$: %.3f'%counts_ara_icecubecombined.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='green',
										linewidth=4)
	
	n, bins, patches= ax_count_pa.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_pa_ahlers,
										label=r'Ahlers 2012: %.2f'%counts_pa_ahlers.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='red',
										linewidth=4)
	n, bins, patches= ax_count_pa.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_pa_icecubethrumu,
										label=r'IceCube Thru-Mu E$^{-2.19}$: %.2f'%counts_pa_icecubethrumu.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='blue',
										linewidth=4)
	n, bins, patches= ax_count_pa.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=counts_pa_icecubecombined,
										label=r'IceCube Combined E$^{-2.5}$: %.3f'%counts_pa_icecubecombined.sum(),
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='green',
										linewidth=4)
	
	n, bins, patches= ax_count_ratio.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=ratio_ahlers,
										label=r'Ahlers 2012',
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='red',
										linewidth=4)
	n, bins, patches= ax_count_ratio.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=ratio_icecubethrumu,
										label=r'IceCube Thru-Mu E$^{-2.19}$',
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='blue',
										linewidth=4)
	n, bins, patches= ax_count_ratio.hist(pa_energy_bins,
										bins=np.power(10.,np.arange(16,22,1)),
										weights=ratio_icecubecombined,
										label=r'IceCube Combined E$^{-2.5}$',
										fill=False, 
										stacked=True, 
										histtype='step', 
										edgecolor='green',
										linewidth=4)


	beautify_counts(ax_count_ara)
	beautify_counts(ax_count_pa)
	ratio_legend = beautify_counts(ax_count_ratio)
	ax_count_ara.set_ylim([0,0.2])
	ax_count_pa.set_ylim([0,0.2])
	ax_count_ratio.set_ylim([0,14])
	ax_count_ratio.set_ylabel('Ratio',size=20) #give it a title
	ratio_legend.set_title('')
	fig_counts.savefig("ara200m_vs_pa_counts.png",edgecolor='none',bbox_inches="tight") #save the figure


	# might plot this later
	# there is a really big gain at the 10e16 bin because of how quickly the sensitivity is falling...
	# ratio = phased_aeff/ara_200m_aeff
	# fig_ratio = plt.figure(figsize=(2.*11,8.5))
	# ax_ratio = fig_ratio.add_subplot(1,1,1)
	# ax_ratio.plot(np.power(10.,ara_logeV),ratio)
	# sizer=20
	# xlow = 1.e15 #the lower x limit
	# xup = 1.e21 #the uppper x limit
	# ax_ratio.set_xlabel('EnergZy  [eV]',size=sizer) #give it a title
	# ax_ratio.set_ylabel('Ratio of PA/ARA200m',size=sizer)
	# ax_ratio.set_xscale('log')
	# ax_ratio.tick_params(labelsize=sizer)
	# ax_ratio.set_xlim([xlow,xup]) #set the x limits of the plot
	# #ax_ratio.set_ylim([ylow,yup]) #set the y limits of the plot
	# ax_ratio.grid()
	# fig_ratio.savefig("ratio_ara200m_to_pa.png",edgecolor='none',bbox_inches="tight") #save the figure


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
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('Events Per Year',size=sizer)
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='upper left',title='Analysis Level Event Counts Per Year')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	return this_legend

main()
