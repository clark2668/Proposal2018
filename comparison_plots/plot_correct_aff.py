# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp


SecPerDay = 86400. #second in a year
KM2toCM2 = 1.e10 #multiply by this to convert from km^2 to cm^2 (divide to go the other way)

def main():

	ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	ara2_limit=ara2_data['limit']
	ara2_limit_1yr = ara2_limit * (2. * 225. / 365.)
	ara2_aperture = 2.44/np.log(10)/0.5/ara2_limit/(2. * 225.*SecPerDay)
	ara2_energy = ara2_data['energy']
	ara2_energy_logev = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])

	# 2.44 is statistical factor
	# ln(10) is conversion from linear to logarithmic bins
	# 0.5 is bin size in log energy space
	# #445. * SecPerDay is removal of time in seconds from limit
	# KM2toCM2 converts cm^2 to km^2

	ara2_interpolator = splrep(ara2_energy_logev, np.log10(ara2_limit),k=2)
	ara2_energy_logev_interp = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	ara2_energy_interp = np.power(10.,ara2_energy_logev_interp)
	ara2_limit_interp = np.power(10.,splev(ara2_energy_logev_interp, ara2_interpolator))
	ara2_aeff_interp = 2.44/np.log(10)/0.5/ara2_limit_interp/(365.*SecPerDay)

	icecube_energy_logev=np.array([16,16.5,17,17.5,18,18.5,19,19.5,20,20.5])
	icecube_energy = np.power(10.,icecube_energy_logev)
	def icecube_thrumu_function(E):
		return 3.03 * ((E/1.e14)**-2.19) * 1e-27
	icecube_thrumu_efe = icecube_thrumu_function(icecube_energy)

	livetime = 37.*365.*86400.

	energy_for_integral_logev = np.arange(16,20.5,0.01) #finely spaced grid
	energy_for_integral = np.power(10.,energy_for_integral_logev)
	icecube_thrumu_efe_for_integral = icecube_thrumu_function(energy_for_integral)
	ara2_aeff_interp_for_integral = 2.44/np.log(10)/0.5/np.power(10.,splev(energy_for_integral_logev, ara2_interpolator))/(365.*SecPerDay)

	counts = np.trapz(icecube_thrumu_efe*ara2_aeff_interp*livetime,icecube_energy)
	#print "Counts Simple Try: ",counts

	counts = np.trapz(icecube_thrumu_efe_for_integral*ara2_aeff_interp_for_integral*livetime,energy_for_integral)
	print "Counts Finer Binning Try: ",counts

	counts_try2=[]
	bins = np.arange(16,21,0.5)
	for bin in bins:
		temp_logev = np.arange(bin,bin+0.5,0.01)
		temp_energy = np.power(10.,temp_logev)
		temp_icecube = icecube_thrumu_function(temp_energy)
		temp_ara2_aeff = 2.44/np.log(10)/0.5/np.power(10.,splev(temp_logev, ara2_interpolator))/(365.*SecPerDay)
		temp_counts = np.trapz(temp_icecube*temp_ara2_aeff*livetime,temp_energy)
		counts_try2.append(temp_counts)
	print "Counts Multi-Bin Method: ",counts


	fig = plt.figure(figsize=(2.*11,8.5))
	ax_limit = fig.add_subplot(1,2,1)
	ax_aff = fig.add_subplot(1,2,2)

	ax_limit.plot(icecube_energy,icecube_thrumu_efe,'-.', linewidth=3.0,color='orange',label=r'IceCube Thru-Mu (E$^{-2.19}$)')
	ax_limit.plot(ara2_energy,ara2_limit,'-s', linewidth=2.0,color='blue',label=r'225 Days of 2 ARA 200m Stations')
	ax_limit.plot(ara2_energy,ara2_limit_1yr,'--s', linewidth=2.0,color='blue',label=r'1 Yr of 1 ARA 200m Staion')
	beautify_limit(ax_limit, 1)

	ax_aff.plot(ara2_data['energy'],ara2_aperture,'-s', linewidth=2.0,color='blue',label="ARA 200m Station")
	beautify_aeff(ax_aff)

	fig.savefig("limit_and_aff.png",edgecolor='none',bbox_inches="tight") #save the figure

def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e3 #the lower x limit
	yup = 1.e10 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
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
	handles, labels = this_ax.get_legend_handles_labels()
	legend1 = this_ax.legend(handles[0:num_theory], labels[0:num_theory], loc='lower left',title='Flux Models')
	this_ax.add_artist(legend1)
	legend2 = this_ax.legend(handles[num_theory:], labels[num_theory:], loc='upper right',title='Limit (Analysis Level)')
	setp(legend1.get_texts(), fontsize=17)
	setp(legend1.get_title(), fontsize=17)
	setp(legend2.get_texts(), fontsize=17)
	setp(legend2.get_title(), fontsize=17)
	this_ax.grid()

main()
