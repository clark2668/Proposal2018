# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import *
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev
from pylab import setp

SecPerDay = 86400. #second in a year
KM2toCM2 = 1.e10 #multiply by this to convert from km^2 to cm^2 (divide to go the other way)

#plots the projected limit from a combination of ARA, ARIANNA, and PA stations
#plots limits and effective areas side by side
#also, correctly accounts for relative livetimes between an autonomous and cabled station

def main():

	'''
		What we want to do is extract an analysis level Effective Area * Steradians
		We are going to do this by backing it out from the published limit curves
		Here are the steps
		1) First, import the data from the limit curve
		2) Then, I rescale the limit curve so that it would represent 1 year
		3) Then, we hvae to interpolate the limit curve to every half decade, to go from logE = 15.5 to 20.5
		4) Then, we can convert from the limit curve to an effective area by multiplying by factors of 2.3, seconds, and converting km^2 to cm^2
	'''


	'''
		Below I comment a line-by-line example of how this works
	'''

	#import data
	ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	#extract the 'limit' column, and rescale to 1 year
	ara2_limit=ara2_data['limit'] * (445./365.) #take out livetime factor (includes 2 station correction)
	#convert to an effective area in units of km^2
	ara2_aeff = 2.3/ara2_limit/(365.*SecPerDay)/KM2toCM2
	#these are the energies for the limit curve in logE
	ara2_energy_logev = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	#which we can unlog to just get the energy
	ara2_energy = np.power(10.,ara2_energy_logev)

	#now, we can set up an interpolator to interpolate across the curve
	#we do this in log-log space, because it's easier to do that interpolation
	ara2_interpolator = splrep(ara2_energy_logev, np.log10(ara2_limit),k=2)
	#we want to define the log(energies) we want have interpolated limits at
	ara2_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	#convert from logE to E
	ara2_energy_interp = np.power(10.,ara2_energy_logev_interp)
	#actually evaluate the spline interpolator at logE
	ara2_limit_interp = np.power(10.,splev(ara2_energy_logev_interp, ara2_interpolator))
	#this will be the actual ara2_aeff at the interpolated energies
	ara2_aeff_interp = 2.3/ara2_limit_interp/(365.*SecPerDay)/KM2toCM2

	'''
		To estimate the sensitivity of the phased array, we are going to cheat just a touch
		We are going to cheat in the sense that we're just going to slide the ARA2 curve to the left
		on an EF(E) plot by a factor of two. Then we repeat the whole business above of interpolating and converting.
	'''

	pa_limit=ara2_limit
	pa_energy = ara2_energy/2 #scoot the curve over
	pa_energy_logev = np.log10(pa_energy)
	pa_aeff = 2.3/pa_limit/(365.*SecPerDay)/KM2toCM2

	pa_interpolator = splrep(pa_energy_logev, np.log10(pa_limit),k=1)
	pa_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	pa_energy_interp = np.power(10.,pa_energy_logev_interp)
	pa_limit_interp = np.power(10.,splev(pa_energy_logev_interp, pa_interpolator))
	pa_aeff_interp = 2.3/pa_limit_interp/(365.*SecPerDay)/KM2toCM2

	testbed_data = np.genfromtxt("limits/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	testbed_limit=testbed_data['limit'][0:-1] * (415./365.) #strip off last point to match ARIANNA, and take out livetime
	testbed_aeff = 2.3/testbed_limit/(365.*SecPerDay)/KM2toCM2
	testbed_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
	testbed_energy = np.power(10.,testbed_energy_logev)

	testbed_interpolator = splrep(testbed_energy_logev[1:], np.log10(testbed_limit[1:]),k=3)
	testbed_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	testbed_energy_interp = np.power(10.,testbed_energy_logev_interp)
	testbed_limit_interp = np.power(10.,splev(testbed_energy_logev_interp, testbed_interpolator))
	testbed_aeff_interp = 2.3/testbed_limit_interp/(365.*SecPerDay)/KM2toCM2	

	hra3_data = np.genfromtxt("limits/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	hra3_limit=hra3_data['limit']/(hra3_data['energy']/1e9) * (170./365.)#need to divide out the GeV, and take out livetime (includes 3 station correction)
	hra3_aeff = 2.3/hra3_limit/(365.*SecPerDay)/KM2toCM2
	hra3_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
	hra3_energy =  np.power(10.,hra3_energy_logev)

	arianna_interpolator = splrep(hra3_energy_logev[1:], np.log10(hra3_limit[1:]),k=3)
	arianna_energy_logev_interp = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	arianna_energy_interp = np.power(10.,arianna_energy_logev_interp)
	arianna_limit_interp = np.power(10.,splev(arianna_energy_logev_interp, arianna_interpolator))
	arianna_aeff_interp = 2.3/arianna_limit_interp/(365.*SecPerDay)/KM2toCM2

	
	'''
		Okay, now we need to import the IceCube measurement
		They have two flux measurements, one from their "thru-mu" and one from "combind-fit"
	'''


	icecube_energy_logev=np.array([15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5])
	icecube_energy = np.power(10.,icecube_energy_logev)

	#E is in eV
	#returns number/eV/cm^2/s/sr

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
		Okay, now we need to pull in some theoretical limits
		We'll do ahlers and kotera
	'''

	ahlers_data = np.genfromtxt("limits/ahlers_2012.csv",delimiter=',',skip_header=1,names=['energy','flux'])
	ahlers_energy_logev = ahlers_data['energy']
	ahlers_energy = np.power(10.,ahlers_energy_logev)
	ahlers_limit_log = ahlers_data['flux']
	ahlers_limit = np.power(10.,ahlers_limit_log) 

	kotera_max_data = np.genfromtxt("limits/kotera_max.csv",delimiter=',',skip_header=1,names=['energy','flux'])
	kotera_max_energy_logev = kotera_max_data['energy']
	kotera_max_energy = np.power(10.,kotera_max_energy_logev)
	kotera_max_limit_log = kotera_max_data['flux']
	kotera_max_limit = np.power(10.,kotera_max_limit_log) 

	kotera_min_data = np.genfromtxt("limits/kotera_min.csv",delimiter=',',skip_header=1,names=['energy','flux'])
	kotera_min_energy_logev = kotera_min_data['energy']
	kotera_min_energy = np.power(10.,kotera_min_energy_logev)
	kotera_min_limit_log = kotera_min_data['flux']
	kotera_min_limit = np.power(10.,kotera_min_limit_log)

	#to make a fill between for the kotera, we need to interpolate the min data to the max data energy points
	kotera_min_interpolator = splrep( kotera_min_energy_logev , kotera_min_limit_log,k=1)
	kotera_min_limit_interp = np.power(10.,splev(kotera_max_energy_logev,kotera_min_interpolator))

	'''
		Okay, now we can design our new "ideal" detectors
		We will want it to be some combination of ARA stations "num_ara"
		phased array enhances ARA stations "num_pa"
		arianna stations "num_arianna"
		and number of years we want them to run for "num_years"
	'''

	num_ara=5.
	num_pa=0
	num_arianna=50.
	num_years=5.

	hybrid_energy_logev = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	hybrid_energy = np.power(10.,hybrid_energy_logev)

	'''
		We want to break up our calculation into two pieces
		An "autonomous" piece, which is composed of just ARIANNA stations
		An "cabled" piece, which is composed of ARA + PA stations
	'''
	
	#so these are the effective areas in km^2 * sr
	#this is the part that sets the total cabled aeff
	hybrid_aeff_cabled = (num_ara*ara2_aeff_interp) + (num_pa*pa_aeff_interp)	
	#this is the part that sets the total autonomous aeff
	hybrid_aeff_autonomous = (num_arianna*arianna_aeff_interp)
	#this is the part that sets the total aeff	
	hybrid_aeff = hybrid_aeff_cabled + hybrid_aeff_autonomous
	
	#to compute the expected limit, we now need to take into accountt the *livetimes*, convert to seconds, and to cm^2
	hybrid_aeff_cabled_with_livetime = hybrid_aeff_cabled*(num_years*365.*SecPerDay)*KM2toCM2
	hybrid_aeff_autonomous_with_livetime = hybrid_aeff_autonomous*(num_years*365.*SecPerDay*0.6)*KM2toCM2
	hybrid_aeff_with_livetime = hybrid_aeff_cabled_with_livetime + hybrid_aeff_autonomous_with_livetime
	
	#now we can actually compute the limit setting power
	hybrid_cabled_limit = 1./hybrid_aeff_cabled_with_livetime
	hybird_autonomous_limit = 1./hybrid_aeff_autonomous_with_livetime
	hybrid_limit = 1./hybrid_aeff_with_livetime

	
	'''
		Now we make plots
	'''

	fig = plt.figure(figsize=(2*11,8.5)) #make a figure object
	ax_efe = fig.add_subplot(1,2,1) #make a subplot for the limit
	ax_aeff = fig.add_subplot(1,2,2) #make a subplot for the effective areas

	ax_efe.plot(ahlers_energy,ahlers_limit,'-', linewidth=3.0,color='gray',label="Ahlers 2012")
	#ax_efe.fill_between(kotera_max_energy,kotera_min_limit_interp,kotera_max_limit,facecolor='lightgray',label="Kotera")
	ax_efe.fill_between(icecube_energy,icecube_thrumu_lower_efe,icecube_thrumu_upper_efe,facecolor='orange',alpha=0.3)
	ax_efe.fill_between(icecube_energy,icecube_combined_lower_efe,icecube_combined_upper_efe,facecolor='magenta',alpha=0.3)	
	ax_efe.plot(icecube_energy,icecube_thrumu_efe,'-.', linewidth=3.0,color='orange',label=r'IceCube Thru-Mu (E$^{-2.19}$)')
	ax_efe.plot(icecube_energy,icecube_combined_efe,':', linewidth=3.0,color='magenta',label='IceCube Combined (E$^{-2.50}$)')	

	#how many theory curves did we add above?
	num_theory=3

	#ax_efe.plot(arianna_energy_interp,arianna_limit_interp/0.68,'--^', linewidth=2.0,color='blue',label=r'0.6 Yr $\times$ 1 ARIANNA')
	#ax_efe.plot(ara2_energy_interp,ara2_limit_interp,'--v', linewidth=2.0,color='green',label=r'1 Yr $\times$ 1 ARA 200m')	
	auto_plot = ax_efe.plot(hybrid_energy,hybird_autonomous_limit,'--^', linewidth=3.0,color='blue',label=r'Auto: (0.6 $\times$ %d Yrs) $\times$ (%d ARIANNA)'%(num_years,num_arianna))
	ax_efe.plot(hybrid_energy,hybrid_cabled_limit,'--v', linewidth=3.0,color='green',label=r'Cable: (%d Yrs) $\times$ (%d ARA + %d PA)'%(num_years,num_ara,num_pa))	
	ax_efe.plot(hybrid_energy,hybrid_limit,'-s', linewidth=3.0,color='black',label='Hybrid: %d ARIANNA + %d ARA + %d PA'%(num_arianna,num_ara,num_pa))	

	beautify_limit(ax_efe,num_theory)

	ax_aeff.plot(hybrid_energy,hybrid_aeff_autonomous,'--^', linewidth=3.0,color='blue',label='Auto: %d ARIANNA'%num_arianna)
	ax_aeff.plot(hybrid_energy,hybrid_aeff_cabled,'--v', linewidth=3.0,color='green',label='Cable: %d ARA + %d PA'%(num_ara,num_pa))
	ax_aeff.plot(hybrid_energy,hybrid_aeff,'-s', linewidth=3.0,color='black',label='Hybrid: %d ARIANNA + %d ARA + %d PA'%(num_arianna,num_ara,num_pa))
	beautify_aeff(ax_aeff)

	fig.savefig("projected_hybrid.png",edgecolor='none',bbox_inches="tight") #save the figure

	'''
		The following two functions just do a bunch of matplotlib beautification things
	'''

def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e-8 #the lower x limit
	yup = 1.e1 #the uppper x limit
	this_ax.set_xlabel('Energy  [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('Analysis Level [A$\Omega]_{eff}$  [km$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_legend = this_ax.legend(loc='lower right',title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	this_ax.set_yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])

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
	legend2 = this_ax.legend(handles[num_theory:], labels[num_theory:], loc='upper right',title='Single Event Sensitivity (Analysis Level)')
	setp(legend1.get_texts(), fontsize=17)
	setp(legend1.get_title(), fontsize=17)
	setp(legend2.get_texts(), fontsize=17)
	setp(legend2.get_title(), fontsize=17)
	this_ax.grid()

main()