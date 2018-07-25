# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d, splrep, splev

SecPerDay = 86400. #second in a year
KM2toCM2 = 1.e10 #multiply by this to convert from km^2 to cm^2 (divide to go the other way)

def main():

	ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	ara2_limit=ara2_data['limit'] * (445./365.) #take out livetime factor (includes 2 station correction)
	ara2_aperture = 2.3/ara2_limit/(365.*SecPerDay)/KM2toCM2
	ara2_energy_logev = np.array([16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	ara2_energy = np.power(10.,ara2_energy_logev)

	ara2_interpolator = splrep(ara2_energy_logev, np.log10(ara2_limit),k=2)
	ara2_energy_logev_interp = np.insert(ara2_energy_logev,0,15.5,axis=0)
	ara2_energy_interp = np.power(10.,ara2_energy_logev_interp)
	ara2_limit_interp = np.power(10.,splev(ara2_energy_logev_interp, ara2_interpolator))
	ara2_aperture_interp = 2.3/ara2_limit_interp/(365.*SecPerDay)/KM2toCM2	

	pa_limit=ara2_limit
	pa_energy = ara2_energy/2
	pa_energy_logev = np.log10(pa_energy)
	pa_aperature = 2.3/pa_limit/(365.*SecPerDay)/KM2toCM2

	pa_interpolator = splrep(pa_energy_logev, np.log10(pa_limit),k=1)
	pa_energy_logev_interp = np.insert(ara2_energy_logev,0,15.5,axis=0)
	pa_energy_interp = np.power(10.,pa_energy_logev_interp)
	pa_limit_interp = np.power(10.,splev(pa_energy_logev_interp, pa_interpolator))
	pa_aperture_interp = 2.3/pa_limit_interp/(365.*SecPerDay)/KM2toCM2

	testbed_data = np.genfromtxt("limits/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	testbed_limit=testbed_data['limit'][0:-1] * (415./365.) #strip off last point to match ARIANNA, and take out livetime
	testbed_aperture = 2.3/testbed_limit/(365.*SecPerDay)/KM2toCM2
	testbed_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
	testbed_energy = np.power(10.,testbed_energy_logev)

	testbed_interpolator = splrep(testbed_energy_logev[1:], np.log10(testbed_limit[1:]),k=3)
	testbed_energy_logev_interp = np.insert(testbed_energy_logev,0,[15.5,16,16.5],axis=0)
	testbed_energy_interp = np.power(10.,testbed_energy_logev_interp)
	testbed_limit_interp = np.power(10.,splev(testbed_energy_logev_interp, testbed_interpolator))
	testbed_aperture_interp = 2.3/testbed_limit_interp/(365.*SecPerDay)/KM2toCM2	

	hra3_data = np.genfromtxt("limits/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	hra3_limit=hra3_data['limit']/(hra3_data['energy']/1e9) * (170./365.)#need to divide out the GeV, and take out livetime (includes 3 station correction)
	hra3_aperture = 2.3/hra3_limit/(365.*SecPerDay)/KM2toCM2
	hra3_energy_logev = np.array([17,17.5,18.,18.5,19,19.5,20,20.5])
	hra3_energy =  np.power(10.,hra3_energy_logev)

	arianna_interpolator = splrep(hra3_energy_logev[1:], np.log10(hra3_limit[1:]),k=3)
	arianna_energy_logev_interp = np.insert(hra3_energy_logev,0,[15.5,16,16.5],axis=0)
	arianna_energy_interp = np.power(10.,arianna_energy_logev_interp)
	arianna_limit_interp = np.power(10.,splev(arianna_energy_logev_interp, arianna_interpolator))
	arianna_aperture_interp = 2.3/arianna_limit_interp/(365.*SecPerDay)/KM2toCM2

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

	#to make a fill between, we need to interpolate the min data to the max data points
	kotera_min_interpolator = splrep( kotera_min_energy_logev , kotera_min_limit_log,k=1)
	kotera_min_limit_interp = np.power(10.,splev(kotera_max_energy_logev,kotera_min_interpolator))


	fig_efe = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_efe = fig_efe.add_subplot(1,1,1) #make a subplot
	ax_efe.plot(testbed_energy_interp,testbed_limit_interp,'-v', linewidth=2.0,color='red',label=r'1 Yr, 1 ARA 30m Testbed (Analysis Level)')
	ax_efe.plot(arianna_energy_interp,arianna_limit_interp,'-^', linewidth=2.0,color='green',label=r'1 Yr, 1 ARIANNA Station (Analysis Level)')
	ax_efe.plot(ara2_energy_interp,ara2_limit_interp,'-s', linewidth=2.0,color='blue',label=r'1 Yr, 1 ARA 200m Station (Analysis Level)')	
	ax_efe.plot(pa_energy_interp,pa_limit_interp,'-o', linewidth=2.0,color='black',label=r'1 Yr, 1 ARA 200m Station + PA ("projected")')
	
	ax_efe.plot(ahlers_energy,ahlers_limit,'--', linewidth=2.0,color='black',label=r'GZK: Ahlers 2012')
	ax_efe.fill_between(kotera_max_energy,kotera_min_limit_interp,kotera_max_limit,facecolor='lightgray',label='GZK: Kotera 2010')

	save_limit(fig_efe,ax_efe,"base_units_limit")

	fig_aff = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_aff = fig_aff.add_subplot(1,1,1) #make a subplot
	ax_aff.plot(pa_energy_interp,pa_aperture_interp,'-o', linewidth=2.0,color='black',label='ARA 200m Station + PA ("projected")')
	ax_aff.plot(ara2_energy_interp,ara2_aperture_interp,'-s', linewidth=2.0,color='blue',label='ARA 200m Station (Analysis Level)')
	ax_aff.plot(arianna_energy_interp,arianna_aperture_interp,'-^', linewidth=2.0,color='green',label='ARIANNA Station (Analysis Level)')
	ax_aff.plot(testbed_energy_interp,testbed_aperture_interp,'-v', linewidth=2.0,color='red',label='ARA 30m Testbed Station (Analysis Level)')
	save_aff(fig_aff,ax_aff,"base_units_aff")

	num_ara=5.
	num_pa=0
	num_arianna=100.
	num_years=5.

	num_ara2=0.
	num_pa2=5.
	num_arianna2=100.
	num_years2=3.

	hybrid_energy_logev = np.array([15.5,16.,16.5,17.,17.5,18.,18.5,19.,19.5,20,20.5])
	hybrid_energy = np.power(10.,hybrid_energy_logev)
	hybrid_aperature = (num_ara*ara2_aperture_interp) + (num_arianna*arianna_aperture_interp) + (num_pa*pa_aperture_interp)
	hybrid_limit = 2.3/hybrid_aperature/(num_years*365.*SecPerDay)/KM2toCM2
	
	hybrid_aperature2 = (num_ara2*ara2_aperture_interp) + (num_arianna2*arianna_aperture_interp) + (num_pa2*pa_aperture_interp)
	hybrid_limit2 = 2.3/hybrid_aperature2/(num_years2*365.*SecPerDay)/KM2toCM2

	fig_efe_hybrid = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_efe_hybrid = fig_efe_hybrid.add_subplot(1,1,1) #make a subplot
	ax_efe_hybrid.plot(hybrid_energy,hybrid_limit,'-o', linewidth=2.0,color='black',label='Hybrid: %d Yr, %d ARA + %d ARA/PA + %d ARIANNA'%(num_years,num_ara,num_pa,num_arianna))
	ax_efe_hybrid.plot(hybrid_energy,hybrid_limit2,'-o', linewidth=2.0,color='red',label='Hybrid: %d Yr, %d ARA + %d ARA/PA + %d ARIANNA'%(num_years2,num_ara2,num_pa2,num_arianna2))
	ax_efe_hybrid.plot(ara2_energy_interp,ara2_limit_interp,'--', linewidth=2.0,color='blue',label=r'1 Yr, 1 ARA 200m Station (Analysis Level)')
	ax_efe_hybrid.plot(arianna_energy_interp,arianna_limit_interp,'--', linewidth=2.0,color='green',label=r'1 Yr, 1 ARIANNA Station (Analysis Level)')
	ax_efe_hybrid.plot(ahlers_energy,ahlers_limit,'--', linewidth=2.0,color='black',label=r'GZK: Ahlers 2012')
	ax_efe_hybrid.fill_between(kotera_max_energy,kotera_min_limit_interp,kotera_max_limit,facecolor='lightgray',label='GZK: Kotera 2012')
	save_limit(fig_efe_hybrid,ax_efe_hybrid,"hybrid")


def save_aff(this_fig, this_ax, title):
	sizer=20
	xlow = 1.e14 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e-8 #the lower x limit
	yup = 1.e1 #the uppper x limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('A$_{eff}$ $\Omega$ [km$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_ax.legend(loc='upper left')
	this_ax.set_yticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1])
	
	this_fig.savefig(title+".pdf",edgecolor='none',bbox_inches="tight") #save the figure
	this_fig.savefig(title+".png",edgecolor='none',bbox_inches="tight") #save the figure


def save_limit(this_fig, this_ax,title):
	sizer=20
	xlow =  1.e14 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-20 #the lower y limit
	yup = 1e-9 #the upper y limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('E F(E) [$cm^{-2} s^{-1} sr^{-1}$]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.legend(loc='lower left')
	this_ax.grid()
	this_fig.savefig(title+".pdf",edgecolor='none',bbox_inches="tight") #save the figure
	this_fig.savefig(title+".png",edgecolor='none',bbox_inches="tight") #save the figure

main()