# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from pylab import setp
import detector_data as detector
import constants as const

def main():

	icecube_logeV, icecube_limit = detector.get_limit('icecube_2018')
	anita_logeV, anita_limit = detector.get_limit('anita_2018')

	icecube_exposure = 2.44/np.log(10)/0.5/icecube_limit
	anita_exposure = 2.44/np.log(10)/0.5/anita_limit
	
	target_exposure = icecube_exposure[:8]
	target_energy = icecube_logeV[:8]

	#this is the exposure we need to get our hands on
	target_exposure = np.append(target_exposure,anita_exposure[4:5])
	target_energy = np.append(target_energy,20)
	
	fig_tobeat = plt.figure(figsize=(2*11,8.5))
	ax_limit_tobeat = fig_tobeat.add_subplot(1,2,1)
	ax_exposure_tobeat = fig_tobeat.add_subplot(1,2,2)

	ax_limit_tobeat.plot(np.power(10.,icecube_logeV),icecube_limit,'-^', linewidth=2.0,color='blue',label=r'IceCube 2018',markersize=10)
	ax_limit_tobeat.plot(np.power(10.,anita_logeV),anita_limit,'-v', linewidth=2.0,color='green',label=r'ANITA 2018',markersize=10)	
	
	ax_exposure_tobeat.plot(np.power(10.,icecube_logeV),icecube_exposure,'-^', linewidth=2.0,color='blue',label=r'IceCube 2018',markersize=10)
	ax_exposure_tobeat.plot(np.power(10.,anita_logeV),anita_exposure,'-v', linewidth=2.0,color='green',label=r'ANITA 2018',markersize=10)
	ax_exposure_tobeat.plot(np.power(10.,target_energy),target_exposure,'--', linewidth=8.0,color='red',label=r'Minimum Target Exposure for NGRA')

	beautify_limit(ax_limit_tobeat)
	beautify_aeff(ax_exposure_tobeat)
	fig_tobeat.savefig("limit_to_beat.png",edgecolor='none',bbox_inches="tight") #save the figure


	total_cost_1_5km=[]
	total_cost_1km=[]
	fixed_costs=[]
	station_number = np.arange(0,100,1)
	for station in station_number:
		total_cost_1_5km.append(cost(1.5,station))
		total_cost_1km.append(cost(1,station))
		fixed_costs.append(543000)

	total_cost_1_5km = np.array(total_cost_1_5km)
	total_cost_1km = np.array(total_cost_1km)
	fixed_costs = np.array(fixed_costs)
	average_cost_1_5km = total_cost_1_5km[1:]/station_number[1:]
	average_cost_1km = total_cost_1km[1:]/station_number[1:]


	num_cabled = 100.
	cost_to_cable = cost(1.5, num_cabled)
	num_extra_arianna = (cost_to_cable/15000.)
	# print "Cost to cable: ", cost_to_cable
	print "Cost to cable: %.2f and cost per cabled station %.2f "%(cost_to_cable,cost_to_cable/num_cabled)
	print "Worth %.2f ARIANNA Stations "%(num_extra_arianna)


	livetime = 5.*const.SecPerYear
	arianna_energy, arianna_aeff = detector.get_aeff('arianna_icrc2017_fromlimit')
	arianna_energy=arianna_energy[1:]
	arianna_energy=arianna_energy[:-2]
	arianna_aeff=arianna_aeff[1:]
	arianna_aeff=arianna_aeff[:-2]
	pa_energy, pa_eff = detector.get_aeff('ara_phased_fromlimit')

	ara_exposure = pa_eff * 5.*const.SecPerYear
	arianna_exposure = arianna_aeff * 5.*const.SecPerYear * 0.58 #arianna only gets 58% uptime
	num_auto = 800.
	num_auto_only = num_extra_arianna+num_auto

	cabled_contribution = (ara_exposure*num_cabled)
	auto_contribution = (arianna_exposure*num_auto)
	auto_only_exposure = arianna_exposure * num_auto_only

	#merged detector
	total = cabled_contribution + auto_contribution

	fig_thiscase = plt.figure(figsize=(11,8.5))
	ax_thiscase = fig_thiscase.add_subplot(1,1,1)

	ax_thiscase.plot(np.power(10.,target_energy),target_exposure,'--', linewidth=4.0,color='red',label=r'Minimum Target Exposure for NGRA')
	ax_thiscase.plot(np.power(10.,pa_energy),total,'--', linewidth=4.0,color='black',label=r'5 years, %d Cabled + %d Autonomous'%(num_cabled,num_auto))

	ax_thiscase.plot(np.power(10.,arianna_energy),auto_contribution,'-v', linewidth=4.0,color='green',label=r'Mixed: %d Auto Contribution'%num_auto,markersize=10)	
	ax_thiscase.plot(np.power(10.,pa_energy),cabled_contribution,'-^', linewidth=4.0,color='blue',label=r'Mixed: %d Cabled Contribution'%num_cabled,markersize=10)
	ax_thiscase.plot(np.power(10.,arianna_energy),auto_only_exposure,'-v', linewidth=2.0,color='magenta',label=r'Auto Only: %d Station'%num_auto_only,markersize=10)	
	
	beautify_aeff(ax_thiscase)
	fig_thiscase.savefig("%dcabled_case.png"%num_cabled,edgecolor='none',bbox_inches="tight") #save the figure



	fig = plt.figure(figsize=(2*11,8.5))
	ax_total = fig.add_subplot(1,2,1)
	ax_average = fig.add_subplot(1,2,2)

	ax_total.plot(station_number,total_cost_1_5km/1000.,label='1.5 km Spacing',linewidth=4)
	ax_total.plot(station_number,total_cost_1km/1000.,label='1 km Spacing',linewidth=4)

	ax_average.plot(station_number[1:],average_cost_1_5km/1000.,label='1.5 km Spacing',linewidth=4)
	ax_average.plot(station_number[1:],average_cost_1km/1000.,label='1 km Spacing',linewidth=4)
	
	ax_total.set_title("Total Cabling Cost",fontsize=24)
	ax_average.set_title("Cabling Cost Per Station",fontsize=24)	

	beautify(ax_total)
	beautify(ax_average)
	sizer=20
	ax_total.set_ylabel('Thousands of USD',size=sizer)
	ax_average.set_ylabel('Thousands of USD',size=sizer)
	ax_average.set_ylim([0,50]) #set the y limits of the plot

	fig.savefig("cost.png",edgecolor='none',bbox_inches="tight") #save the figure

def cost(spacing, num):
	fixed_cost = 543000
	
	cost_per_km = 3800+700+125
	cable_cost_per_station = cost_per_km * spacing
	electronics_per_station=75+40
	cost_per_station = (electronics_per_station+cable_cost_per_station)
	scaled_cost = cost_per_station * num
	
	num_bundles = (num/10)+1
	bundle_cost = (1385 + cable_cost_per_station) * num_bundles

	total = fixed_cost + scaled_cost + bundle_cost
	return total


def beautify(this_ax):
	sizer=20
	this_ax.set_xlabel('Number of Cabled Stations',size=sizer) #give it a title
	# this_ax.set_yscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.grid()
	this_legend = this_ax.legend(loc='upper left')#,title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=20)
	setp(this_legend.get_title(), fontsize=20)

#pass it the axes and the number of theory curves you're including
#always pass the theory curves first
def beautify_limit(this_ax):
	sizer=20
	xlow =  1e16 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-19 #the lower y limit
	yup = 1e-12 #the upper y limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('E F(E) [$cm^{-2} s^{-1} sr^{-1}$]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_legend = this_ax.legend(loc='upper right')#,title='Analysis Level')
	setp(this_legend.get_texts(), fontsize=17)
	setp(this_legend.get_title(), fontsize=17)
	this_ax.grid()

def beautify_aeff(this_ax):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e13 #the lower x limit
	yup = 1.e20 #the uppper x limit
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

main()