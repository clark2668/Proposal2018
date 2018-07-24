# -*- coding: utf-8 -*-
#from __future__ import unicode_literals
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
#import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math
from scipy.interpolate import interp1d

SecPerDay = 86400. #second in a year
KM2toCM2 = 1.e10 #multiply by this to convert from km^2 to cm^2 (divide to go the other way)

def main():

	ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	ara2_limit=ara2_data['limit']
	ara2_aperture = 2.3/ara2_limit/(445.*SecPerDay)/KM2toCM2/2

	testbed_data = np.genfromtxt("limits/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	testbed_limit=testbed_data['limit']
	testbed_aperture = 2.3/testbed_limit/(415.*SecPerDay)/KM2toCM2

	hra3_data = np.genfromtxt("limits/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	hra3_limit=hra3_data['limit']/(hra3_data['energy']/1e9) #need to divide out the GeV
	hra3_aperture = 2.3/hra3_limit/(170.*SecPerDay)/KM2toCM2/3

	fig_efe = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_efe = fig_efe.add_subplot(1,1,1) #make a subplot
	ax_efe.plot(testbed_data['energy'],testbed_data['limit'],'-v', linewidth=2.0,color='red',label=r'415 Days of 1 ARA 30m Testbed')
	ax_efe.plot(hra3_data['energy'],hra3_limit,'-^', linewidth=2.0,color='green',label=r'170 Days of 3 ARIANNA Stations')
	ax_efe.plot(ara2_data['energy'],ara2_data['limit'],'-s', linewidth=2.0,color='blue',label=r'222 Days of 2 ARA 200m Stations')
	save_limit(fig_efe,ax_efe,"current_limit")

	fig_aff = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_aff = fig_aff.add_subplot(1,1,1) #make a subplot
	ax_aff_diffu = ax_aff.twinx()
	ax_aff.plot(ara2_data['energy'],ara2_aperture,'-s', linewidth=2.0,color='blue',label="ARA 200m Station (Analysis Level)")
	ax_aff.plot(hra3_data['energy'],hra3_aperture,'-^', linewidth=2.0,color='green',label="ARIANNA Station (Analysis Level)")
	ax_aff.plot(testbed_data['energy'],testbed_aperture,'-v', linewidth=2.0,color='red',label="ARA 30 Testbed Station (Analysis Level)")
	ax_aff_diffu.plot(testbed_data['energy'],testbed_aperture*1.e10,'-v', linewidth=2.0,color='red')
	save_aff(fig_aff,ax_aff,ax_aff_diffu)\


def save_aff(this_fig, this_ax, this_ax2):
	sizer=20
	xlow = 1.e15 #the lower x limit
	xup = 1.e21 #the uppper x limit
	ylow = 1.e-7 #the lower x limit
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
	
	this_ax2.set_ylabel('A$_{eff}$ $\Omega$ [cm$^2$sr]',size=sizer)
	this_ax2.set_yscale('log')
	this_ax2.set_xscale('log')
	this_ax2.tick_params(labelsize=sizer)
	this_ax2.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax2.set_ylim([ylow*KM2toCM2,yup*KM2toCM2]) #set the y limits of the plot
	this_ax2.set_yticks([1e-7*KM2toCM2,1e-6*KM2toCM2,1e-5*KM2toCM2,1e-4*KM2toCM2,1e-3*KM2toCM2,1e-2*KM2toCM2,1e-1*KM2toCM2,1e0*KM2toCM2,1e1*KM2toCM2])
	
	this_fig.savefig('current_aff.pdf',edgecolor='none',bbox_inches="tight") #save the figure
	this_fig.savefig('current_aff.png',edgecolor='none',bbox_inches="tight") #save the figure



def save_limit(this_fig, this_ax,title):
	xlow =  1.00e14 #the lower x limit
	xup = 1e21 #the uppper x limit
	ylow = 1e-17 #the lower y limit
	yup = 1e-10 #the upper y limit
	sizer=20
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('E F(E) [$cm^{-2} s^{-1} sr^{-1}$]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.legend(loc='lower left',fontsize=11)
	this_ax.grid()
	this_fig.savefig(title+".pdf",edgecolor='none',bbox_inches="tight") #save the figure
	this_fig.savefig(title+".png",edgecolor='none',bbox_inches="tight") #save the figure

main()
