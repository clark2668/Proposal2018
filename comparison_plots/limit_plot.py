# -*- coding: utf-8 -*-
#from __future__ import unicode_literals
import numpy as np #import numpy
import matplotlib.pyplot as plt #import matplotlib
#import matplotlib
from matplotlib.pyplot import rcParams
rcParams['mathtext.default'] = 'regular'
import math


#plt.rc('text', usetex=True)
#plt.rc('font', family='sans-serif')

#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True

def main():

	ara2_data = np.genfromtxt("limits/ara2_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	ara2_limit=ara2_data['limit']
	ara2_aperture = 2.3/ara2_limit/(445.*86400.)/1e10/2

	testbed_data = np.genfromtxt("limits/testbed_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	testbed_limit=testbed_data['limit']
	testbed_aperture = 2.3/testbed_limit/(415.*86400.)/1e10

	hra3_data = np.genfromtxt("limits/hra3_limit.csv",delimiter=',',skip_header=1,names=['energy','limit'])
	hra3_limit=hra3_data['limit']/(hra3_data['energy']/1e9) #need to divide out the GeV
	hra3_aperture = 2.3/hra3_limit/(170.*86400.)/1e10/3
	#this is 170 days of data

	fig_efe = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_efe = fig_efe.add_subplot(1,1,1) #make a subplot
	#ax_efe.plot(testbed_data['energy'],testbed_data['limit'],'-v', linewidth=2.0,color='magenta',label=r'1 $\times$ (415 Days of 1 ARA 30m Testbed)')
	ax_efe.plot(testbed_data['energy'],testbed_data['limit']*(415./365.),'-v', linewidth=2.0,color='red',label=r'365 Days of 1 ARA 30m Testbed (Analysis Level)')
	#ax_efe.plot(ara2_data['energy'],ara2_data['limit'],'-s', linewidth=2.0,color='blue',label=r'2 $\times$ (222 Days of 1 ARA 200m Station)')
	ax_efe.plot(hra3_data['energy'],hra3_limit*(3.*170./365.),'-o', linewidth=2.0,color='green',label=r'365 Days of 1 ARIANNA Station (Analysis Level)')
	ax_efe.plot(ara2_data['energy'],ara2_data['limit']*(2.*222./365.),'-s', linewidth=2.0,color='blue',label=r'365 Days of 1 ARA 200m Station (Analysis Level)')
	save_limit(fig_efe,ax_efe)

	fig_aff = plt.figure(figsize=(11,8.5)) #make a figure object
	ax_aff = fig_aff.add_subplot(1,1,1) #make a subplot
	ax_aff.plot(ara2_data['energy'],ara2_aperture,'-s', linewidth=2.0,color='blue',label="ARA Deep Station (Analysis Level)")
	ax_aff.plot(hra3_data['energy'],hra3_aperture,'-o', linewidth=2.0,color='green',label="ARIANNA Station (Analysis Level)")
	ax_aff.plot(testbed_data['energy'],testbed_aperture,'-v', linewidth=2.0,color='red',label="ARA Testbed Station (Analysis Level)")
	
	save_aff(fig_aff,ax_aff)


def ara2(axis,color):
	#this is the limit from Thomas' PRD for the analysis level, 2 stations after 7.5 month of livetime
	print "unused function"


def save_aff(this_fig, this_ax):
	sizer=20
	xlow = 1.00e14 #the lower x limit
	xup = 1e21 #the uppper x limit
	this_ax.set_xlabel('Energy [eV]',size=sizer) #give it a title
	this_ax.set_ylabel('A$_{eff}$ $\Omega$ [km$^2$sr]',size=sizer)
	this_ax.set_yscale('log')
	this_ax.set_xscale('log')
	this_ax.tick_params(labelsize=sizer)
	this_ax.set_xlim([xlow,xup]) #set the x limits of the plot
	#this_ax.set_ylim([ylow,yup]) #set the y limits of the plot
	this_ax.grid()
	this_ax.legend(loc='upper left')
	this_fig.savefig('aff.pdf',edgecolor='none',bbox_inches="tight") #save the figure



def save_limit(this_fig, this_ax):
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
	this_ax.legend(loc='lower left')
	this_ax.grid()
	this_fig.savefig('limit.pdf',edgecolor='none',bbox_inches="tight") #save the figure


main()
