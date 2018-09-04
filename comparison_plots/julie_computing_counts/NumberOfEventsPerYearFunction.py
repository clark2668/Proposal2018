# Function to calculate the number of neutrino events seen by a detector in one year
# Julie Rolla
# 8/7/2018

# Function based on the equation for number of events below:
# N = int F(E) dE * int dt * int dA_eff * int dOmega

# NumberOfEventsPerYear(BinEdges, DetectorAcceptance, Flux)
# Inputs:
#       BinEdges = list with endpoints of the energy bins 
#            for example:
#                      Bin_Edges = np.array([3.162e15,1e16,3.16228e16,1e17,3.16228e17,1e18,3.16228e18, 1e19,3.16228e19,1e20,3.16228e20])
#       DetectorAcceptance = array of acceptances ([A*Omega]_effective) of the detector at each energy bin edge in cm^2 sr (use factor of 1.e10 to convert from km^2 to cm^2)
#       Flux = the estimated neutrino flux values used (F(E), not EF(E)) in cm^-2 sr^-1 s^-1
#           for example:   
#                     IC_upgoingmuon_limit=np.array([4.970987E-15, 1.263114E-15, 3.209539E-16, 8.155350E-17, 2.072252E-17, 5.265537E-18,1.337958E-18, 3.399716E-19, 8.638585E-20, 2.195041E-20, 5.577539E-21])
#                     IC_upgoingmuon_flux = IC_upgoingmuon_limit/Bin_Edges
#Output:
#       N = number of event the detector will see per year assuming the given flux
#
# To run in anther program, put this file in working directory, and type from filename import functionname. ie:
# from NumberOfEventsPerYearFunction import NumberOfEventsPerYear
# 
# If you need to edit this function, use the following three lines instead:
# import sys, importlib
# importlib.reload(sys.modules['NumberOfEventsPerYearFunction'])
# from NumberOfEventsPerYearFunction import NumberOfEventsPerYear



import numpy as np #import numpy
SecPerDay = 86400. #second in a day
SecPerYear = 365.*SecPerDay #second in a year

def NumberOfEventsPerYear(BinEdges, DetectorAcceptance, Flux): # function name and inputs
    # Create empty arrays of the integral of the acceptance, integral of flux, and number of events
    Int_Acc = []
    Int_Flux = []
    N=[]
    
    NumberOfBins=len(BinEdges)-1 # calculates number of energy bins
    
    # for loop that uses the trapazoid rule to integrate the acceptance and flux and calculate N for each energy bin
    for x in range(0,NumberOfBins):
        # First performs trapazoid rule on acceptance data and store it
        acc = np.trapz([DetectorAcceptance[x],DetectorAcceptance[x+1]], x=[BinEdges[x],BinEdges[x+1]])
        Int_Acc.append(acc)
        
        # Next performs trapazoid rule on flux data and stores it
        fluxint = np.trapz([Flux[x],Flux[x+1]], x=[BinEdges[x],BinEdges[x+1]]) 
        Int_Flux.append(fluxint)
    
        # Finally, the loop calculates the number of events in the bin
        N.append(Int_Flux[x]*SecPerYear*Int_Acc[x]/BinEdges[x])
        # **** I am not 100% sure why I need to divide by BinEdges. I think the two integrals end up counting the energy bins twice...???
   
    N = np.array(N) # Turns list into array
    #print("N =", N) #uncomment to print when testing
    return N
