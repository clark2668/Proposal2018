import numpy as np #import numpy
SecPerDay = 86400. #second in a day
KM2toCM2 = 1.e10


import sys, importlib
#importlib.reload(sys.modules['NumberOfEventsPerYearFunction'])
from NumberOfEventsPerYearFunction import NumberOfEventsPerYear



# Example 1 attempted with acceptances based loosely on hybrid plot
Bin_Edges = np.array([3.162e15,1e16,3.16228e16,1e17,3.16228e17,1e18,3.16228e18, 1e19,3.16228e19,1e20,3.16228e20])
AcceptanceExample = KM2toCM2 * np.array([1e-8, 1e-5, 3e-3, 3e-2, 1e-1, 2e-1, 3e-1, 8e-1, 1e0, 2e0, 4e0])
print "Acceptance ",AcceptanceExample
IC_upgoingmuon_limit=np.array([4.970987E-15, 1.263114E-15, 3.209539E-16, 8.155350E-17, 2.072252E-17, 5.265537E-18,1.337958E-18, 3.399716E-19, 8.638585E-20, 2.195041E-20, 5.577539E-21])
IC_upgoingmuon_flux = IC_upgoingmuon_limit/Bin_Edges

nExample=NumberOfEventsPerYear(Bin_Edges, AcceptanceExample, IC_upgoingmuon_flux)
print "nExample is ",nExample
print "Total is ",nExample.sum()

# Example 2 of ARIANNA1296. 5yr acceptance adjust to 1 year
AcceptanceARIANNA1296_5yr = [1734083.197, 14914089.43, 87249881.38, 284353000.5, 1544137142, 5349746464, 13755185405, 29944050061, 59155337434, 110706058196, 191168819784]
nArianna1296_1yr=NumberOfEventsPerYear(Bin_Edges, AcceptanceARIANNA1296_5yr, IC_upgoingmuon_flux)

# Example 3 of ANITA. 17.4 days acceptance adjust to 1 year. had to use new bin size for available acceptance data
AccAnita_17p4days=[3.84E+06, 1.61E+08, 3.11E+09, 2.44E+10, 1.38E+11, 4.43E+11, 1.06E+12]
Bin_Edges_Anita=np.array([1e18,3.16228e18, 1e19,3.16228e19,1e20,3.16228e20, 1e21])
IC_upgoingmuon_limit_ANITA=np.array([5.265537E-18,1.337958E-18, 3.399716E-19, 8.638585E-20, 2.195041E-20, 5.577539E-21, 1.417237E-21])
IC_upgoingmuon_flux_ANITA = IC_upgoingmuon_limit_ANITA/Bin_Edges_Anita

nANITA_1yr=NumberOfEventsPerYear(Bin_Edges_Anita, AccAnita_17p4days, IC_upgoingmuon_flux_ANITA)