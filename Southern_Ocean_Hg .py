#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:44:57 2018

THIS SCRIPT PLOTS A COMPARISON OF THE 2 OCEAN OPTIONS FOR GEOS-CHEM:
1) Dynamic Ocean vs Offline Ocean

@author: ncp532
"""

# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates           # 
import matplotlib.ticker as ticker

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# Define the datasets

# SIMULATIONS
dataset1a = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Dynamic Ocean Hg Data
dataset1b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.StateMet.2013*.nc4') # Dynamic Ocean Met Data
dataset2a = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Offline Ocean Hg Data
dataset2b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.StateMet.2013*.nc4') # Offline Ocean Met Data
dataset3  = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.2013.nc')  # Default
dataset4  = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.2013.nc') # InvOcean

# OBSERVATIONS
# Retrieve the observational data
#x1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(0), invalid_raise=False) # day number
#y1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

#------------------------------------------------------------------------------
# Define the GEOS-Chem variables
# Note1: To convert from pptv to ng/m3 = pptv*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
# Note2: To convert from mol/mol dry to ng/m3 = mol/mol*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in kg/m3)
    
# Dynamic Ocean
airden1 = dataset1b.variables['Met_AIRDEN'][:,0,7,29] # dry air density (kg/m3)
T1      = dataset1b.variables['Met_T'][:,0,7,29] # Temperature (K)
P1      = dataset1b.variables['Met_DELP'][:,0,7,29] # Delta-pressure across grid box (wet air) (hPa)
P2      = dataset1b.variables['Met_PMID'][:,0,7,29] # Pressure (w/r/t moist air) at level centers (hPa)
P3      = dataset1b.variables['Met_PMIDDRY'][:,0,7,29] # Pressure (w/r/t dry air) at level centers (hPa)
R       = 8.3145 # universal gas constant (J/mol K)
MM_Hg   = 200.59 # Molar Mass Hg (g/mol)
MM_Air  = 28.97 # Molar Mass Dry Air (g/mol)
Avo     = 6.0221*1e23 # Avogadro's number
Hg1     = (dataset1a.variables['SpeciesConc_Hg0'][:,0,7,29]) # dry mixing ratio of Hg0 (mol/mol dry)
Hg1     = Hg1*1e12*(MM_Hg/MM_Air)*airden1

# Offline Ocean
airden2 = dataset2b.variables['Met_AIRDEN'][:,0,7,29] # dry air density (kg/m3)
Hg2     = (dataset2a.variables['SpeciesConc_Hg0'][:,0,7,29]) # dry mixing ratio of Hg0 (mol/mol dry)
Hg2     = Hg2*1e12*(MM_Hg/MM_Air)*airden2

# Default
airden3 = dataset3.variables['BXHGHT_S__AIRNUMDE'][:,0,7,29] # dry air number density (molecules air/cm3)
Hg3     = (dataset3.variables['IJ_AVG_S__Hg0'][:,0,7,29]) # Hg0 (pptv)
Hg3     = Hg3*1e3*MM_Hg/Avo*airden3 # cm3 to m3 is a factor of 1e-6
time    = dataset3.variables['time'][:] # in minutes

# InvOcean
airden4 = dataset4.variables['BXHGHT_S__AIRNUMDE'][:,0,7,29] # dry air number density (molecules air/cm3)
Hg4     = (dataset4.variables['IJ_AVG_S__Hg0'][:,0,7,29]) # Hg0 (pptv)
Hg4     = Hg4*1e3*MM_Hg/Avo*airden4

#------------------------------------------------------------------------------
## Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
datea = []
for t in time:
    hrs=timedelta(hours=int(t))
    datea.append(d0+hrs)

# Default
# Convert datea to monthsa
monthsa = [int(d.strftime('%Y%m')) for d in datea]
monthsa = set(monthsa)
monthsa = np.sort([str(m)+"01" for m in monthsa])
monthsa = [datetime.strptime(m,"%Y%m%d") for m in monthsa]
#------------------------------------------------------------------------------
# Plot the axis for each graph
# Graph 1
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax=plt.subplot(311) # options graph 1 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the observations for the location
#plt.plot(datea,y1, 'k', linewidth=2, label="Observations")

# SIMULATIONS
# Plot the modelling simulations for each location
plt.plot(monthsa, Hg1, ls='--', color='red', linewidth=2, label="Dynamic Ocean") # Dynamic Ocean
plt.plot(monthsa, Hg2, ls='--', color='blue', linewidth=2, label="Offline Ocean") # Offline Ocean
plt.plot(monthsa, Hg3, ls='--', color='green', linewidth=2, label="Default") # Default simulation in GEOS-Chem v11.01
plt.plot(monthsa, Hg4, ls='--', color='orange', linewidth=2, label="InvOcean") # Inverse Ocean simulation in GEOS-Chem v11.01

plt.xlim(datetime(2013,1,1),datetime(2013,12,31)) # set the date limits
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b/%Y') # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')

# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.title("Daily Average Hg$^0$ concentrations (Amsterdam Island, TAAF)", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
