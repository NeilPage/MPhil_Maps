#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:23:17 2018

@author: ncp532
"""

# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# Define the datasets

# SIMULATIONS
# Retrieve the model simulation data

#dataset1a = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Dynamic Ocean Hg Data
#dataset1b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.StateMet.2013*.nc4') # Dynamic Ocean Met Data
#dataset2a = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Offline Ocean Hg Data
#dataset2b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.StateMet.2013*.nc4') # Offline Ocean Met Data
#dataset3  = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.2013.nc')  # Default
#dataset4  = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.2013.nc') # InvOcean

# OBSERVATIONS
# Retrieve the observational data

#SIPEX = np.genfromtxt('/Users/ncp532/Documents/Data/SIPEX_II_BrO/all_BrO_vmr_prof_20120917.dat')

SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_BrO/all_BrO_vmr_prof_20121113.dat', sep='\s+', index_col=0) # BrO data for SIPEXII (2012)
#V1 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # BrO data for CAMPCANN V1 (2017/18)
#V2 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # BrO data for CAMPCANN V2 (2017/18)
#V3 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # BrO data for CAMPCANN V3 (2017/18)
#V4 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # BrO data for CAMPCANN V4 (2017/18)

#------------------------------------------------------------------------------
# SET THE DATE AND TIME
d0 = datetime(2012,9,17,0,10,0)
end = datetime(2012,9,17,23,39,59)
step = timedelta(minutes=20)

date = []

while d0 < end:
    date.append(d0.strftime('%Y-%m-%d %H:%M:%S'))
    d0 += step

#CONVERT TO DATETIME FROM STRING
date1=[]
for i in range(len(date)):
    date1.append(datetime.strptime(date[i],'%Y-%m-%d %H:%M:%S'))

date2 = mdates.date2num(date1)

#------------------------------------------------------------------------------
# SET UP THE VALUES TO PLOT

y = SIPEXII.index # set the values for the y-axis
x = np.array(SIPEXII.dtypes.index) # set the values for the x-axis

z = SIPEXII.copy() # identify the matrix containing the z-values (BrO in ppMv)
z[z==-9999]=np.nan # set the erroneous values as NaN 
z = z.loc[:]*1e3 # change from ppMv to ppBv

# when you plot colormaps it changes Nan values to 0.
# Need to mask the array so that Nan values are plotted as grey and not a concentration of 0.
mz=np.ma.masked_where(np.isnan(z),z) 

#------------------------------------------------------------------------------
# PLOT A COLORMAP OF BrO CONCENTRATIONS

fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax=plt.subplot(111)

plt.pcolormesh(date2, y, mz, cmap='RdBu')
plt.title('BrO Concentration (15/10/2012)')
plt.ylabel('Altitude (km)', fontsize=10)
plt.xlabel('Time (UT)', fontsize=10)
# set the limits of the plot to the limits of the data
#plt.axis([x.min(), x.max(), y.min(), y.max()])
clb = plt.colorbar()
clb.set_label(r"BrO VCD (ppbv)", fontsize =10)
# Setup the DateFormatter for the x axis
date_format = mdates.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(date_format)
ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))

# Rotates the labels to fit
fig.autofmt_xdate()

plt.show()