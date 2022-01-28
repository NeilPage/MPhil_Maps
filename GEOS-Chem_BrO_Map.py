#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 10:55:05 2020

@author: ncp532
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:22:07 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
#from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
from netCDF4 import Dataset, MFDataset
import xarray as xr

# DATA HANDLING PACKAGES
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm                   # imports the colormap function from matplotlib
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates           # 
import matplotlib.ticker as ticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# GEOS-Chem Simulations
DS1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/trac_avg.geosfp_2x25_Hg.2013*.nc') # Dynamic Ocean
DS2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_LPOLARBR/trac_avg.geosfp_2x25_Hg.2013*.nc') # LPOLARBR

# Observations
V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V1_17_Data_SI.csv',header=0,encoding = 'unicode_escape')
#SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/SIPEXII_Data.csv',header=0,encoding = 'unicode_escape')

#------------------------------------------------------------------------------
# FIX UP THE BrO OBSERVATIONS

# Set the date
V1_17['DateTime'] = pd.to_datetime(V1_17['DateTime'])

# V1_17 Davis (07:00 to 18:00)
start_time = '07:00:00'
end_time = '18:00:00'
Midday = (V1_17['Time'] >= start_time) & (V1_17['Time'] < end_time)
V1_17_MM = V1_17[Midday]

# Filter dataframe for when filter is less than 60%
V1_17F = (V1_17_MM['Filter'] < 0.6)
V1_17T = V1_17_MM[V1_17F]

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

R       = 8.3145 # universal gas constant (J/mol K)
MM_Air  = 28.97 # Molar Mass Dry Air (g/mol)
Avo     = 6.0221 * 1e23 # Avogadro's number
MM_BrO   = 95.904 # Molar Mass BrO (g/mol)

#----------------------
# OBSERVATIONS
#----------------------
# BrO VMR
BrO_Obs = np.array(V1_17T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to ppbv

#----------------------
# SIMULATIONS (DS1)
#----------------------
# Air Density
airden1 = DS1.variables[u'BXHGHT_S__AIRNUMDE'][:] # dry air number density (molec air/cm3) 

# BrO (Convert from molec/cm3 to ppbv)
#BrO_Sim1     = DS1.variables[u'PL_HG2_S__BrO'][:] # BrO concentration (molec/cm3) 
#BrO_Sim1     = BrO_Sim1 * (1e-3) * MM_BrO / (Avo) * airden1
BrO_Sim1      = DS1.variables['PL_HG2_S__BrO_pol'] # Polar BrO concentration (pptv)

# Time
time1         = DS1.variables['time'][:]# in minutes

#----------------------
# SIMULATIONS (DS2)
#----------------------
# Air Density
airden2 = DS2.variables[u'BXHGHT_S__AIRNUMDE'][:] # dry air number density (molec air/cm3) 

# BrO
#BrO_Sim2     = DS2.variables[u'PL_HG2_S__BrO'][:] # BrO concentration (molec/cm3) 
#BrO_Sim2     = BrO_Sim2 * (1e-3) * MM_BrO / (Avo) * airden1
BrO_Sim2      = DS2.variables['PL_HG2_S__BrO_pol'] # Polar BrO concentration (pptv)

# Time
time2         = DS2.variables['time'][:]# in minutes

#------------------------------------------------------------------------------
# SET THE DATE AND TIME

#----------------------
# OBSERVATIONS
#----------------------
dat = np.array(V1_17T['Date'])
tim = np.array(V1_17T['Time'])
dattim = dat+' '+tim

#CONVERT TO DATETIME FROM STRING
date1=[]
for i in range(len(dattim)):
    date1.append(datetime.strptime(dattim[i],'%d/%m/%Y %H:%M:%S'))

#----------------------
# SIMULATIONS
#----------------------
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)

# DS1
date2 = []
for t in time1:
    hrs=timedelta(hours=int(t))
    date2.append(d0+hrs)

# DS2
date3 = []
for t in time2:
    hrs=timedelta(hours=int(t))
    date3.append(d0+hrs)

#------------------------------------------------------------------------------    
# CALCULATE THE DAILY MEAN
def daily(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('D').mean()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 
# OBSERVATIONS  
BrO_Obs_DM, date1_DM = daily(BrO_Obs,date1) # V1_17

#------------------------------------------------------------------------------
# CALCULATE THE DAILY MEDIAN
def daily(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('D').median()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 
# OBSERVATIONS  
BrO_Obs_DMed, date1_DMed = daily(BrO_Obs,date1) # V1_17

#------------------------------------------------------------------------------
# PLOT THE GRAPH
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

#------------------------------
# GRAPH 1
#------------------------------
ax=plt.subplot(111) # options graph 1 (vertical no, horizontal no, graph no)

# PLOT LINES FOR SIM & OBS
ax.plot(date1_DM, BrO_Obs_DM, marker='+', c='blue', markersize = 3.0, linestyle='None', label="Observations (Daily Mean)") # Daily Mean
ax.plot(date1_DMed, BrO_Obs_DMed, marker='+', c='green', markersize = 3.0, linestyle='None', label="Observations (Daily Median)") # Daily Median
ax.plot(date2, BrO_Sim1[:,0,65,101], marker='+', c='red', markersize = 3.0, linestyle='None', label="DS1") # DS1
ax.plot(date3, BrO_Sim2[:,0,65,101], marker='+', c='black', markersize = 3.0, linestyle='None', label="DS2") #DS2

# FORMAT THE X-AXIS
#plt.xlim(datetime(2017,11,13),datetime(2017,11,23)) # at Davis station
plt.xticks(rotation=40)
date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')

# FORMAT THE Y-AXIS
ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_ticks_position('both')

# PLOT THE AXIS LABELS, LEGEND AND TITLE
ax.set_ylabel('BrO VMR\n(pptv)', fontsize=15)
ax.set_xlabel('Date', fontsize=15)
plt.title('Observed and simulated BrO concentrations for V1 (CAMMPCAN 2017-18)', fontsize=25, y=1.05)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)