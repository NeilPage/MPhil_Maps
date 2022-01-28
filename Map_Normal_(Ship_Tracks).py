#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:19:13 2017

@author: ncp532
"""
# FILE SYSTEM PACKAGES
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# DRAWING PACKAGES
from mpl_toolkits.basemap import Basemap	# function to plot maps
import matplotlib.pyplot as plt			# imports pyplot (provides a MATLAB-like plotting function)
import matplotlib as mpl				    # import matplotlib package as shorter nickname
import matplotlib.colors as colors			# for visual representation of matplotlib colormaps
from matplotlib import cm                   # imports the colormap function from matplotlib

# DATE AND TIME HANDLING PACKAGES
from datetime import datetime,timedelta		# functions to handle date and time

# DATA HANDLING PACKAGES
import math						         # package to do mathematical calculations (such as log, exponential)
import numpy as np					    # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.
import pandas as pd

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# SIMULATIONS
#dataset1a = MFDataset('/Users/macbook/Documents/UOW/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Dynamic Ocean Hg Data
#dataset1b = MFDataset('/Users/macbook/Documents/UOW/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.StateMet.2013*.nc4') # Dynamic Ocean Met Data
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_LPOLARBR/trac_avg.geosfp_2x25_Hg.2013*.nc') # Offline Ocean Hg Data
#dataset2b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.StateMet.2013*.nc4') # Offline Ocean Met Data

# OBSERVATIONS

#SIPEXII
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for SIPEXII (2012)
SIPEXII_Long = np.array(SIPEXII['LONGITUDE']) # Longitude for SIPEXII
SIPEXII_Lat = np.array(SIPEXII['LATITUDE']) # Latitude for SIPEXII
SIPEXII_Hg0 = np.array(SIPEXII['ng/m3']) # Hg0 for SIPEXII

#PCAN
RVInv = pd.read_csv('/Users/ncp532/Documents/Data/RVInvestigator_Hg_2017/RVInvestigator_Hg0_Lat_Long.csv') # Hg0 data for RVInv (2017)
RVInv_Long = np.array(RVInv['LONGITUDE']) # Longitude for RVInv
RVInv_Lat = np.array(RVInv['LATITUDE']) # Latitude for RVInv
RVInv_Hg0 = np.array(RVInv['ng/m3']) # Hg0 for RVInv

#RV Xuelong
RVXue = pd.read_csv('/Users/ncp532/Documents/Data/RVXuelong_Hg_2014/RVXuelong_Hg0_Lat_Long.csv') # Hg0 data for RVXue (2014/15)
RVXue_Long = np.array(RVXue['LONGITUDE']) # Longitude for RVXue
RVXue_Lat = np.array(RVXue['LATITUDE']) # Latitude for RVXue
RVXue_Hg0 = np.array(RVXue['ng/m3']) # Hg0 for RVXue

#CHINARE
CHINARE = pd.read_csv('/Users/ncp532/Documents/Data/CHINARE_2012_13/CHINARE_Hg0_Lat_Long.csv') # Hg0 data for CHINARE
CHINARE_Long = np.array(CHINARE['LONGITUDE']) # Longitude for CHINARE
CHINARE_Lat = np.array(CHINARE['LATITUDE']) # Latitude for CHINARE
CHINARE_Hg0 = np.array(CHINARE['ng/m3']) # Hg0 for CHINARE

#OSO V1 (2010)
OSO_V1 = pd.read_csv('/Users/ncp532/Documents/Data/OSO_2010_11/OSO_V1_Hg0_Lat_Long.csv') # Hg0 data for OSO
OSO_V1_Long = np.array(OSO_V1['LONGITUDE']) # Longitude for OSO
OSO_V1_Lat = np.array(OSO_V1['LATITUDE']) # Latitude for OSO
#OSO_V1_Hg0 = np.array(OSO_V1['ng/m3']) # Hg0 for OSO

#OSO V2 (2011)
OSO_V2 = pd.read_csv('/Users/ncp532/Documents/Data/OSO_2010_11/OSO_V2_Hg0_Lat_Long.csv') # Hg0 data for OSO
OSO_V2_Long = np.array(OSO_V2['LONGITUDE']) # Longitude for OSO
OSO_V2_Lat = np.array(OSO_V2['LATITUDE']) # Latitude for OSO
#OSO_V2_Hg0 = np.array(OSO_V2['ng/m3']) # Hg0 for OSO

#OSO V3 (2011)
OSO_V3 = pd.read_csv('/Users/ncp532/Documents/Data/OSO_2010_11/OSO_V3_Hg0_Lat_Long.csv') # Hg0 data for OSO
OSO_V3_Long = np.array(OSO_V3['LONGITUDE']) # Longitude for OSO
OSO_V3_Lat = np.array(OSO_V3['LATITUDE']) # Latitude for OSO
#OSO_V2_Hg0 = np.array(OSO_V3['ng/m3']) # Hg0 for OSO

#ANTXXIX_V6
ANTXXIX_V6 = pd.read_csv('/Users/ncp532/Documents/Data/ANTXXIX/ANTXXIX_V6.txt',delimiter="\t") # Hg0 data for ANTXXIX_V6
ANTXXIX_V6_Long = np.array(ANTXXIX_V6['LONGITUDE']) # Longitude for ANTXXIX_V6
ANTXXIX_V6_Lat = np.array(ANTXXIX_V6['LATITUDE']) # Latitude for ANTXXIX_V6
#ANTXXIX_V6_Hg0 = np.array(ANTXXIX_V6['ng/m3']) # Hg0 for ANTXXIX_V6

#ANTXXIX_V7
ANTXXIX_V7 = pd.read_csv('/Users/ncp532/Documents/Data/ANTXXIX/ANTXXIX_V7.txt',delimiter="\t") # Hg0 data for ANTXXIX_V7
ANTXXIX_V7_Long = np.array(ANTXXIX_V7['LONGITUDE']) # Longitude for ANTXXIX_V7
ANTXXIX_V7_Lat = np.array(ANTXXIX_V7['LATITUDE']) # Latitude for ANTXXIX_V7
#ANTXXIX_V7_Hg0 = np.array(ANTXXIX_V7['ng/m3']) # Hg0 for ANTXXIX_V7

#CAMMPCAN (2017/18)
#V1
V1_17_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V1_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V1 (2017/18)
V1_17_18_Long = np.array(V1_17_18['LONGITUDE']) # Longitude for V1(17_18)
V1_17_18_Lat = np.array(V1_17_18['LATITUDE']) # Latitude for V1(17_18)
V1_17_18_Hg0 = np.array(V1_17_18['ng/m3']) # Hg0 for V1(17_18)

#V2
V2_17_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V2_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V2 (2017/18)
V2_17_18_Long = np.array(V2_17_18['LONGITUDE']) # Longitude for V2(17_18)
V2_17_18_Lat = np.array(V2_17_18['LATITUDE']) # Latitude for V2(17_18)
V2_17_18_Hg0 = np.array(V2_17_18['ng/m3']) # Hg0 for V2(17_18)

#V3
V3_17_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V3_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V3 (2017/18)
V3_17_18_Long = np.array(V3_17_18['LONGITUDE']) # Longitude for V3(17_18)
V3_17_18_Lat = np.array(V3_17_18['LATITUDE']) # Latitude for V3(17_18)
V3_17_18_Hg0 = np.array(V3_17_18['ng/m3']) # Hg0 for V3(17_18)

#V4
V4_17_18 = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2017_18/CAMMPCAN_V4_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V4 (2017/18)
V4_17_18_Long = np.array(V4_17_18['LONGITUDE']) # Longitude for V4(17_18)
V4_17_18_Lat = np.array(V4_17_18['LATITUDE']) # Latitude for V4(17_18)
V4_17_18_Hg0 = np.array(V4_17_18['ng/m3']) # Hg0 for V4(17_18)

#CAMMPCAN (2018/19)
#V1
#V1(18_19) = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1(18/19)_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V1 (2017/18)
#V1(18_19)_Long = np.array(V4(17_18)['LONGITUDE']) # Longitude for V1(18_19)
#V1(18_19)_Lat = np.array(V4(17_18)['LATITUDE']) # Latitude for V1(18_19)
#V1(18_19)_Hg0 = np.array(V4(17_18)['ng/m3']) # Hg0 for V1(18_19)

#V2
#V2(18_19) = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1(18/19)_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V2 (2017/18)
#V2(18_19)_Long = np.array(V4(17_18)['LONGITUDE']) # Longitude for V2(18_19)
#V2(18_19)_Lat = np.array(V4(17_18)['LATITUDE']) # Latitude for V2(18_19)
#V2(18_19)_Hg0 = np.array(V4(17_18)['ng/m3']) # Hg0 for V2(18_19)

#V3
#V3(18_19) = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1(18/19)_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V3 (2017/18)
#V3(18_19)_Long = np.array(V4(17_18)['LONGITUDE']) # Longitude for V3(18_19)
#V3(18_19)_Lat = np.array(V4(17_18)['LATITUDE']) # Latitude for V3(18_19)
#V3(18_19)_Hg0 = np.array(V4(17_18)['ng/m3']) # Hg0 for V3(18_19)

#V4
#V4(18_19) = pd.read_csv('/Users/ncp532/Documents/Data/CAMMPCAN_2018_19/V1(18/19)_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V4 (2017/18)
#V4(18_19)_Long = np.array(V4(17_18)['LONGITUDE']) # Longitude for V4(18_19)
#V4(18_19)_Lat = np.array(V4(17_18)['LATITUDE']) # Latitude for V4(18_19)
#V4(18_19)_Hg0 = np.array(V4(17_18)['ng/m3']) # Hg0 for V4(18_19)

#------------------------------------------------------------------------------
# GET THE VARIABLES YOU WOULD LIKE FROM THE GEOS-CHEM SIMULATION
# you can display the variable names by typing "dataset1" into the python console
# you can display more information about individual variables by typing "dataset1.variables['INSERT THE VARIABLES NAME']" into the python console
# e.g. for Hg0 type "dataset1.variables['IJ_AVG_S__Hg0']" into the python console

# Dynamic Ocean
#airden1 = dataset1b.variables['Met_AIRDEN'][:] # dry air density (kg/m3)
#T1      = dataset1b.variables['Met_T'][:] # Temperature (K)
#P1      = dataset1b.variables['Met_DELP'][:] # Delta-pressure across grid box (wet air) (hPa)
#P2      = dataset1b.variables['Met_PMID'][:] # Pressure (w/r/t moist air) at level centers (hPa)
#P3      = dataset1b.variables['Met_PMIDDRY'][:] # Pressure (w/r/t dry air) at level centers (hPa)
#R       = 8.3145 # universal gas constant (J/mol K)
#MM_Hg   = 200.59 # Molar Mass Hg (g/mol)
#MM_Air  = 28.97 # Molar Mass Dry Air (g/mol)
#Avo     = 6.0221*1e23 # Avogadro's number
#Hg1     = (dataset1a.variables['SpeciesConc_Hg0'][:]) # dry mixing ratio of Hg0 (mol/mol dry)
#Hg1     = Hg1*1e12*(MM_Hg/MM_Air)*airden1
#
#Hg_Av   = (Hg1[8,0,:,:] + Hg1[9,0,:,:]) / 2 # Mean simulated Hg0 concentration for Sept and Oct

# Offline Ocean
#airden2 = dataset2b.variables['Met_AIRDEN'][:] # dry air density (kg/m3)
#Hg2     = (dataset2a.variables['SpeciesConc_Hg0'][:]) # dry mixing ratio of Hg0 (mol/mol dry)
#Hg2     = Hg2*1e12*(MM_Hg/MM_Air)*airden2

# LPOLARBR
airden1 = dataset1.variables[u'BXHGHT_S__AIRNUMDE'][:] # dry air number density (molec air/cm3) 
R       = 8.3145 # universal gas constant (J/mol K)
MM_Air  = 28.97 # Molar Mass Dry Air (g/mol)
Avo     = 6.0221*1e23 # Avogadro's number
# Hg0
MM_Hg   = 200.59 # Molar Mass Hg (g/mol)
Hg0     = dataset1.variables[u'IJ_AVG_S__Hg0'][:] # Hg0 (pptv) 
Hg0     = Hg0*(1e-3)*MM_Hg/(Avo)*airden1*1e6 
Hg0_Av   = (Hg0[8,0,:,:] + Hg0[9,0,:,:]) / 2 # Mean simulated Hg0 concentration for Sept and Oct

# BrO
MM_BrO   = 95.904 # Molar Mass BrO (g/mol)
BrO     = dataset1.variables[u'PL_HG2_S__BrO'][:] # BrO concentration (molec/cm3) 
BrO     = BrO*(1e-3)*MM_BrO/(Avo)*airden1
BrO_Av   = (BrO[8,0,:,:] + BrO[9,0,:,:]) / 2 # Mean simulated BrO concentration for Sept and Oct

# Br
MM_Br   = 79.904 # Molar Mass Br (g/mol)
Br      = dataset1.variables[u'PL_HG2_S__Br'][:] # BrO concentration (molec/cm3) 
Br      = Br*(1e-3)*MM_Br/(Avo)*airden1
Br_Av   = (Br[8,0,:,:] + Br[9,0,:,:]) / 2 # Mean simulated Br concentration for Sept and Oct

#------------------------------------------------------------------------------
# FIX THE LATITUDE AND LONGITUDE FOR PLOTTING (based on a size of 2x2.5 for each gridbox) 
# NOTE: GEOS-Chem defines the latitude and longitude based on the bottom-left corner of each grid-box
#       we have to correct the latitude and longitude to represent the middle of each grid-box
#       if we dont do this the map will be slightly off centre when plotted

# Values for the longitudinal edge of each gridbox
lons_e = np.array([ -181.25, -178.75, -176.25, -173.75, -171.25, -168.75,
                    -166.25, -163.75, -161.25, -158.75, -156.25, -153.75,
                    -151.25, -148.75, -146.25, -143.75, -141.25, -138.75,
                    -136.25, -133.75, -131.25, -128.75, -126.25, -123.75,
                    -121.25, -118.75, -116.25, -113.75, -111.25, -108.75,
                    -106.25, -103.75, -101.25,  -98.75,  -96.25,  -93.75,
                    -91.25,  -88.75,  -86.25,  -83.75, -81.25,  -78.75,
                    -76.25,  -73.75,  -71.25,  -68.75,  -66.25,  -63.75,
                    -61.25,  -58.75,  -56.25,  -53.75,  -51.25,  -48.75,
                    -46.25,  -43.75,  -41.25,  -38.75,  -36.25,  -33.75,
                    -31.25,  -28.75,  -26.25,  -23.75, -21.25,  -18.75,
                    -16.25,  -13.75,  -11.25,   -8.75,   -6.25,   -3.75,
                    -1.25,    1.25,    3.75,    6.25,    8.75,   11.25,
                    13.75,   16.25,   18.75,   21.25,   23.75,   26.25,
                      28.75,   31.25,   33.75,   36.25,  38.75,   41.25,
                      43.75,   46.25,   48.75,   51.25,   53.75,   56.25,
                      58.75,   61.25,   63.75,   66.25,   68.75,   71.25,
                      73.75,   76.25,   78.75,   81.25,   83.75,   86.25,
                      88.75,   91.25,   93.75,   96.25,   98.75,  101.25,
                      103.75,  106.25,  108.75,  111.25,  113.75,  116.25,
                      118.75,  121.25,  123.75,  126.25,  128.75,  131.25,
                      133.75,  136.25,  138.75,  141.25,  143.75,  146.25,
                      148.75,  151.25,  153.75,  156.25,  158.75,  161.25,
                      163.75,  166.25,  168.75,  171.25,  173.75,  176.25,
                      178.75,])
# Values for the longitudinal middle of each gridbox
lons_m= np.array([-180.0, -177.5, -175.0, -172.5, -170.0, -167.5, -165.0, -162.5,
                       -160.0, -157.5, -155.0, -152.5, -150.0, -147.5, -145.0, -142.5,
                       -140.0, -137.5, -135.0, -132.5, -130.0, -127.5, -125.0, -122.5,
                       -120.0, -117.5, -115.0, -112.5, -110.0, -107.5, -105.0, -102.5,
                       -100.0,  -97.5,  -95.0,  -92.5,  -90.0,  -87.5,  -85.0,  -82.5,
                       -80.0,  -77.5,  -75.0,  -72.5,  -70.0,  -67.5,  -65.0,  -62.5,
                       -60.0,  -57.5,  -55.0,  -52.5,  -50.0,  -47.5,  -45.0,  -42.5,
                       -40.0,  -37.5,  -35.0,  -32.5,  -30.0,  -27.5,  -25.0,  -22.5,
                       -20.0,  -17.5,  -15.0,  -12.5,  -10.0,   -7.5,   -5.0,   -2.5,
                       0.0,    2.5,    5.0,    7.5,   10.0,   12.5,   15.0,   17.5,
                       20.0,   22.5,   25.0,   27.5,   30.0,   32.5,   35.0,   37.5,
                       40.0,   42.5,   45.0,   47.5,   50.0,   52.5,   55.0,   57.5,
                       60.0,   62.5,   65.0,   67.5,   70.0,   72.5,   75.0,   77.5,
                       80.0,   82.5,   85.0,   87.5,   90.0,   92.5,   95.0,   97.5,
                       100.0,  102.5,  105.0,  107.5,  110.0,  112.5,  115.0,  117.5,
                       120.0,  122.5,  125.0,  127.5,  130.0,  132.5,  135.0,  137.5,
                       140.0,  142.5,  145.0,  147.5,  150.0,  152.5,  155.0,  157.5,
                       160.0,  162.5,  165.0,  167.5,  170.0,  172.5,  175.0,  177.5,
                       ])
# Values for the latitudinal edge of each gridbox
lats_e=np.array([  -90.,  -89.,  -87.,  -85.,  -83.,  -81.,  -79.,  -77.,
                 -75.,  -73.,  -71.,  -69.,  -67.,  -65.,  -63.,  -61.,
                 -59.,  -57.,  -55.,  -53.,  -51.,  -49.,  -47.,  -45.,
                 -43.,  -41.,  -39.,  -37.,  -35.,  -33.,  -31.,  -29.,
                 -27.,  -25.,  -23.,  -21.,  -19.,  -17.,  -15.,  -13.,
                 -11.,   -9.,   -7.,   -5.,   -3.,   -1.,    1.,    3.,
                 5.,    7.,    9.,   11.,   13.,   15.,   17.,   19.,
                 21.,   23.,   25.,   27.,   29.,   31.,   33.,   35.,
                 37.,   39.,   41.,   43.,   45.,   47.,   49.,   51.,
                 53.,   55.,   57.,   59.,   61.,   63.,   65.,   67.,
                 69.,   71.,   73.,   75.,   77.,   79.,   81.,   83.,
                 85.,   87.,   89.,   90., ])
# Values for the latitudinal edge of each gridbox
lats_m=np.array([  -89.5, -88.,  -86.,  -84.,  -82.,  -80.,  -78.,  -76.,
                -74.,  -72.,  -70.,  -68.,  -66.,  -64.,  -62.,  -60.,
                -58.,  -56.,  -54.,  -52.,  -50.,  -48.,  -46.,  -44.,
                -42.,  -40.,  -38.,  -36.,  -34.,  -32.,  -30.,  -28.,
                -26.,  -24.,  -22.,  -20.,  -18.,  -16.,  -14.,  -12.,
                -10.,   -8.,   -6.,   -4.,   -2.,    0.,    2.,    4.,
                6.,    8.,   10.,   12.,   14.,   16.,   18.,   20.,
                22.,   24.,   26.,   28.,   30.,   32.,   34.,   36.,
                38.,   40.,   42.,   44.,   46.,   48.,   50.,   52.,
                54.,   56.,   58.,   60.,   62.,   64.,   66.,   68.,
                70.,   72.,   74.,   76.,   78.,   80.,   82.,   84.,
                86.,   88.,   89.5, ])

#------------------------------------------------------------------------------
# DEFFINE THE FUNCTION FOR PLOTTING YOUR MAP
def plotkg(x, bar, tit):  # (x = the value to be plotted, bar = colorbar (Y/N), tit = title for map)
    lons, lats = np.meshgrid(lons_e,lats_e)
    
    # NORMAL MAP PLOT
    m=Basemap(llcrnrlat=-90,  urcrnrlat=-20,          # set the latitude limits
              llcrnrlon=-180, urcrnrlon=180,         # set the longitude limits
              resolution='c',projection='gall')      # choose the resolution and projection type
    
    # SOUTH POLAR STEROGRAPHIC PROJECTION          
   #m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='c')          
              
    mlons,mlats = m(lons, lats) # convert lons/lat from middle of grid-box to to corner of grid-box
    #cs=m.pcolormesh(mlons, mlats, x, latlon=False, vmin=0.4, vmax=1.6, cmap=cm.plasma) # define the colorbar attributes (colors used, add "_r" to flip the color scale)    
    m.drawcoastlines() # draw coastlines on the map
    m.fillcontinents(color = 'gray', zorder = 0)
    
#    # FORMAT THE COLORBAR
#    if bar=='Y': # Y if you want a colorbar, N if you don't
#        cb=m.colorbar(cs,"bottom",size="5%", pad="20%") # Position the colorbar next to your map and define its size
#        cb.set_label(r"$ng$ $m$$^-$$^3$", fontsize =20) # add text for the colourbar label and define fontsize
#        cb.ax.tick_params(labelsize=20) # add numbers to the colorbar
#        cb.ax.minorticks_on()           # add ticks to the colorbar
#    
    # PLOT PARALLELS AND MERIDANS
    parallels = np.arange(-90,10, 10.0)  # set the latitude limits and how often you wish parallels to be plotted
    m.drawparallels(parallels,labels=[True,False,True,False], fontsize=15) # define the fontsize for parallels
    meridians = np.arange(-180,179,30.0) # set the longitude limits and how often you wish meridans to be plotted
    m.drawmeridians(meridians,labels=[True,False,False,True], fontsize=15) # define the fontsize for meridans
    
    # PLOT LOCATION MARKERS FOR MONITORING SITES
    lons1 = [144.6832, 77.3423, 18.4897, 140.0007, 123.1956, 2.3206, -71.2512, 166.4006]    # [Cape Grim, Amsterdam Island, Cape Point, Dumont d'Urville, Concordia, Troll, Bariloche, McMurdo]
    lats1 = [-40.6832, -37.4748, -34.3535, -66.3946, -75.0559, -72.0041, -41.7438, -77.5047]  # [Cape Grim, Amsterdam Island, Cape Point, Dumont d'Urville, Concordia, Troll, Bariloche, McMurdo]
    x1,y1 = m(lons1, lats1)
    
    # OPTION 1: if you want plain location markers use this line
    m.scatter(x1, y1, c='black', edgecolors='black', s=50, zorder = 10) # define the location markers (color, edgecolor, size)
    
    # OPTION2: if you want the location markers to represent observed Hg concentrations use these two lines
    #obs = [0.909, 0.974, 0.808, 1.032, 1.092] # define the observed Hg concentration for each site
    
    # PLOT THE LOCATION MARKERS FOR AAD STATIONS
    lons2 = [62.5227, 77.5803, 110.3136, 158.5609]    # [Mawson, Davis, Casey, Macquarie Island]
    lats2 = [-67.3612, -68.3436, -66.1657, -54.2959]  # [Mawson, Davis, Casey, Macquarie Island]
    x2,y2 = m(lons2, lats2)
    
    # OPTION 1: if you want plain location markers use this line
    m.scatter(x2, y2, c='black', edgecolors='black', s=50) # define the location markers (color, edgecolor, size)
    
    #---------------------------
    # PLOT THE SHIP TRACKS
    # SIPEXII
    lons3 = SIPEXII_Long
    lats3 = SIPEXII_Lat
    x3,y3 = m(lons3, lats3)
    obs = SIPEXII_Hg0
    m.scatter(x3, y3, c='green', s=10, label="SIPEXII") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
        
    # PCAN
    lons4 = RVInv_Long
    lats4 = RVInv_Lat
    x4,y4 = m(lons4, lats4)
    obs = RVInv_Hg0
    m.scatter(x4, y4, c='blue', s=10, label="PCAN") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    
    # CAMMPCAN (2017/18)
    #V1
    lons5 = V1_17_18_Long
    lats5 = V1_17_18_Lat
    x5,y5 = m(lons5, lats5)
    obs = V1_17_18_Hg0
    m.scatter(x5, y5, c='yellow', s=10, label="CAMMPCAN V1") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V2
    lons6 = V2_17_18_Long
    lats6 = V2_17_18_Lat
    x6,y6 = m(lons6, lats6)
    obs = V2_17_18_Hg0
    m.scatter(x6, y6, c='red', s=10, label="CAMMPCAN V2") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V3
    lons7 = V3_17_18_Long
    lats7 = V3_17_18_Lat
    x7,y7 = m(lons7, lats7)
    obs = V3_17_18_Hg0
    m.scatter(x7, y7, c='magenta', s=10, label="CAMMPCAN V3") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V4
    lons8 = V4_17_18_Long
    lats8 = V4_17_18_Lat
    x8,y8 = m(lons8, lats8)
    obs = V4_17_18_Hg0
    m.scatter(x8, y8, c='cyan', s=10, label="CAMMPCAN V4") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    
    # CAMMPCAN (2018/19)
    #V1
    #lons9 = V1(18_19)_Long
    #lats9 = V1(18_19)_Lat
    #x9,y9 = m(lons9, lats9)
    #obs = V1(18_19)_Hg0
    #m.scatter(x9, y9, c=obs, vmin=0.4, vmax=1.6, cmap=cm.plasma, s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V2
    #lons10 = V1(18_19)_Long
    #lats10 = V1(18_19)_Lat
    #x10,y10 = m(lons10, lats10)
    #obs = V1(18_19)_Hg0
    #m.scatter(x10, y10, c=obs, vmin=0.4, vmax=1.6, cmap=cm.plasma, s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V3
    #lons11 = V1(18_19)_Long
    #lats11 = V1(18_19)_Lat
    #x11,y11 = m(lons11, lats11)
    #obs = V1(18_19)_Hg0
    #m.scatter(x11, y11, c=obs, vmin=0.4, vmax=1.6, cmap=cm.plasma, s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    #V4
    #lons12 = V1(18_19)_Long
    #lats12 = V1(18_19)_Lat
    #x12,y12 = m(lons12, lats12)
    #obs = V1(18_19)_Hg0
    #m.scatter(x12, y12, c=obs, vmin=0.4, vmax=1.6, cmap=cm.plasma, s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
   
    # RV Xuelong
    lons13 = RVXue_Long
    lats13 = RVXue_Lat
    x13,y13 = m(lons13, lats13)
    obs = RVXue_Hg0
    m.scatter(x13, y13, c='brown', s=10, label="RV Xuelong ") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  

    # CHINARE
    lons14 = CHINARE_Long
    lats14 = CHINARE_Lat
    x14,y14 = m(lons14, lats14)
    obs = CHINARE_Hg0
    m.scatter(x14, y14, c='orange', s=10, label="CHINARE") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    
    # OSO_V1
    lons15 = OSO_V1_Long
    lats15 = OSO_V1_Lat
    x15,y15 = m(lons15, lats15)
    #obs = OSO_V1_Hg0
    m.scatter(x15, y15, c='purple', s=10, label="OSO 10/11") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  

    # OSO_V2
    lons16 = OSO_V2_Long
    lats16 = OSO_V2_Lat
    x16,y16 = m(lons16, lats16)
    #obs = OSO_V2_Hg0
    m.scatter(x16, y16, c='purple', s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)

    # OSO_V3
    lons17 = OSO_V3_Long
    lats17 = OSO_V3_Lat
    x17,y17 = m(lons17, lats17)
    #obs = OSO_V3_Hg0
    m.scatter(x17, y17, c='purple', s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)

    # ANTXXIX_V6
    lons18 = ANTXXIX_V6_Long
    lats18 = ANTXXIX_V6_Lat
    x18,y18 = m(lons18, lats18)
    #obs = ANTXXIX_V6_Hg0
    m.scatter(x18, y18, c='darkgoldenrod', s=10, label="ANTXXIX_V6-7") # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)

    # ANTXXIX_V7
    lons19 = ANTXXIX_V7_Long
    lats19 = ANTXXIX_V7_Lat
    x19,y19 = m(lons19, lats19)
    #obs = ANTXXIX_V7_Hg0
    m.scatter(x19, y19, c='darkgoldenrod', s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)
     
#    # PLOT THE CAMPCANN VOYAGES
#    # V1 (Hobart to Davis)
#    lons5 = [147.1930,77.5803]
#    lats5 = [-42.5250,-68.3436]
#    # V2 (Hobart to Casey)
#    lons6 = [147.1930,110.3136]
#    lats6 = [-42.5250,-66.1657]
#    # V3 (Hobart to Mawson)
#    lons7 = [147.1930,62.5227]
#    lats7 = [-42.5250,-67.3612]
#    # V4 (Hobart to Macquarie Island)
#    lons8 = [147.1930,158.5609]
#    lats8 = [-42.5250,-54.2959]
#    
#    x5, y5 = m(lons5,lats5)
#    x6, y6 = m(lons6,lats6)
#    x7, y7 = m(lons7,lats7) 
#    x8, y8 = m(lons8,lats8)
#    plt.plot(x5,y5,x6,y6,x7,y7,x8,y8, color='grey')
    
    #---------------------------
    # PLOT LABELS FOR THE MONITORING SITES
    labels = ['Cape Grim', 'Amsterdam \n Island', 'Cape Point', 'Dumont d\'Urville', 'Concordia', 'Troll', 'Bariloche', 'McMurdo'] # text you want displayed
    x_offsets = [-200000, -450000, -450000, -450000, -450000, 450000, -450000, -450000] # Adjust the label positions on the x-axis
    y_offsets = [175000,  -900000, -350000, -350000, -350000, -350000, -350000, -350000] # Adjust the label positions on the y-axis
    for label, xpt, ypt, x_offset, y_offset in zip(labels, x1, y1, x_offsets, y_offsets):
        plt.text(xpt+x_offset, ypt+y_offset, label, fontsize=10, color='black', fontweight='bold') # define how want the labels displayed (color, fontsize, etc...)   
    
    # PLOT THE LABELS FOR THE AAD STATIONS
    labels = ['Mawson', 'Davis', 'Casey', 'Macquarie \n Island'] # text you want displayed
    x_offsets = [-1000000, -450000, -450000, -450000] # Adjust the label positions on the x-axis
    y_offsets = [-350000, -350000, -350000, -500000] # Adjust the label positions on the y-axis
    for label, xpt, ypt, x_offset, y_offset in zip(labels, x2, y2, x_offsets, y_offsets):
        plt.text(xpt+x_offset, ypt+y_offset, label, fontsize=10, color='black', fontweight='bold') # define how want the labels displayed (color, fontsize, etc...)   
    
    # PLOT THE LABELS FOR THE X-AXIS and Y-AXIS
    plt.xlabel('Longitude [$^\circ$]', fontsize=15, labelpad=30) # (text displayed, fontsize, distance between map and label) 
    plt.ylabel('Latitude [$^\circ$]', fontsize=15, labelpad=50)  # (text displayed, fontsize, distance between map and label) 
    
    # DEFINE FONTSIZE FOR THE MAP TITLE
    plt.title(tit, fontsize=25)
    return m

#------------------------------------------------------------------------------
# START PLOTTING THE MAP
fig = plt.figure()

# CHOOSE HOW MANY FIGURES YOU WISH TO PLOT
ax1 = plt.subplot(111) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
ax1 = plotkg(Hg0_Av[:,:], 'Y', "Study area and ship tracks") # Hg$^0$ Note that Hg0 is an array Hg0['time','lev','lat','lon']
ax1 = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#ax1 = plt.annotate("SIPEXII", xy=(-110,-20), color='red', fontweight='bold')
#ax1 = plt.subplot(111) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
#ax1 = plotkg(BrO_Av[:,:], 'Y', "Simulated BrO Concentrations (Sept & Oct 2012)") # Note that Hg0 is an array Hg0['time','lev','lat','lon']
#ax1 = plt.subplot(111) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
#ax1 = plotkg(Br_Av[:,:], 'Y', "Simulated Br Concentrations (Sept & Oct 2012)") # Note that Hg0 is an array Hg0['time','lev','lat','lon']
#ax1 = plt.subplot(111) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
#ax1 = plotkg(BrO_Av[:,:]+Br_Av[:,:], 'Y', "Simulated BrO + Br Concentrations (Sept & Oct 2012)") # Note that Hg0 is an array Hg0['time','lev','lat','lon']