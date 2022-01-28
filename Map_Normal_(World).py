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

# Simulations
#dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Dynamic_Ocean/trac_avg.geosfp_2x25_Hg.2013*.nc') # Dynamic Ocean Hg Data
#dataset1b = MFDataset('/Users/macbook/Documents/UOW/Some_scripts_/nc_Dynamic_Ocean/GEOSChem.StateMet.2013*.nc4') # Dynamic Ocean Met Data
#dataset2a = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.SpeciesConc.2013*.nc4') # Offline Ocean Hg Data
#dataset2b = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Offline_Ocean/GEOSChem.StateMet.2013*.nc4') # Offline Ocean Met Data

# Observations
#SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/SIPEXII_2012/SIPEXII_Hg_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for SIPEXII (2012)
#V1 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V1 (2017/18)
#V2 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V2 (2017/18)
#V3 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V3 (2017/18)
#V4 = pd.read_csv('/Users/ncp532/Documents/Data/SIPEX_II_Mercury_Air/SIPEXII_Hg0_Lat_Long.csv') # Hg0 data for CAMPCANN V4 (2017/18)

# SIPEXII_Long = np.array(SIPEXII['LONGITUDE']) # Longitude for SIPEXII
# SIPEXII_Lat = np.array(SIPEXII['LATITUDE']) # Latitude for SIPEXII
# SIPEXII_Hg0 = np.array(SIPEXII['ng/m3']) # Hg0 for SIPEXII

#------------------------------------------------------------------------------
# GET THE VARIABLES YOU WOULD LIKE FROM THE GEOS-CHEM SIMULATION
# you can display the variable names by typing "dataset1" into the python console
# you can display more information about individual variables by typing "dataset1.variables['INSERT THE VARIABLES NAME']" into the python console
# e.g. for Hg0 type "dataset1.variables['IJ_AVG_S__Hg0']" into the python console

# Dynamic Ocean
# airden1 = dataset1.variables['BXHGHT_S__AIRNUMDE'][:] # Air Density
# lat1    = dataset1.variables['lat'][:] # in degrees
# lon1    = dataset1.variables['lon'][:] # in degrees

# R       = 8.3145 # universal gas constant (J/mol K)
# MM_Hg   = 200.59 # Molar Mass Hg (g/mol)
# MM_Air  = 28.97 # Molar Mass Dry Air (g/mol)
# Avo     = 6.0221*1e23 # Avogadro's number

# Hg0     = (dataset1.variables['IJ_AVG_S__Hg0'][:]) # dry mixing ratio of Hg0 (mol/mol dry)
# #Hg0     = Hg0 * 1e12 * (MM_Hg / MM_Air) * airden1
# Hg0     = Hg0 * 1e3 * MM_Hg / Avo * airden1

# HgII    = (dataset1.variables['IJ_AVG_S__Hg2'][:]) # dry mixing ratio of Hg0 (mol/mol dry)
# HgII    = HgII*1e12*(MM_Hg/MM_Air)*airden1

# HgP     = (dataset1.variables['IJ_AVG_S__HgP'][:]) # dry mixing ratio of Hg0 (mol/mol dry)
# HgP     = HgP*1e12*(MM_Hg/MM_Air)*airden1


#Hg_Av   = (Hg1[8,0,:,:] + Hg1[9,0,:,:]) / 2 # Mean simulated Hg0 concentration for Sept and Oct

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
    m=Basemap(llcrnrlat=-90,  urcrnrlat=90,          # set the latitude limits
              llcrnrlon=-180, urcrnrlon=180,         # set the longitude limits
              resolution='c',projection='gall')      # choose the resolution and projection type
    
    # SOUTH POLAR STEROGRAPHIC PROJECTION          
    #m = Basemap(projection='spstere',boundinglat=-10,lon_0=90,resolution='c')          
              
    mlons,mlats = m(lons, lats) # convert lons/lat from middle of grid-box to to corner of grid-box
    cs=m.pcolormesh(mlons, mlats, x, latlon=False, vmin=0.4, vmax=1.6, cmap=cm.plasma)   # define the colorbar attributes (colors used, add "_r" to flip the color scale)    
    m.drawcoastlines() # draw coastlines on the map
    
    # FORMAT THE COLORBAR
    if bar=='Y': # Y if you want a colorbar, N if you don't
        cb=m.colorbar(cs,"bottom",size="5%", pad="20%") # Position the colorbar next to your map and define its size
        cb.set_label(r"$ng$ $m$$^-$$^3$", fontsize =20) # add text for the colourbar label and define fontsize
        cb.ax.tick_params(labelsize=20) # add numbers to the colorbar
        cb.ax.minorticks_on()           # add ticks to the colorbar
    
    # PLOT PARALLELS AND MERIDANS
    parallels = np.arange(-90,10, 10.0)  # set the latitude limits and how often you wish parallels to be plotted
    m.drawparallels(parallels,labels=[True,False,True,False], fontsize=15) # define the fontsize for parallels
    meridians = np.arange(-180,179,30.0) # set the longitude limits and how often you wish meridans to be plotted
    m.drawmeridians(meridians,labels=[True,False,False,True], fontsize=15) # define the fontsize for meridans
    
    # PLOT LOCATION MARKERS FOR MONITORING SITES
    lons1 = [144.6832, 77.3423, 18.4897, 140.0007, 123.1956, 2.3206, -71.2512, 166.4006]    # [Cape Grim, Amsterdam Island, Cape Point, Mawson, Davis, Casey, Macquarie Island, Dumont d'Urville, Concordia, Troll, Bariloche, McMurdo]
    lats1 = [-40.6832, -37.4748, -34.3535, -66.3946, -75.0559, -72.0041, -41.7438, -77.5047]  # [Cape Grim, Amsterdam Island, Cape Point, Mawson, Davis, Casey, Macquarie Island, Dumont d'Urville, Concordia, Troll, Bariloche, McMurdo]
    x1,y1 = m(lons1, lats1)
    
    # OPTION 1: if you want plain location markers use this line
    m.scatter(x1, y1, c='white', edgecolors='black', s=50) # define the location markers (color, edgecolor, size)
    
    # OPTION2: if you want the location markers to represent observed Hg concentrations use these two lines
    #obs = [0.909, 0.974, 0.808, 1.032, 1.092] # define the observed Hg concentration for each site
    
    # PLOT THE LOCATION MARKERS FOR AAD STATIONS
    lons2 = [62.5227, 77.5803, 110.3136, 158.5609]    # [Mawson, Davis, Casey, Macquarie Island]
    lats2 = [-67.3612, -68.3436, -66.1657, -54.2959]  # [Mawson, Davis, Casey, Macquarie Island]
    x2,y2 = m(lons2, lats2)
    
    # OPTION 1: if you want plain location markers use this line
    m.scatter(x2, y2, c='black', edgecolors='black', s=50) # define the location markers (color, edgecolor, size)
    
    # # PLOT THE TRACK FOR SIPEXII
    # lons3 = SIPEXII_Long
    # lats3 = SIPEXII_Lat
    # x3,y3 = m(lons3, lats3)
    
    # # OPTION 1: if you want plain location markers use this line
    # #m.scatter(x3, y3, c='white', edgecolors='black', s=200) # define the location markers (color, edgecolor, size)
    
    # # OPTION2: if you want the location markers to represent observed Hg concentrations use these two lines
    # #obs = [0.909, 0.974, 0.808, 1.032, 1.092] # define the observed Hg concentration for each site
    # obs = SIPEXII_Hg0
    # m.scatter(x3, y3, c=obs, vmin=0.4, vmax=1.6, cmap=cm.plasma, s=10) # define the location markers (color=obs, you need to set vmin/vmax and the cmap to same values as the colorbar above)  
    
    # PLOT LABELS FOR THE MONITORING SITES
    labels = ['Cape Grim', 'Amsterdam \n Island', 'Cape Point', 'Dumont d\'Urville', 'Concordia', 'Troll', 'Bariloche', 'McMurdo'] # text you want displayed
    x_offsets = [-200000, -450000, -450000, -450000, -450000, 450000, -450000, -450000] # Adjust the label positions on the x-axis
    y_offsets = [175000,  -500000, -350000, -350000, -350000, -350000, -350000, -350000] # Adjust the label positions on the y-axis
    for label, xpt, ypt, x_offset, y_offset in zip(labels, x1, y1, x_offsets, y_offsets):
        plt.text(xpt+x_offset, ypt+y_offset, label, fontsize=10, color='white', fontweight='bold') # define how want the labels displayed (color, fontsize, etc...)   
    
    # PLOT THE LABELS FOR THE AAD STATIONS
    labels = ['Mawson', 'Davis', 'Casey', 'Macquarie \n Island'] # text you want displayed
    x_offsets = [-1000000, -450000, -450000, -450000] # Adjust the label positions on the x-axis
    y_offsets = [-350000, -350000, -350000, -500000] # Adjust the label positions on the y-axis
    for label, xpt, ypt, x_offset, y_offset in zip(labels, x2, y2, x_offsets, y_offsets):
        plt.text(xpt+x_offset, ypt+y_offset, label, fontsize=10, color='black', fontweight='bold') # define how want the labels displayed (color, fontsize, etc...)   
    
    # PLOT THE CAMPCANN VOYAGES
    # V3 (Hobart to Mawson)
    lons4 = [147.1930,62.5227]
    lats4 = [-42.5250,-67.3612]
    # V1 (Hobart to Davis)
    lons5 = [147.1930,77.5803]
    lats5 = [-42.5250,-68.3436]
    # V2 (Hobart to Casey)
    lons6 = [147.1930,110.3136]
    lats6 = [-42.5250,-66.1657]
    # V4 (Hobart to Macquarie Island)
    lons7 = [147.1930,158.5609]
    lats7 = [-42.5250,-54.2959]
    
    x4, y4 = m(lons4,lats4)
    x5, y5 = m(lons5,lats5)
    x6, y6 = m(lons6,lats6)
    x7, y7 = m(lons7,lats7) 
    plt.plot(x4,y4,x5,y5,x6,y6,x7,y7, color='grey')
    
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
#ax1 = plotkg(Hg0[0,:,:], 'Y', "Simulated & Observed GEM Concentrations (Sept & Oct 2012)") # Hg$^0$ Note that Hg0 is an array Hg0['time','lev','lat','lon']
#ax1 = plt.subplot(312) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
#ax1 = plotkg(Hg1[8,0,:,:], 'Y', "Hg0 Sept") # Note that Hg0 is an array Hg0['time','lev','lat','lon']
#ax1 = plt.subplot(313) # subplot(number of figures vertical, number of figures horizontal, total number of figures)
#ax1 = plotkg(Hg1[9,0,:,:], 'Y', "Hg0 Oct") # Note that Hg0 is an array Hg0['time','lev','lat','lon']