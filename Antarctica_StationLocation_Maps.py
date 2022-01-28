#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:32:21 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
import xarray as xr

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# DATA HANDLING PACKAGES
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# DRAWING PACKAGES
import cartopy.crs as ccrs
import cartopy
from matplotlib import cm                   # imports the colormap function from matplotlib
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt             
import matplotlib.dates as mdates            
from matplotlib.ticker import MaxNLocator

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
# START PLOTTING THE MAP

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

# We will view the map on a polar stereo coordinate system
projection = ccrs.SouthPolarStereo(central_longitude=0)

plt.close()
fig = plt.figure(figsize=(10, 10))

#--------------------------------------------
# GRAPH 1 (ANTARCTICA MAP)
ax1 = plt.subplot(111, projection=projection) # (vertical no, horizontal no, graph no)

#ax.set_global()
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax1.coastlines()
gl = ax1.gridlines(color='gray',alpha=0.5,)
ax1.set_extent([-180, 180, -59, -90], crs=ccrs.PlateCarree())

# # PLOT PARALLELS AND MERIDANS
# parallels = np.arange(-90,10, 10.0)  # set the latitude limits and how often you wish parallels to be plotted
# m.drawparallels(parallels,labels=[True,False,True,False], fontsize=15) # define the fontsize for parallels
# meridians = np.arange(-180,179,30.0) # set the longitude limits and how often you wish meridans to be plotted
# m.drawmeridians(meridians,labels=[True,False,False,True], fontsize=15) # define the fontsize for meridans

#---------------------------
# Station location (Lon, Lat)
Davis_lat,     Davis_lon     = -68.5762, 77.9696
Mawson_lat,    Mawson_lon    = -67.6033, 62.8742
Casey_lat,     Casey_lon     = -66.2821, 110.5285
DuUrv_lat,     DuUrv_lon     = -66.3946, 140.0007
Concor_lat,    Concor_lon    = -75.0559, 123.1956
Troll_lat,     Troll_lon     = -72.0041, 2.3206
McMurdo_lat,   McMurdo_lon   = -77.5047, 166.4006
Neumayer_lat,  Neumayer_lon  = -70.6770, -8.2720
TerraNB_lat,   TerraNB_lon   = -74.6940, 164.114
Halley_lat,    Halley_lon    = -75.5674, -25.5165
Zhongshan_lat, Zhongshan_lon = -69.3700, 76.3800
Vostok_lat,    Vostok_lon    = -78.4645, 106.8339
Amun_Scot_lat, Amun_Scot_lon = -90.0000, -139.2667
Belgrano_lat,  Belgrano_lon  = -77.8736, -34.6267
Marambio_lat,  Marambio_lon  = -64.2501, -56.6420

lons = [Davis_lon, Mawson_lon, Casey_lon, DuUrv_lon, Concor_lon, Troll_lon, McMurdo_lon, Neumayer_lon, TerraNB_lon, Halley_lon, Zhongshan_lon, Vostok_lon, Amun_Scot_lon, Belgrano_lon, Marambio_lon]
lats = [Davis_lat, Mawson_lat, Casey_lat, DuUrv_lat, Concor_lat, Troll_lat, McMurdo_lat, Neumayer_lat, TerraNB_lat, Halley_lat, Zhongshan_lat, Vostok_lat, Amun_Scot_lat, Belgrano_lat, Marambio_lat]

ax1.scatter(lons, lats, c='red', transform=data_crs, zorder=2, vmin=0.4, vmax = 1.6, cmap=cm.get_cmap('jet',12), edgecolors='black', s=50)

# Plot the marker labels
ax1.text(Davis_lon + 2,     Davis_lat + 2.5,   'DV',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Mawson_lon + 2,    Mawson_lat - 1,    'MW',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Casey_lon + 2,     Casey_lat - 1,     'CY',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(DuUrv_lon,         DuUrv_lat +2.5,    'DDU', transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Concor_lon + 2,    Concor_lat - 1,    'DC',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Troll_lon + 2,     Troll_lat - 2,     'TR',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(McMurdo_lon + 4,   McMurdo_lat + 1,   'MM',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Neumayer_lon - 2,  Neumayer_lat,      'NM',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(TerraNB_lon +2,    TerraNB_lat,       'TNB', transform=data_crs, horizontalalignment='right', fontsize=12)

ax1.text(Halley_lon  -2,    Halley_lat      ,  'HA',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Zhongshan_lon + 2, Zhongshan_lat - 1, 'ZG',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Vostok_lon + 1,    Vostok_lat - 1,    'VK',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Amun_Scot_lon + 4, Amun_Scot_lat + 1, 'AS',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Belgrano_lon+2 ,  Belgrano_lat +1,  'BE',  transform=data_crs, horizontalalignment='right', fontsize=12)
ax1.text(Marambio_lon +2,   Marambio_lat - 3,  'MB',  transform=data_crs, horizontalalignment='right', fontsize=12)
