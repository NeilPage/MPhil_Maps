#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:23:17 2018

@author: ncp532
"""

# File system packages
from netCDF4 import Dataset, MFDataset
import xarray as xr

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates            
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from matplotlib import gridspec
import cmocean
import matplotlib.image as mpimg
from matplotlib import cm                   # imports the colormap function from matplotlib
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy is great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats
import os
import glob

# Date and Time handling package
import datetime as dt
from datetime import datetime,time, timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# OBSERVATIONS

#-------------
# Sea Ice
#-------------
SeaIce_V1_17A   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181107.nc', decode_cf=False, engine='netcdf4')
SeaIce_V1_17D   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181115.nc', decode_cf=False, engine='netcdf4')

SeaIce_V2_17A   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181215.nc', decode_cf=False, engine='netcdf4')
SeaIce_V2_17D   = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20181230.nc', decode_cf=False, engine='netcdf4')

SeaIce_V3_17A  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20190130.nc', decode_cf=False, engine='netcdf4')
SeaIce_V3_17D  = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/OISST_SeaIce/2018-19/oisst-avhrr-v02r01.20190209.nc', decode_cf=False, engine='netcdf4')

SeaIce_PCAN    = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20170126_median5day.nc', decode_cf=False, engine='netcdf4')

#------------------------------------------------------------------------------
# SET UP THE VALUES TO PLOT

# Sea ice
ice_V1_17A = SeaIce_V1_17A.variables['ice'][0,0,:,:]    # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)
ice_V1_17D = SeaIce_V1_17A.variables['ice'][0,0,:,:]    # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)

ice_V2_17A = SeaIce_V2_17A.variables['ice'][0,0,:,:]    # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)
ice_V2_17D = SeaIce_V2_17D.variables['ice'][0,0,:,:]    # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)

ice_V3_17A = SeaIce_V3_17A.variables['ice'][0,0,:,:]   # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)
ice_V3_17D = SeaIce_V3_17D.variables['ice'][0,0,:,:]   # Sea Ice concentration (time,lev,y,x)(1,1,720,1440)

land_PCAN  = SeaIce_PCAN.variables['land'][0,:,:] 

# Latitudes
lats_V1_17A = SeaIce_V1_17A.lat
lats_V1_17D = SeaIce_V1_17D.lat

lats_V2_17A = SeaIce_V2_17A.lat
lats_V2_17D = SeaIce_V2_17D.lat

lats_V3_17A = SeaIce_V3_17A.lat
lats_V3_17D = SeaIce_V3_17D.lat

lats_PCAN   = SeaIce_PCAN.latitude

# Longitudes
lons_V1_17A = SeaIce_V1_17A.lon
lons_V1_17D = SeaIce_V1_17D.lon

lons_V2_17A = SeaIce_V2_17A.lon
lons_V2_17D = SeaIce_V2_17D.lon

lons_V3_17A = SeaIce_V3_17A.lon
lons_V3_17D = SeaIce_V3_17D.lon

lons_PCAN   = SeaIce_PCAN.longitude

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

#------------------------------------------------------------------------------
# PLOT THE FIGURE
#------------------------------------------------------------------------------
fig = plt.figure()
#fig.suptitle('3 January 2018', fontsize=20, y=0.95)
gs = fig.add_gridspec(ncols=2, nrows=3)
plt.subplots_adjust(wspace=0.05)

#-------------------------------------
# Graph 1 (V1 Arrival)

ax = plt.subplot(gs[0,0], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)
#ax = plt.subplot(gs[0,0], projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.))

# SET UP THE PLOT
#polar_extent = [0, 180, -90, -55]
#ax.set_extent(polar_extent, crs=ccrs.PlateCarree())
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V1_17A = np.ma.masked_where(ice_V1_17A==-999,ice_V1_17A)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V1_17A, lats_V1_17A, seaice_data_V1_17A, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)

# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)')#, rotation=90)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.75)
ax.text(-0.35, 0.5, "        V1 (Davis)        ", transform=ax.transAxes, fontsize=14, verticalalignment='center', bbox=props, rotation=90)

props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "7 Nov 2018", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

# Plot title for arrival column
plt.title('Arrival', fontsize=20, y=1.10)

#-------------------------------------
# Graph 2 (V1 Departure)

ax = plt.subplot(gs[0,1], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)

# SET UP THE PLOT
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V1_17D = np.ma.masked_where(ice_V1_17D==-999,ice_V1_17D)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V1_17D, lats_V1_17D, seaice_data_V1_17D, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)
# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)')#, rotation=90)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "15 Nov 2018", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

# Plot title for departure column
plt.title('Departure', fontsize=20, y=1.10)

#-------------------------------------
# Graph 3 (V2 Arrival)

ax = plt.subplot(gs[1,0], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)

# SET UP THE PLOT
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V2_17A = np.ma.masked_where(ice_V2_17A==-999,ice_V2_17A)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V2_17A, lats_V2_17A, seaice_data_V2_17A, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)

# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)')#, rotation=90)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.75)
ax.text(-0.35, 0.5, "        V2 (Casey)        ", transform=ax.transAxes, fontsize=14, verticalalignment='center', bbox=props, rotation=90)

props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "15 Dec 2018", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

#-------------------------------------
# Graph 4 (V2 Departure)

ax = plt.subplot(gs[1,1], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)

# SET UP THE PLOT
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V2_17D = np.ma.masked_where(ice_V2_17D==-999,ice_V2_17D)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V2_17D, lats_V2_17D, seaice_data_V2_17D, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)

# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)')#, rotation=90)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "30 Dec 2018", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

#-------------------------------------
# Graph 5 (V3 Arrival)

ax = plt.subplot(gs[2,0], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)

# SET UP THE PLOT
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V3_17A = np.ma.masked_where(ice_V3_17A==-999,ice_V3_17A)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V3_17A, lats_V3_17A, seaice_data_V3_17A, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)

# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)')#, rotation=90)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='wheat', alpha=0.75)
ax.text(-0.35, 0.5, "V3 (Mawson & Davis)", transform=ax.transAxes, fontsize=14, verticalalignment='center', bbox=props, rotation=90)

props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "30 Jan 2019", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

#-------------------------------------
# Graph 6 (V3 Departure)

ax = plt.subplot(gs[2,1], projection=ccrs.SouthPolarStereo())#PlateCarree()) # options graph 1 (vertical no, horizontal no)

# SET UP THE PLOT
ax.set_extent([-45, 135, -42.5, -90])#, crs=ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

# PLOT THE DATA (SEA ICE CONCENTRATION) 
seaice_data_V3_17D = np.ma.masked_where(ice_V3_17D==-999,ice_V3_17D)
#cmap=cm.get_cmap('viridis')
cmap=cmocean.cm.ice
cmap.set_bad(color='lightgrey')
cs = ax.pcolormesh(lons_V3_17D, lats_V3_17D, seaice_data_V3_17D, transform=data_crs, vmin=0, vmax=100, cmap=cmap) #, bins=np.arange(0,100, 10))
#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cmocean.cm.ice) #, bins=np.arange(0,100, 10))

# PLOT LAND
land_PCAN = np.ma.masked_where(land_PCAN==0,land_PCAN)
cmap=cmocean.cm.ice_r
ax.pcolormesh(lons_PCAN, lats_PCAN, land_PCAN,transform=data_crs,cmap=cmap)

# PLOT THE BACK TRAJECTORIES
cmap = plt.cm.autumn_r
#norm = BoundaryNorm(np.arange(-120,0,10), cmap.N)
norm = BoundaryNorm(np.arange(0,90,10), cmap.N)

#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['Traj Age'], label='Traj Height (m)')
#ax.scatter(Traj['Traj Lon'], Traj['Traj Lat'], zorder=2, cmap=cmap, transform=data_crs, marker='o', s=1, norm=norm, c=Traj['IceContact_100m'], label='Traj Height (m)')

# PLOT THE AAD STATIONS (Lon, Lat)
Davis_lon,   Davis_lat   = -68.5766, 77.9674
Mawson_lon,  Mawson_lat  = -67.6027, 62.8738
Casey_lon,   Casey_lat   = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='r', markeredgecolor='k', ms=10, marker='*')
ax.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='y', markeredgecolor='k', ms=10, marker='*')
#ax.plot(SIPEXII_lat, SIPEXII_lon, transform=data_crs, color='k', marker='o')

# Plot the marker labels
ax.text(Davis_lat + 3,  Davis_lon  - 2, 'Davis',  transform=data_crs, horizontalalignment='right')
ax.text(Mawson_lat + 3, Mawson_lon - 2, 'Mawson', transform=data_crs, horizontalalignment='right')
ax.text(Casey_lat + 3,  Casey_lon  - 2, 'Casey',  transform=data_crs, horizontalalignment='right')
#ax.text(SIPEXII_lat + 3, SIPEXII_lon - 2, 'SIPEXII',horizontalalignment='right')

# PLOT THE MAP GRIDLINES
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='-')
#gl.xlocator   = ticker.FixedLocator([-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,110,120,130,140,150,160,170,180])
#gl.ylocator   = ticker.FixedLocator([-35,-40,-45,-50,-55,-60,-65,-70,-75,-80,-85])
gl.xlocator   = ticker.FixedLocator([-180,-90,0,90,180])
gl.ylocator   = ticker.FixedLocator([-60,-90])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (15/12/2018)", y=1.1, fontsize=20)

# ax.text(-0.12, 0.55, 'latitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='vertical', rotation_mode='anchor',
#          transform=ax.transAxes)

# ax.text(0.5, -0.4, 'longitude [$^\circ$]', fontsize=10, va='bottom', ha='center',
#          rotation='horizontal', rotation_mode='anchor',
#          transform=ax.transAxes)

# adjust the axis labels and ticks
ax.xaxis.labelpad = 30
ax.yaxis.labelpad = 30
ax.tick_params(labelsize=10)

# Text box in upper left
props = dict(boxstyle='round', facecolor='white', alpha=1.0)
ax.text(0.05, 0.105, "9 Feb 2019", transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)


#-------------------------------------
# Add a colorbar for the figure
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
cb = fig.colorbar(cs, cax=cbar_ax, ticks=[0,10,20,30,40,50,60,70,80,90,100], pad = 0.08, shrink=.995)#, orientation="horizontal")
cb.set_label('Sea ice cover (%)')#, rotation=90)

custom_lines = [Line2D([0], [0], color='lightgrey', markeredgecolor='black', lw=10)]
fig.legend(custom_lines, ['No data'], loc='upper left', bbox_to_anchor=(0.8455, 0.130), fontsize=10, frameon=False)

plt.tight_layout()