#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:22:07 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
#from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
from netCDF4 import Dataset,MFDataset
import xarray as xr

# DATA HANDLING PACKAGES
import numpy as np
import matplotlib.pyplot as plt

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm                   # imports the colormap function from matplotlib
from matplotlib.patches import Rectangle, Polygon
import matplotlib.ticker as mticker
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# HAMBURG
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20180112_median5day.nc')

# NSIDC
# 2017
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2017-18/seaice_conc_daily_sh_f17_20180103_v03r01.nc')
# 2018
#cubes = Dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2018-19/seaice_conc_daily_icdr_sh_f18_20181230_v01r00.nc')

# XARRAY
#cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/Hamburg_ICDC/20180215_median5day.nc',decode_cf=False,engine='netcdf4')
cubes = xr.open_dataset('/Users/ncp532/Documents/Data/Sea_Ice_Cover/NSIDC/2012/seaice_conc_daily_sh_f17_20121111_v03r01.nc',decode_cf=False,engine='netcdf4')

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

# NETCDF
#lats = cubes.variables['latitude'][:] # (y,x) (664,632)
#lons = cubes.variables['longitude'][:] # (y,x) (664,632)
# NSIDC
#seaice_data = cubes.variables['seaice_conc_cdr'][0,:,:]*100 # Sea Ice concentration (time,y,x)(1,664,632)
# HAMBURG
#seaice_data = cubes.variables['sea_ice_area_fraction'][0,:,:] # Sea Ice concentration (time,y,x)(1,664,632)

# XARRAY
lats = cubes.latitude
lons = cubes.longitude
#land = cubes.land
#seaice_data = cubes.sea_ice_area_fraction[0,:,:] # Hamburg
seaice_data = cubes.seaice_conc_cdr[0,:,:] # Hamburg

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

#------------------------------------------------------------------------------
# START PLOTTING THE MAP
fig = plt.figure(figsize=(10, 6))
plt.subplots_adjust()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([40, 140, -55, -75])#, crs=ccrs.PlateCarree())
#ax2 = ax.twinx()
#ax2.set_extent([40, 120, -50, -75])#, crs=ccrs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='-')
gl.xlocator = mticker.FixedLocator([40,50,60,70,80,90,100,110,110,120,130,140,150,160])
gl.ylocator = mticker.FixedLocator([-55,-60,-65,-70,-75,-80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 15, 'color': 'black'}
gl.ylabel_style = {'size': 15, 'color': 'black'}

#------------------------------------
# PLOT THE DATA (SEA ICE CONCENTRATION) 

#cs = ax.pcolormesh(lons, lats, seaice_data, transform=data_crs, cmap=cm.get_cmap('jet',20)) #, bins=np.arange(0,100, 10))
#cb = fig.colorbar(cs,ticks=[0,10,20,30,40,50,60,70,80,90,100],shrink=.50)
#cb = fig.colorbar(cs,ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],shrink=.50)
#------------------------------------
# PLOT THE AAD STATIONS

# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
SIPEXII_lon, SIPEXII_lat = -61.5205, 121.1855

# Plot the station markers
ax.plot(Davis_lat, Davis_lon, color='k', marker='o')
ax.plot(Mawson_lat, Mawson_lon, color='k', marker='o')
ax.plot(Casey_lat, Casey_lon, color='k', marker='o')
ax.plot(SIPEXII_lat, SIPEXII_lon, color='k', marker='o')

# Plot the mareker labels
ax.text(Davis_lat + 8, Davis_lon - 3.5, 'Davis (V1)',fontsize=25,horizontalalignment='right')
ax.text(Mawson_lat + 5, Mawson_lon - 2.5, 'Mawson (V3)',fontsize=25,horizontalalignment='right')
ax.text(Casey_lat + 6, Casey_lon - 2.5, 'Casey (V2)',fontsize=25,horizontalalignment='right')
ax.text(SIPEXII_lat + 7, SIPEXII_lon - 2.5, 'SIPEXII',fontsize=25,horizontalalignment='right')

#------------------------------------
# PLOT A GRID BOX

#-----------
# ACTUAL (LON/LAT)
#-----------
# Davis
rectDA = Rectangle((77.5 - 1.25 ,-68 -1),2.5,2,linewidth=1,edgecolor='r',facecolor='r',alpha=0.75)
ax.add_patch(rectDA)

# Casey
rectCA = Rectangle((110.0 - 1.25 ,-66 -1),2.5,2,linewidth=1,edgecolor='r',facecolor='r',alpha=0.75)
ax.add_patch(rectCA)

# Mawson
rectMA = Rectangle((62.5 - 1.25 ,-68 -1),2.5,2,linewidth=1,edgecolor='r',facecolor='r',alpha=0.75)
ax.add_patch(rectMA)

# SIPEXII
rectSIPA = Rectangle((120.0 - 1.25 ,-62 -1),2.5,2,linewidth=1,edgecolor='r',facecolor='r',alpha=0.75)
ax.add_patch(rectSIPA)

#-----------
# SUGGESTED (LON/LAT)
#-----------
# Davis
rectDS = Rectangle((77.5 - 1.25 ,-66 -1),2.5,2,linewidth=1,edgecolor='g',facecolor='g',alpha=0.75)
ax.add_patch(rectDS)

# Casey
rectCS = Rectangle((110.0 - 1.25 ,-64 -1),2.5,2,linewidth=1,edgecolor='g',facecolor='g',alpha=0.75)
ax.add_patch(rectCS)

# Mawson
rectMS = Rectangle((62.5 - 1.25 ,-66 -1),2.5,2,linewidth=1,edgecolor='g',facecolor='g',alpha=0.75)
ax.add_patch(rectMS)

# SIPEXII
rectSIPS = Rectangle((122.5 - 1.25 ,-62 -1),2.5,2,linewidth=1,edgecolor='g',facecolor='g',alpha=0.75)
ax.add_patch(rectSIPS)

#------------------------------------
# PLOT TITLE, AXIS LABEL & LEGEND TITLE
#plt.title("Sea Ice Cover (11/11/2012)", y=1.1, fontsize=20)
#cb.set_label('Concentration (%)', rotation=90)

ax.text(-0.05, 0.55, 'latitude [$^\circ$]', fontsize=25, va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

ax.text(0.5, -0.2, 'longitude [$^\circ$]', fontsize=25, va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)
plt.show()