#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 15:08:51 2019

@author: ncp532
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:32:21 2019

@author: ncp532
"""
# FILE SYSTEM PACKAGES
from netCDF4 import Dataset,MFDataset				# function used to open multiple netcdf files
import xarray as xr

# DATA HANDLING PACKAGES
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# DRAWING PACKAGES
import cartopy.crs as ccrs
from matplotlib import cm                   # imports the colormap function from matplotlib
import cartopy
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

#------------------------------------------------------------------------------
# DEFINE THE DATASET

# CAMMPCAN 2017-18
V1_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V1_17_Data.csv',header=0,encoding = 'unicode_escape')
V2_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V2_17_Data.csv',header=0,encoding = 'unicode_escape')
V3_17 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V3_17_Data.csv',header=0,encoding = 'unicode_escape')

# CAMMPCAN 2018-19
V1_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V1_18_Data.csv',header=0,encoding = 'unicode_escape')
V2_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V2_18_Data.csv',header=0,encoding = 'unicode_escape')
V3_18 = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/V3_18_Data.csv',header=0,encoding = 'unicode_escape')

# SIPEXII 2012
SIPEXII = pd.read_csv('/Users/ncp532/Documents/Data/V1_17_Apriori/SIPEXII_Data.csv',header=0,encoding = 'unicode_escape')

#------------------------------------------------------------------------------
# FIX UP THE BrO DATASET

# Set the date
V1_17['DateTime'] = pd.to_datetime(V1_17['DateTime'])
V2_17['DateTime'] = pd.to_datetime(V2_17['DateTime'])
V3_17['DateTime'] = pd.to_datetime(V3_17['DateTime'])

V1_18['DateTime'] = pd.to_datetime(V1_18['DateTime'])
V2_18['DateTime'] = pd.to_datetime(V2_18['DateTime'])
V3_18['DateTime'] = pd.to_datetime(V3_18['DateTime'])

SIPEXII['DateTime'] = pd.to_datetime(SIPEXII['DateTime'])

#-------------------------
# Filter for Midday hours only

# V1_17 Davis (07:00 to 18:00)
start_time = '07:00:00'
end_time   = '18:00:00'
Midday     = (V1_17['Time'] >= start_time) & (V1_17['Time'] < end_time)
V1_17_MM   = V1_17[Midday]

# V2_17 Casey (08:00 to 16:00)
start_time = '08:00:00'
end_time   = '16:00:00'
Midday     = (V2_17['Time'] >= start_time) & (V2_17['Time'] < end_time)
V2_17_MM   = V2_17[Midday]

# V3_17 Mawson (08:00 to 18:00)
start_time = '08:00:00'
end_time   = '18:00:00'
Midday     = (V3_17['Time'] >= start_time) & (V3_17['Time'] < end_time)
V3_17_MM   = V3_17[Midday]

# V1_18 Davis (07:00 to 18:00)
start_time = '07:00:00'
end_time   = '18:00:00'
Midday     = (V1_18['Time'] >= start_time) & (V1_18['Time'] < end_time)
V1_18_MM   = V1_18[Midday]

# V2_18 Casey (08:00 to 16:00)
start_time = '08:00:00'
end_time   = '16:00:00'
Midday     = (V2_18['Time'] >= start_time) & (V2_18['Time'] < end_time)
V2_18_MM   = V2_18[Midday]

# V3_18 Mawson (08:00 to 18:00)
start_time = '08:00:00'
end_time   = '18:00:00'
Midday     = (V3_18['Time'] >= start_time) & (V3_18['Time'] < end_time)
V3_18_MM   = V3_18[Midday]

# SIPEXII (07:00 to 18:00)
start_time = '07:00:00'
end_time   = '18:00:00'
Midday     = (SIPEXII['Time'] >= start_time) & (SIPEXII['Time'] < end_time)
SIPEXII_MM = SIPEXII[Midday]

#-------------------------
# Filter dataframe for when filter is less than 60%

V1_17F = (V1_17_MM['Filter'] < 0.6)
V1_17T = V1_17_MM[V1_17F]

V2_17F = (V2_17_MM['Filter'] < 0.6)
V2_17T = V2_17_MM[V2_17F]

V3_17F = (V3_17_MM['Filter'] < 0.6)
V3_17T = V3_17_MM[V3_17F]

V1_18F = (V1_18_MM['Filter'] < 0.6)
V1_18T = V1_18_MM[V1_18F]

V2_18F = (V2_18_MM['Filter'] < 0.6)
V2_18T = V2_18_MM[V2_18F]

V3_18F = (V3_18_MM['Filter'] < 0.6)
V3_18T = V3_18_MM[V3_18F]

SIPEXIIF = (SIPEXII_MM['Filter'] < 0.6)
SIPEXIIT = SIPEXII_MM[SIPEXIIF]

#------------------------------------------------------------------------------
# DEFINE THE VARIABLES

# BrO
BrO_V1_17 = np.array(V1_17T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_V2_17 = np.array(V2_17T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_V3_17 = np.array(V3_17T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_V1_18 = np.array(V1_18T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_V2_18 = np.array(V2_18T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_V3_18 = np.array(V3_18T['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv
BrO_SIPEXII = np.array(SIPEXIIT['surf_vmr(ppmv)']) * 1e6 # convert from ppmv to pptv

# Latitude
Lat_V1_17 = np.array(V1_17T['LATITUDE'])
Lat_V2_17 = np.array(V2_17T['LATITUDE'])
Lat_V3_17 = np.array(V3_17T['LATITUDE'])
Lat_V1_18 = np.array(V1_18T['LATITUDE'])
Lat_V2_18 = np.array(V2_18T['LATITUDE'])
Lat_V3_18 = np.array(V3_18T['LATITUDE'])
Lat_SIPEXII = np.array(SIPEXIIT['LATITUDE'])

# Longitude
Long_V1_17 = np.array(V1_17T['LONGITUDE'])
Long_V2_17 = np.array(V2_17T['LONGITUDE'])
Long_V3_17 = np.array(V3_17T['LONGITUDE'])
Long_V1_18 = np.array(V1_18T['LONGITUDE'])
Long_V2_18 = np.array(V2_18T['LONGITUDE'])
Long_V3_18 = np.array(V3_18T['LONGITUDE'])
Long_SIPEXII = np.array(SIPEXIIT['LONGITUDE'])

#------------------------------------------------------------------------------
# START PLOTTING THE MAP

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

# We will view the map on a polar stereo coordinate system
projection = ccrs.SouthPolarStereo(central_longitude=0)

plt.close()
#fig = plt.figure(figsize=(6, 3))

fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
gs = gridspec.GridSpec(1,3,
                       width_ratios=[3, 3, 3],
                       height_ratios=[3])
                       #height_ratios=[3, 3, 3])

#--------------------------------------------
# GRAPH 1
ax1 = plt.subplot(gs[0], projection=projection) # (vertical no, horizontal no, graph no)

ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax1.coastlines()
gl = ax1.gridlines(color='gray',alpha=0.5,)
#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
#gl.ylocator = mticker.FixedLocator([-90,-45,0,45,90])
ax1.set_extent([0, 180, -45, -90], crs=ccrs.PlateCarree())

#---------------------------
# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
MacIsl_lon, MacIsl_lat = -54.2959, 158.5609
Hobart_lon, Hobart_lat = -42.8821, 147.3272

# Plot the station markers
ax1.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax1.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax1.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax1.plot(MacIsl_lat, MacIsl_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax1.plot(Hobart_lat, Hobart_lon, transform=data_crs, color='k', marker='*', markersize=10)

# Plot the marker labels
ax1.text(Davis_lat + 2, Davis_lon - 1, 'Davis', transform=data_crs, horizontalalignment='right')
ax1.text(Mawson_lat + 2, Mawson_lon - 1, 'Mawson', transform=data_crs, horizontalalignment='right')
ax1.text(Casey_lat + 2, Casey_lon - 1, 'Casey', transform=data_crs, horizontalalignment='right')
ax1.text(MacIsl_lat - 2, MacIsl_lon - 2, 'Macquarie\nIsland', transform=data_crs, horizontalalignment='right')
ax1.text(Hobart_lat + 3, Hobart_lon +2, 'Hobart', transform=data_crs, horizontalalignment='right')

# Scale marker size by BrO concentration (lower = smaller, higher = larger) 
s1 = BrO_SIPEXII * BrO_SIPEXII

#---------------------------
# PLOT THE SHIP TRACKS
cs = ax1.scatter(Long_SIPEXII, Lat_SIPEXII, transform=data_crs, c=BrO_SIPEXII, s=s1, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")

ax1.set_aspect('auto')

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("SIPEXII (2012)", fontweight = "bold", fontsize = 15, pad = 10, color = 'green')

#--------------------------------------------
# GRAPH 2
ax2 = plt.subplot(gs[1], projection=projection) # (vertical no, horizontal no, graph no)

ax2.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax2.coastlines()
gl = ax2.gridlines(color='gray',alpha=0.5,)
#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
#gl.ylocator = mticker.FixedLocator([-90,-45,0,45,90])
ax2.set_extent([0, 180, -45, -90], crs=ccrs.PlateCarree())

#---------------------------
# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
MacIsl_lon, MacIsl_lat = -54.2959, 158.5609
Hobart_lon, Hobart_lat = -42.8821, 147.3272

# Plot the station markers
ax2.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax2.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax2.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax2.plot(MacIsl_lat, MacIsl_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax2.plot(Hobart_lat, Hobart_lon, transform=data_crs, color='k', marker='*', markersize=10)

# Plot the marker labels
ax2.text(Davis_lat + 2, Davis_lon - 1, 'Davis', transform=data_crs, horizontalalignment='right')
ax2.text(Mawson_lat + 2, Mawson_lon - 1, 'Mawson', transform=data_crs, horizontalalignment='right')
ax2.text(Casey_lat + 2, Casey_lon - 1, 'Casey', transform=data_crs, horizontalalignment='right')
ax2.text(MacIsl_lat - 2, MacIsl_lon - 2, 'Macquarie\nIsland', transform=data_crs, horizontalalignment='right')
ax2.text(Hobart_lat + 3, Hobart_lon + 2, 'Hobart', transform=data_crs, horizontalalignment='right')

# Scale marker size by BrO concentration (lower = smaller, higher = larger) 
s1 = BrO_V1_17 * BrO_V1_17
s2 = BrO_V2_17 * BrO_V2_17
s3 = BrO_V3_17 * BrO_V3_17

#---------------------------
# PLOT THE SHIP TRACKS

ax2.scatter(Long_V1_17, Lat_V1_17, transform=data_crs, c=BrO_V1_17, s=s1, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")
ax2.scatter(Long_V2_17, Lat_V2_17, transform=data_crs, c=BrO_V2_17, s=s2, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")
ax2.scatter(Long_V3_17, Lat_V3_17, transform=data_crs, c=BrO_V3_17, s=s3, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")

ax2.set_aspect('auto')

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("CAMMPCAN (2017-18)", fontweight = "bold", fontsize = 15, pad = 10, color = 'blue')

#--------------------------------------------
# GRAPH 3
ax3 = plt.subplot(gs[2], projection=projection) # (vertical no, horizontal no, graph no)

ax3.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax3.coastlines()
gl = ax3.gridlines(color='gray',alpha=0.5,)
#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
#gl.ylocator = mticker.FixedLocator([-90,-45,0,45,90])
ax3.set_extent([0, 180, -45, -90], crs=ccrs.PlateCarree())

#---------------------------
# Station location (Lon, Lat)
Davis_lon, Davis_lat = -68.5766, 77.9674
Mawson_lon, Mawson_lat = -67.6027, 62.8738
Casey_lon, Casey_lat = -66.2818, 110.5276
MacIsl_lon, MacIsl_lat = -54.2959, 158.5609
Hobart_lon, Hobart_lat = -42.8821, 147.3272

# Plot the station markers
ax3.plot(Davis_lat,  Davis_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax3.plot(Mawson_lat, Mawson_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax3.plot(Casey_lat,  Casey_lon,  transform=data_crs, color='k', marker='*', markersize=10)
ax3.plot(MacIsl_lat, MacIsl_lon, transform=data_crs, color='k', marker='*', markersize=10)
ax3.plot(Hobart_lat, Hobart_lon, transform=data_crs, color='k', marker='*', markersize=10)

# Plot the marker labels
ax3.text(Davis_lat + 2,  Davis_lon - 1,  'Davis',  transform=data_crs, horizontalalignment='right')
ax3.text(Mawson_lat + 2, Mawson_lon - 1, 'Mawson', transform=data_crs, horizontalalignment='right')
ax3.text(Casey_lat + 2,  Casey_lon - 1,  'Casey',  transform=data_crs, horizontalalignment='right')
ax3.text(MacIsl_lat - 2, MacIsl_lon - 2, 'Macquarie\nIsland', transform=data_crs, horizontalalignment='right')
ax3.text(Hobart_lat + 3, Hobart_lon + 2, 'Hobart', transform=data_crs, horizontalalignment='right')

# Scale marker size by BrO concentration (lower = smaller, higher = larger) 
s1 = BrO_V1_18 * BrO_V1_18
s2 = BrO_V2_18 * BrO_V2_18
s3 = BrO_V3_18 * BrO_V3_18

#---------------------------
# PLOT THE SHIP TRACKS

ax3.scatter(Long_V1_18, Lat_V1_18, transform=data_crs, c=BrO_V1_18, s=s1, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")
ax3.scatter(Long_V2_18, Lat_V2_18, transform=data_crs, c=BrO_V2_18, s=s2, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")
ax3.scatter(Long_V3_18, Lat_V3_18, transform=data_crs, c=BrO_V3_18, s=s3, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")

# Plot the colorbar
cb = fig.colorbar(cs,ticks=[2,2.5,3,3.5,4,4.5,5,5.5,6])
ax3.set_aspect('auto')
cb.set_label('BrO VMR (pptv)', rotation=90,fontsize=15)

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("CAMMPCAN (2018-19)", fontweight = "bold", fontsize = 15, pad = 10, color = 'red')
