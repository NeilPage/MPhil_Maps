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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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
# START PLOTTING THE MAP (BrO Maximum)

# The data are defined in lat/lon coordinate system, so PlateCarree()
# is the appropriate coordinate system:
data_crs = ccrs.PlateCarree()

# We will view the map on a polar stereo coordinate system
projection = ccrs.SouthPolarStereo(central_longitude=0)

plt.close()
#fig = plt.figure(figsize=(6, 3))

fig1 = plt.figure()
plt.subplots_adjust(hspace=0.2)
gs = gridspec.GridSpec(nrows=2,
                       ncols=2,
                       width_ratios=[1,1],
                       height_ratios=[1,1])
                       #height_ratios=[3, 3, 3])

#--------------------------------------------
# GRAPH 1
ax1 = plt.subplot(gs[1,0], projection=projection) # (vertical no, horizontal no, graph no)

#ax.set_global()
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax1.coastlines()
gl = ax1.gridlines(color='gray',alpha=0.5,)
#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
#gl.ylocator = mticker.FixedLocator([-90,-45,0,45,90])
ax1.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())

#---------------------------
# Station location (Lon, Lat)
Davis_lon,        Davis_lat        = -68.5766, 77.9674
Mawson_lon,       Mawson_lat       = -67.6027, 62.8738
Casey_lon,        Casey_lat        = -66.2818, 110.5276
#MacIsl_lon,       MacIsl_lat       = -54.2959, 158.5609
#Hobart_lon,       Hobart_lat       = -42.8821, 147.3272
SIPEXII_lon,      SIPEXII_lat      = -61.5205, 121.1855
Halley_lon,       Halley_lat       = -75.583,  -26.567
Scott_lon,        Scott_lat        = -77.848,  166.760
McMurdoSound_lon, McMurdoSound_lat = -77.820,  166.580
Marambio_lon,     Marambio_lat     = -64.241,  -56.627
BelgranoII_lon,   BelgranoII_lat   = -77.874,  -34.628
Neumayer_lon,     Neumayer_lat     = -70.645,  -8.264
DumontDUrv_lon,   DumontDUrv_lat   = -66.663,  140.002
DomeC_lon,        DomeC_lat        = -75.0559, 123.1956
ArrivalHei_lon,   ArrivalHei_lat   = -77.83,   166.66

# Max BrO Observation
Davis_max        = 14.8
Mawson_max       = 14.3
Casey_max        = 21.5
SIPEXII_max      = 15.0
Halley_max       = 20.2
Scott_max        = 10.9
McMurdoSound_max = 14.4
Marambio_max     = 26.0
BelgranoII_max   = 8.1
Neumayer_max     = 13.0
DumontDUrv_max   = 2.0
DomeC_max        = 2.5
ArrivalHei_max   = 30.0

# Plot the station markers
#ax1.plot(Davis_lat,        Davis_lon,        transform=data_crs, color='k', marker='o')
#ax1.plot(Mawson_lat,       Mawson_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(Casey_lat,        Casey_lon,        transform=data_crs, color='k', marker='o')
##ax1.plot(MacIsl_lat,       MacIsl_lon,       transform=data_crs, color='k', marker='o')
##ax1.plot(Hobart_lat,       Hobart_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(SIPEXII_lat,      SIPEXII_lon,      transform=data_crs, color='k', marker='o')
#ax1.plot(Halley_lat,       Halley_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(Scott_lat,        Scott_lon,        transform=data_crs, color='k', marker='o')
##ax1.plot(McMurdoSound_lat, McMurdoSound_lon, transform=data_crs, color='k', marker='o')
#ax1.plot(Marambio_lat,     Marambio_lon,     transform=data_crs, color='k', marker='o')
#ax1.plot(BelgranoII_lat,   BelgranoII_lon,   transform=data_crs, color='k', marker='o')
#ax1.plot(Neumayer_lat,     Neumayer_lon,     transform=data_crs, color='k', marker='o')
#ax1.plot(DumontDUrv_lat,  DumontDUrv_lon,   transform=data_crs, color='k', marker='o')

# Main Map
mylons = [Davis_lon, Mawson_lon, Casey_lon, SIPEXII_lon]
mylats = [Davis_lat, Mawson_lat, Casey_lat, SIPEXII_lat]
myobs  = [Davis_max, Mawson_max, Casey_max, SIPEXII_max]
cs = ax1.scatter(mylats, mylons, c=myobs, transform=data_crs, vmin=0, vmax = 30, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=200, zorder=2, marker='*')

lons = [Halley_lon, Scott_lon, Marambio_lon, BelgranoII_lon, Neumayer_lon, DumontDUrv_lon, McMurdoSound_lon, DomeC_lon, ArrivalHei_lon]
lats = [Halley_lat, Scott_lat, Marambio_lat, BelgranoII_lat, Neumayer_lat, DumontDUrv_lat, McMurdoSound_lat, DomeC_lat, ArrivalHei_lat]
obs  = [Halley_max, Scott_max, Marambio_max, BelgranoII_max, Neumayer_max, DumontDUrv_max, McMurdoSound_max, DomeC_max, ArrivalHei_max]
ax1.scatter(lats, lons, c=obs, transform=data_crs, vmin=0, vmax = 30, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=100, zorder=2)

# Inset Map
lons2 = [Scott_lon, McMurdoSound_lon, ArrivalHei_lon]
lats2 = [Scott_lat, McMurdoSound_lat, ArrivalHei_lat]
obs2  = [Scott_max, McMurdoSound_max, ArrivalHei_max]

# Plot the marker labels
ax1.text(Davis_lat +2,        Davis_lon -2,         'Davis',             transform=data_crs, horizontalalignment='right')
ax1.text(Mawson_lat +2,       Mawson_lon -2,        'Mawson',            transform=data_crs, horizontalalignment='right')
ax1.text(Casey_lat +2,        Casey_lon -1.5,         'Casey',             transform=data_crs, horizontalalignment='right')
#ax1.text(MacIsl_lat -2,       MacIsl_lon -2,        'Macquarie\nIsland', transform=data_crs, horizontalalignment='right')
#ax1.text(Hobart_lat +3,       Hobart_lon +2,        'Hobart',            transform=data_crs, horizontalalignment='right')
ax1.text(SIPEXII_lat +5,      SIPEXII_lon +2,       'SIPEXII',           transform=data_crs, horizontalalignment='right')
ax1.text(Halley_lat +3,       Halley_lon +1,        'Halley V',          transform=data_crs, horizontalalignment='right')
#ax1.text(Scott_lat +2,        Scott_lon -1,         'Scott',             transform=data_crs, horizontalalignment='right')
#ax1.text(McMurdoSound_lat +2, McMurdoSound_lon -1,  'McMurdo Sound',     transform=data_crs, horizontalalignment='right')
ax1.text(Marambio_lat + 4,    Marambio_lon,         'Marambio',          transform=data_crs, horizontalalignment='right')
ax1.text(BelgranoII_lat -4,   BelgranoII_lon,       'Belgrano II',       transform=data_crs, horizontalalignment='right')
ax1.text(Neumayer_lat -1,     Neumayer_lon +1,      'Neumayer',          transform=data_crs, horizontalalignment='right')
ax1.text(DumontDUrv_lat +2,   DumontDUrv_lon +4,    "Dumont d'Urville",  transform=data_crs, horizontalalignment='right')
ax1.text(DomeC_lat +5,        DomeC_lon +4,         'Dome-C',            transform=data_crs, horizontalalignment='right')
#ax1.text(ArrivalHei_lat +3,   ArrivalHei_lon +3,    'Arrival Heights',   transform=data_crs, horizontalalignment='right')

#---------------------------
# PLOT THE SHIP TRACKS
#cs = ax1.scatter(Long_SIPEXII, Lat_SIPEXII, transform=data_crs, c=BrO_SIPEXII, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")

ax1.set_aspect('auto')

# Plot the colorbar
cb = fig1.colorbar(cs,ticks=[0,3,6,9,12,15,18,21,24,27,30])
ax1.set_aspect('auto')
cb.set_label('BrO VMR (pptv)', rotation=90,fontsize=15)

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("Maximum BrO", fontweight = "bold", fontsize = 15, pad = 10, color = 'black')

#---------------------------
# PLOT THE MAP INSET
#ax2 = ax1.twinx()

# Plot an inset over the main map
axins = inset_axes(ax1, width="45%", height="30%", loc=3, axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(map_projection=projection))#cartopy.crs.PlateCarree()))
axins.add_feature(cartopy.feature.COASTLINE)
axins.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
axins.coastlines()
#axins.gridlines(color='gray',alpha=0.5,)
#axins.set_extent([166.0, 167.5, -77.5, -78], crs=ccrs.PlateCarree())
axins.set_extent([165.0, 169.0, -77.5, -78.5], crs=ccrs.PlateCarree())
axins.scatter(lats2, lons2, c=obs2, transform=data_crs, vmin=0, vmax = 30, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=100, zorder=2)

axins.text(Scott_lat +0.5,        Scott_lon +0.02,        'Scott',             transform=data_crs, horizontalalignment='right', color = 'black')
axins.text(McMurdoSound_lat -2.5, McMurdoSound_lon +0.02, 'McMurdo Sound',     transform=data_crs, horizontalalignment='right', color = 'black')
axins.text(ArrivalHei_lat -0.5,   ArrivalHei_lon +0.2,   'Arrival Heights',   transform=data_crs, horizontalalignment='right', color = 'black')

mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# #------------------------------------------------------------------------------
# # START PLOTTING THE MAP (BrO Median)

# # The data are defined in lat/lon coordinate system, so PlateCarree()
# # is the appropriate coordinate system:
# data_crs = ccrs.PlateCarree()

# # We will view the map on a polar stereo coordinate system
# projection = ccrs.SouthPolarStereo(central_longitude=0)

# plt.close()
# #fig = plt.figure(figsize=(6, 3))

# fig1 = plt.figure()
# plt.subplots_adjust(hspace=0.5)
# gs = gridspec.GridSpec(nrows=1,
#                        ncols=1,
#                        width_ratios=[1],
#                        height_ratios=[1])
#                        #height_ratios=[3, 3, 3])

#--------------------------------------------
# GRAPH 1
ax1 = plt.subplot(gs[0,0], projection=projection) # (vertical no, horizontal no, graph no)

#ax.set_global()
ax1.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
ax1.coastlines()
gl = ax1.gridlines(color='gray',alpha=0.5,)
#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180])
#gl.ylocator = mticker.FixedLocator([-90,-45,0,45,90])
ax1.set_extent([-180, 180, -55, -90], crs=ccrs.PlateCarree())

#---------------------------
# Station location (Lon, Lat)
Davis_lon,        Davis_lat        = -68.5766, 77.9674
Mawson_lon,       Mawson_lat       = -67.6027, 62.8738
Casey_lon,        Casey_lat        = -66.2818, 110.5276
#MacIsl_lon,       MacIsl_lat       = -54.2959, 158.5609
#Hobart_lon,       Hobart_lat       = -42.8821, 147.3272
SIPEXII_lon,      SIPEXII_lat      = -61.5205, 121.1855
Halley_lon,       Halley_lat       = -75.583,  -26.567
Scott_lon,        Scott_lat        = -77.848,  166.760
McMurdoSound_lon, McMurdoSound_lat = -77.820,  166.580
Marambio_lon,     Marambio_lat     = -64.241,  -56.627
BelgranoII_lon,   BelgranoII_lat   = -77.874,  -34.628
Neumayer_lon,     Neumayer_lat     = -70.645,  -8.264
DumontDUrv_lon,   DumontDUrv_lat   = -66.663,  140.002
DomeC_lon,        DomeC_lat        = -75.0559, 123.1956
ArrivalHei_lon,   ArrivalHei_lat   = -77.83,   166.66

# Median BrO Observation
Davis_median        = 4.0
Mawson_median       = 3.4
Casey_median        = 3.8
SIPEXII_median      = 2.4
Halley_median       = 3.0
Scott_median        = 3.8
McMurdoSound_median = 6.35
Marambio_median     = 1.6
BelgranoII_median   = 1.6
DumontDUrv_median   = 2.0
DomeC_median        = 2.5

# Plot the station markers
#ax1.plot(Davis_lat,        Davis_lon,        transform=data_crs, color='k', marker='o')
#ax1.plot(Mawson_lat,       Mawson_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(Casey_lat,        Casey_lon,        transform=data_crs, color='k', marker='o')
##ax1.plot(MacIsl_lat,       MacIsl_lon,       transform=data_crs, color='k', marker='o')
##ax1.plot(Hobart_lat,       Hobart_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(SIPEXII_lat,      SIPEXII_lon,      transform=data_crs, color='k', marker='o')
#ax1.plot(Halley_lat,       Halley_lon,       transform=data_crs, color='k', marker='o')
#ax1.plot(Scott_lat,        Scott_lon,        transform=data_crs, color='k', marker='o')
##ax1.plot(McMurdoSound_lat, McMurdoSound_lon, transform=data_crs, color='k', marker='o')
#ax1.plot(Marambio_lat,     Marambio_lon,     transform=data_crs, color='k', marker='o')
#ax1.plot(BelgranoII_lat,   BelgranoII_lon,   transform=data_crs, color='k', marker='o')
#ax1.plot(Neumayer_lat,     Neumayer_lon,     transform=data_crs, color='k', marker='o')
#ax1.plot(DumontDUrv_lat,  DumontDUrv_lon,   transform=data_crs, color='k', marker='o')

# Main Map
mylons = [Davis_lon, Mawson_lon, Casey_lon, SIPEXII_lon]
mylats = [Davis_lat, Mawson_lat, Casey_lat, SIPEXII_lat]
myobs  = [Davis_median, Mawson_median, Casey_median, SIPEXII_median]
cs = ax1.scatter(mylats, mylons, c=myobs, transform=data_crs, vmin=0, vmax = 7, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=200, zorder=2, marker='*')

lons = [Halley_lon, Marambio_lon, BelgranoII_lon, DumontDUrv_lon, McMurdoSound_lon, DomeC_lon, Scott_lon]
lats = [Halley_lat, Marambio_lat, BelgranoII_lat, DumontDUrv_lat, McMurdoSound_lat, DomeC_lat, Scott_lat]
obs  = [Halley_median, Marambio_median, BelgranoII_median, DumontDUrv_median, McMurdoSound_median, DomeC_median, Scott_median,]
ax1.scatter(lats, lons, c=obs, transform=data_crs, vmin=0, vmax = 7, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=100, zorder=2)

# Inset Map
lons2 = [Scott_lon, McMurdoSound_lon]
lats2 = [Scott_lat, McMurdoSound_lat]
obs2  = [Scott_median, McMurdoSound_median]

# Plot the marker labels
ax1.text(Davis_lat +2,        Davis_lon -2,         'Davis',             transform=data_crs, horizontalalignment='right')
ax1.text(Mawson_lat +2,       Mawson_lon -2,        'Mawson',            transform=data_crs, horizontalalignment='right')
ax1.text(Casey_lat +2,        Casey_lon -2,         'Casey',             transform=data_crs, horizontalalignment='right')
#ax1.text(MacIsl_lat -2,       MacIsl_lon -2,        'Macquarie\nIsland', transform=data_crs, horizontalalignment='right')
#ax1.text(Hobart_lat +3,       Hobart_lon +2,        'Hobart',            transform=data_crs, horizontalalignment='right')
ax1.text(SIPEXII_lat +5,      SIPEXII_lon +2,       'SIPEXII',           transform=data_crs, horizontalalignment='right')
ax1.text(Halley_lat +3,       Halley_lon +1,        'Halley V',          transform=data_crs, horizontalalignment='right')
#ax1.text(Scott_lat +2,        Scott_lon -1,         'Scott',             transform=data_crs, horizontalalignment='right')
#ax1.text(McMurdoSound_lat +2, McMurdoSound_lon -1,  'McMurdo Sound',     transform=data_crs, horizontalalignment='right')
ax1.text(Marambio_lat + 4,    Marambio_lon,         'Marambio',          transform=data_crs, horizontalalignment='right')
ax1.text(BelgranoII_lat -5,   BelgranoII_lon,       'Belgrano II',       transform=data_crs, horizontalalignment='right')
#ax1.text(Neumayer_lat -1,     Neumayer_lon +1,      'Neumayer',          transform=data_crs, horizontalalignment='right')
ax1.text(DumontDUrv_lat +2,   DumontDUrv_lon +4,    "Dumont d'Urville",  transform=data_crs, horizontalalignment='right')
ax1.text(DomeC_lat +5,        DomeC_lon +4,         'Dome-C',            transform=data_crs, horizontalalignment='right')
#ax1.text(ArrivalHei_lat +3,   ArrivalHei_lon +3,    'Arrival Heights',   transform=data_crs, horizontalalignment='right')

#---------------------------
# PLOT THE SHIP TRACKS
#cs = ax1.scatter(Long_SIPEXII, Lat_SIPEXII, transform=data_crs, c=BrO_SIPEXII, vmin=2.0, vmax=6.0, cmap=cm.get_cmap('jet',16), label="SIPEXII")

ax1.set_aspect('auto')

# Plot the colorbar
cb = fig1.colorbar(cs,ticks=[0,1,2,3,4,5,6,7])
ax1.set_aspect('auto')
cb.set_label('BrO VMR (pptv)', rotation=90,fontsize=15)

# PLOT TITLE, AXIS LABEL & LEGEND TITLE
plt.title("Medium BrO", fontweight = "bold", fontsize = 15, pad = 10, color = 'black')

#---------------------------
# PLOT THE MAP INSET
#ax2 = ax1.twinx()

# Plot an inset over the main map
axins = inset_axes(ax1, width="45%", height="30%", loc=3, axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(map_projection=projection))#cartopy.crs.PlateCarree()))
axins.add_feature(cartopy.feature.COASTLINE)
axins.add_feature(cartopy.feature.LAND, zorder=1, edgecolor='k',color='grey')
axins.coastlines()
#axins.gridlines(color='gray',alpha=0.5,)
#axins.set_extent([166.0, 167.5, -77.5, -78], crs=ccrs.PlateCarree())
axins.set_extent([165.0, 169.0, -77.5, -78.5], crs=ccrs.PlateCarree())
axins.scatter(lats2, lons2, c=obs2, transform=data_crs, vmin=0, vmax = 7, cmap=cm.get_cmap('viridis',30), edgecolors='black', s=100, zorder=2)

axins.text(Scott_lat +0.5,        Scott_lon +0.02,        'Scott',             transform=data_crs, horizontalalignment='right', color = 'black')
axins.text(McMurdoSound_lat -2.5, McMurdoSound_lon +0.02, 'McMurdo Sound',     transform=data_crs, horizontalalignment='right', color = 'black')
#axins.text(ArrivalHei_lat -1.0,   ArrivalHei_lon +0.175,   'Arrival Heights',   transform=data_crs, horizontalalignment='right', color = 'black')

mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

