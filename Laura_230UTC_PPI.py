#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 10:41:41 2021

@author: noahbrauer
"""

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

import pyart
import cartopy
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader



file = 'ncswp_SR1_20200827_023021.671_v005_s004_004.0_SUR_.nc'

nc = Dataset(file, 'r')

azimuth = nc.variables['Azimuth'][:]
z = nc.variables['DBZ'][:]
zdr = nc.variables['ZDR'][:]-2.5
radar_range = np.arange(z.shape[1])*nc.variables['Cell_Spacing'][:]+ nc.variables['Range_to_First_Cell'][:]
RadarLatitude = nc.variables['Latitude'][:]
RadarLongitude = nc.variables['Longitude'][:]


#Now compute to cartesian coordinates

x = radar_range*np.sin(np.deg2rad(azimuth))[:,None]
y = radar_range*np.cos(np.deg2rad(azimuth))[:,None]


print(np.nanmax(zdr))
print(np.nanmin(zdr))


#%%

#Now we plot

radar_data_refl = np.ma.array(z, mask = np.isnan(z))
proj = cartopy.crs.LambertConformal(central_longitude = RadarLongitude, central_latitude = RadarLatitude)
    
state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())    



fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = 0.; cmax = 60.; cint = 5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)



mesh = ax.pcolormesh(x,y,radar_data_refl, cmap = 'pyart_NWSRef', vmin = 0, vmax = 60)

ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='black')

distance_in_degrees = 1.2

ax.set_extent([RadarLongitude - distance_in_degrees, RadarLongitude + distance_in_degrees, RadarLatitude - distance_in_degrees, RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(18)
    
  
    
plt.title(r'SMART-R1 $Z_{H}$ 8/27 0230 UTC', size = 24) 


#Convert x,y distance in km to degrees

def dist_to_deg(distance):
    
    degree_distance = distance/111.
    
    return degree_distance


x_sect_distance = dist_to_deg(70.)

#Convert lat and longitudinal distance to degrees 

lat_dist = np.cos(200)*x_sect_distance*-1
lon_dist = np.sin(200)*x_sect_distance


#Compute end location of cross-section location


lat_end = RadarLatitude+lat_dist
lon_end = RadarLongitude-lon_dist


#Now add SR-1 RHI location at 160 degree azimuth

start = plt.plot([RadarLongitude,lon_end], [RadarLatitude,lat_end], marker = 'o', color = 'k', markersize = 12, linewidth = 4,transform=ccrs.Geodetic(),)





#end = plt.plot([lon_end,lat_end], marker = 'o', color = 'k', markersize = 12)





plt.show()



#%%
#Now for Zdr

radar_data_zdr = np.ma.array(zdr, mask = np.isnan(zdr))
proj = cartopy.crs.LambertConformal(central_longitude = RadarLongitude, central_latitude = RadarLatitude)
    
state_borders = cartopy.feature.NaturalEarthFeature(category='cultural', name = 'admin_1_states_provinces_lakes', scale = '50m', facecolor = 'none')

reader = shpreader.Reader('countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())    



fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1,projection = proj)
cmin = -1; cmax = 2.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)



mesh = ax.pcolormesh(x,y,radar_data_zdr, cmap = 'pyart_NWSRef', vmin = -1, vmax = 2)

ax.add_feature(state_borders, edgecolor = 'black', linewidth = 1.5, zorder = 2)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='black')

distance_in_degrees = 1.2

ax.set_extent([RadarLongitude - distance_in_degrees, RadarLongitude + distance_in_degrees, RadarLatitude - distance_in_degrees, RadarLatitude + distance_in_degrees]) 

cbar = plt.colorbar(mesh, shrink = 0.9)
cbar.ax.set_ylabel('[dB]',name='Calibri',size=22)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)

    cbar.set_ticks(clevs[::2])
    cbar.set_ticklabels(cticks[::2])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(18)
    
  
    
plt.title(r'SMART-R1 $Z_{DR}$ 8/27 0230 UTC', size = 24) 


#Convert x,y distance in km to degrees

def dist_to_deg(distance):
    
    degree_distance = distance/111.
    
    return degree_distance


x_sect_distance = dist_to_deg(70.)

#Convert lat and longitudinal distance to degrees 

lat_dist = np.cos(200(np.pi/180))*x_sect_distance
lon_dist = np.sin(200(np.pi/180))*x_sect_distance


#Compute end location of cross-section location


lat_end = RadarLatitude-lat_dist
lon_end = RadarLongitude-lon_dist


#Now add SR-1 RHI location at 160 degree azimuth

start = plt.plot([RadarLongitude,lon_end], [RadarLatitude,lat_end], marker = 'o', color = 'k', markersize = 12, linewidth = 4,transform=ccrs.Geodetic(),)





#end = plt.plot([lon_end,lat_end], marker = 'o', color = 'k', markersize = 12)





plt.show()



