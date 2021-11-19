#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 17:30:14 2020

@author: noahbrauer
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from netCDF4 import Dataset, num2date
import os


import pyart
from mpl_toolkits.basemap import Basemap
import math


file = 'KLCH_N0R_20200827_055300.nc'

nc = Dataset(file, 'r')

lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

lat2,lon2 = np.meshgrid(lat,lon)
time = nc.variables['time'][:]
z = nc.variables['bref'][:]



plt.figure(figsize=(14,14))

min_value = -5
max_value = 75
value_interval = 1
title_font_size = 26
label_size = 18
text_size = 12
cbar_label_size = 22

cmin = min_value; cmax = max_value; cint = value_interval; clevs = np.round(np.arange(cmin,cmax,cint),2)
xlim = np.array([-95,-91]); ylim = np.array([29,32])
   
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(linewidth = 2); m.drawcountries()

cs = m.contourf(lon2,lat2,z[0,:,:].T, clevs, cmap = 'pyart_NWSRef', extend = 'both')

m.drawcounties(linewidth = 0.25)

cbar = plt.colorbar(fraction = 0.035)
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = cbar_label_size)
plt.title(r'KLCH 8/27 0553 UTC $Z_{H}$', size = title_font_size)


label = 'PIPS 1A'
x2star,y2star = m(-93.518746,30.224910)
m.plot(x2star,y2star,'ro',markersize=12, color = 'k')
plt.text(x2star-0.14,y2star-0.1,label, weight = 'bold', size = text_size)

label = 'PIPS 1B'
x2star,y2star = m(-93.359979,30.266642)
m.plot(x2star,y2star,'ro',markersize=12, color = 'k')
plt.text(x2star+0.06,y2star-0.02,label, weight = 'bold', size = text_size)

label = 'PIPS 2A'
x2star,y2star = m(-92.928744,30.245168)
m.plot(x2star,y2star,'ro',markersize=12, color = 'k')
plt.text(x2star+0.03,y2star+0.05,label, weight = 'bold', size = text_size)

label = 'PIPS 2B'
x2star,y2star = m(-92.739657,30.246693)
m.plot(x2star,y2star,'ro',markersize=12, color = 'k')
plt.text(x2star+0.05,y2star-0.05,label, weight = 'bold', size = text_size)

label = 'KLCH'
x2star,y2star = m(-93.2211,30.1230)
m.plot(x2star,y2star,'r*',markersize=12, color = 'b')
plt.text(x2star+0.03,y2star+0.03,label, weight = 'bold', size = text_size)


label = 'SR1-P'
x2star,y2star = m(-93.22,30.36)
m.plot(x2star,y2star,'r*',markersize=12, color = 'b')
plt.text(x2star+0.03,y2star+0.03,label, weight = 'bold', size = text_size)

x2star,y2star = m(-93.28,29.7)
m.plot(x2star,y2star, "ro", markersize = 12, color = 'fuchsia')


label = '100 km'
circle = Circle(xy = m(-93.28,29.7), radius = 0.9009, fill = False, linewidth = 2)
plt.gca().add_patch(circle)
plt.text(-93.35,30.615, label, weight = 'bold', size = 14)

label = '250 km'
circle = Circle(xy = m(-93.28,29.7), radius = 2.252, fill = False, linewidth = 2)
plt.gca().add_patch(circle)
plt.text(-93.35,31.85, label, weight = 'bold', size = 14)


plt.show()