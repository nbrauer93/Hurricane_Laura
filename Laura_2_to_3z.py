#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 21:12:58 2021

@author: noahbrauer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 20:59:35 2021

@author: noahbrauer
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from netCDF4 import Dataset, num2date, MFDataset
import glob







import pyart

#Read in file

file = 'ncswp_SR1_20200827_025927.829_v155_s000_160.0_RHI_.nc'


#Import file

nc = Dataset(file, 'r')


#Extract attributes in file

zh = nc.variables['DBZ'][:]
zdr = nc.variables['ZDR'][:]
kdp = nc.variables['KDP'][:]
rho_hv = nc.variables['RHOHV'][:]


###Account for bird-bath scan; Zdr calibration

zdr_calibration = nc.variables['zdr_correction_db'][:]


elevation = nc.variables['Elevation'][:]


#Now do this for multiple files:

files = glob.glob('*.nc')

zh_all = np.ones((len(files), zh.shape[0], zh.shape[1]))*np.nan
zdr_all = np.ones((len(files), zdr.shape[0], zdr.shape[1]))*np.nan
rhohv_all = np.ones((len(files), rho_hv.shape[0], rho_hv.shape[1]))*np.nan


for n,j in enumerate(files):
    
    data = Dataset(j, 'r')
    print(j)
    
    zh_all[n,:,:] = data.variables['DBZ'][:]
    zdr_all[n,:,:] = data.variables['ZDR'][:]
    rhohv_all[n,:,:] = data.variables['RHOHV'][:]


#Compute mean values from 8/27 0117 UTC - 0435 UTC
    
zh_mean = np.nanmean(zh_all, axis = 0)
zdr_mean = np.nanmean(zdr_all, axis = 0)
rhohv_mean = np.nanmean(rhohv_all, axis = 0)    


#NaN out bad data


zh_mean[zh_mean<=0] = np.nan
zdr_mean[zdr_mean<=-4] = np.nan
rhohv_mean[rhohv_mean<=0.8] = np.nan

#Now calculate beam height (height) and distance (x) from the radar using Doviak and Zrnic:


deg_to_rad = np.pi/180. #Convert from degrees to radians

range_gate = np.zeros_like(zh[0,:])
range_gate = np.arange(0,len(range_gate))
range_gate = nc.variables['Range_to_First_Cell'].getValue() + range_gate*nc.variables['Cell_Spacing'].getValue()
a_e = (8494.66667)*1000 


height = np.ones(zh.shape)*np.nan
s = np.ones(zh.shape)*np.nan  

for i in range(0,len(elevation)):

    a = range_gate**2 + a_e**2 + (2.0*range_gate*a_e*np.sin(elevation[i]*deg_to_rad))
    height[i,:] =((a**0.5)-a_e)+2
    s[i,:] = a_e*np.arcsin((range_gate*np.cos(elevation[i]*deg_to_rad))/(a_e+height[i,:]))



#Now let's plot les data   
    #%%
    

label_size = 24
tick_label_declutter = 200



    
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, zh, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(5,70,step = 2.5))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(0,175,25)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)

plt.title(r'SMART-R1 $160^{o}$ Azimuth $Z_{H}$ 8/27 0434 UTC', size = label_size)
plt.show()


#%%
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, zdr-2.5, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0,6,step = 0.25))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dB]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(0,175,25)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)



plt.title(r'SMART-R1 $160^{o}$ Azimuth $Z_{DR}$ 8/27 0434 UTC', size = label_size)
plt.show()

'''

fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, kdp, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0,4,step = 0.2))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[deg/km]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(0,175,25)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)



plt.title(r'SMART-R1 $160^{o}$ Azimuth $K_{DP}$ 8/27 0434 UTC', size = label_size)
plt.show()


'''

    
fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, rho_hv, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0.9,1.05,step = 0.01))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[$\rho_{hv}$]',size = label_size)
plt.xlabel('Distance from Radar (km)', size = label_size)
plt.ylabel('Altitude (km)', size = label_size)

xlabels = np.arange(0,175,25)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = label_size)
plt.yticks(ylabels, size = label_size)



plt.title(r'SMART-R1 $160^{o}$ Azimuth $\rho_{hv}$ 8/27 0434 UTC', size = label_size)
plt.show()





#%%

#Plot time composites


tick_label_size = 32

fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, zh_mean, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(5,70,step = 2.5))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = '[dBZ]',size = tick_label_size)
plt.xlabel('Distance from Radar (km)', size = tick_label_size)
plt.ylabel('Altitude (km)', size = tick_label_size)

xlabels = np.arange(0,175,10)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = tick_label_size)
plt.yticks(ylabels, size = tick_label_size)
plt.xlim(0,70)
plt.title(r'SR1-P $160^{o}$ Azimuth Mean $Z_{H}$ 8/27 0211-0259 UTC', size = label_size)
plt.show()


fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, zdr_mean-1.7, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(-1,3,step = 0.1))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = '[dB]',size = tick_label_size)
plt.xlabel('Distance from Radar (km)', size = tick_label_size)
plt.ylabel('Altitude (km)', size = tick_label_size)

xlabels = np.arange(0,175,10)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = tick_label_size)
plt.yticks(ylabels, size = tick_label_size)
plt.xlim(0,70)
plt.title(r'SR1-P $160^{o}$ Azimuth Mean $Z_{DR}$ 8/27 0211-0259 UTC', size = label_size)
plt.show()


fig,ax = plt.subplots(figsize=(14,14))
plt.contourf(s/1000,height/1000, rhohv_mean, extend = 'both', cmap = 'pyart_NWSRef', levels = np.arange(0.9,1,step = 0.01))
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = r'[$\rho_{hv}$]',size = tick_label_size)
plt.xlabel('Distance from Radar (km)', size = tick_label_size)
plt.ylabel('Altitude (km)', size = tick_label_size)

xlabels = np.arange(0,175,10)
ylabels = np.arange(0,11,1)

plt.ylim(0,10)

plt.xticks(xlabels, size = tick_label_size)
plt.yticks(ylabels, size = tick_label_size)
plt.xlim(0,70)
plt.title(r'SR1-P $160^{o}$ Azimuth Mean $\rho_{hv}$ 8/27 0211-0259 UTC', size = label_size)
plt.show()




 



