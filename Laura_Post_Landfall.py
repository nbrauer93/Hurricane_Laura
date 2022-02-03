import h5py 
import numpy as np
import matplotlib.pyplot as plt
#import cartopy
from pyproj import Proj
import scipy
from scipy import stats





from mpl_toolkits.basemap import Basemap
import pyart


file = '2A.GPM.DPR.V8-20180723.20200827-S114415-E131648.036913.V06A.HDF5'


#Read in file, extract lat, lons, PSD

DPR = h5py.File(file, 'r')

#nscan is the nummber of scans in each granule 
#nray is the number of angle bins in each  scan; think of these as footprint scans (5 km in diameter for each footprint)
#nbin is the number of range bins in each ray (angle bins)
#nDSD: Parameters are N0 (number concentration) and D0 (mean drop diameter)


lat = DPR['NS']['Latitude'][:,:]    
lon = DPR['NS']['Longitude'][:,:]
z = DPR['NS']['SLV']['zFactorCorrectedNearSurface'][:] #nscan x nray (7934,49)
dsd = DPR['NS']['SLV']['paramDSD'][:] #nscan x nray x nbin x DSDmoment (7934,49,176,2)
freezing = DPR['NS']['VER']['heightZeroDeg'][:]
lwp = DPR['NS']['SLV']['precipWaterIntegrated'][:,:,0] #liquid water path
iwp = DPR['NS']['SLV']['precipWaterIntegrated'][:,:,1] #ice water path
precip_type = DPR['NS']['CSF']['typePrecip'][:]/10**7 #Convective vs stratiform 


#Add NaNs for zero values of iwp and lwp

lwp[lwp==0] = np.nan
iwp[iwp==0] = np.nan


#Compute the ratio of liquid water path to ice water path

ratio = lwp/iwp





###Define cross-section length

ind1 = np.where((lon[:,0]>=-94)) #Where are the DPR lons >= -100
ind2 = np.where((lon[:,0])<=-89) #Where are lons <= -85
ind3 = np.intersect1d(ind1,ind2) #Both conditions need to be satisfied here


#Change the ray to change cross-section location (i.e. the "27" in this case)
#%%

###Setup to 2D grid for plotting

x = 2.* 17 #48 degrees (from -17 to 17)
re = 6378. #radius of the earth
theta = -1*(x/2.) + (x/48.)*np.arange(0,49) #Split into equal degrees (from -17 to 17)
theta2  = np.ones(theta.shape[0]+1)*np.nan #Define an empty array (NaNs) with same shape as ray dimension
theta = theta - 0.70833333/2. #Shift degrees for plotting pruposes
theta2[:-1] = theta #remove last array dimension in theta so python won't get angry about shape
theta2[-1] = theta[-1] + 0.70833333
theta = theta2*(np.pi/180.) #Convert from degrees to radians

prh = np.ones((177,50))*np.nan #Define empty grid

for i in range(prh.shape[0]): #Loop over each range gate
    for j in range(prh.shape[1]): #Loop over each scan
            a = np.arcsin(((re+407)/re)*np.sin(theta[j]))-theta[j] #Orbit height of 407 km 
            
            prh[i,j] = (176-(i))*0.125*np.cos(theta[j]+a)
            
h2 = prh #where prh is the (range bins,ray)-space
h3 =np.ones((177,50))*np.nan

for i in range(h3.shape[1]):
    h3[:,i] = h2[::-1,i] #This reverses the vertical dimension so as indices increase, height increases

#%%
    
ku = DPR['NS']['SLV']['zFactorCorrected'][ind3,:,:] #Read in ku-band reflectivity; nscan x nray (554,49,176)
n0 = dsd[ind3,:,:,0]/10 #Read in the number concentration
d0 = dsd[ind3,:,:,1] #Read in the mean drop diameter  #Both have dimensions nscan x nray x nbin (554,49,176)
zeroDeg = freezing[ind3,:]

#Cut all parameters so they are at same ray as above

ray = 10

ku = ku[:,ray,:]
n0 = n0[:,ray,:]
d0 = d0[:,ray,:]
zero_deg_isotherm = zeroDeg[:,ray]/1000 #Convert from meters to km

#Take lats and lons along same ray
lons = DPR['NS']['Longitude'][ind3,ray]
lats = DPR['NS']['Latitude'][ind3,ray]




#Choose a starting point, then calculate distance
lat0 = lats[0]
lon0 = lons[0]


p = Proj(proj='laea', zone=10, ellps='WGS84',lat_0=lat0,lon_0=lon0) #Define a projection and set starting lat an lon to same point as above

#Define a 2D array for plotting purposes

lat_3d = np.ones(ku.shape)*np.nan
lon_3d = np.ones(ku.shape)*np.nan

for i in range(ku.shape[0]):
    lat_3d[i,:] = lats[i] 
    lon_3d[i,:] = lons[i]  
        

x,y = p(lon_3d,lat_3d) #Now convert degrees to distance (in km)
R_gpm = np.sqrt(x**2 + y**2)*np.sign(x) #Keeps sign of number; converts to radial distance 

#Reverse range gate order for all parameters

ku = ku[:,::-1]
n0 = n0[:,::-1]
d0 = d0[:,::-1]



ku = np.ma.masked_where(ku<=12, ku) #Mask all the bad points in ku data
y = np.ones([ku.shape[0], ku.shape[1]]) #Define an empty array

#Define the appropriate range bins
h4 = h3[:,ray] #This number MUST match the same ray being used
for i in range(y.shape[1]):
    y[:,i] = h4[i]
    
    
    
##Now compute 2D histogram for ku-band reflectivity moving along track
    
#%%  
   
def create_histogram(input_values, num_bins, min_value, max_value):
    """Creates a histogram with uniform bin-spacing.
    N = number of input values
    K = number of bins
    :param input_values: length-N numpy array of input values (to be binned).
    :param num_bins: Number of bins.
    :param min_value: Minimum value to include in histogram.  Any input value <
        `min_value` will be assigned to the first bin.
    :param max_value: Maximum value to include in histogram.  Any input value >
        `max_value` will be assigned to the last bin.
    :return: input_to_bin_indices: length-N numpy array of bin indices.  If
        input_values[i] = j, the [i]th input value belongs in the [j]th bin.
    :return: num_examples_by_bin: length-K numpy array, where the [j]th value is
        the number of inputs assigned to the [j]th bin.
    """

    bin_cutoffs = np.linspace(min_value, max_value, num=num_bins + 1)
    input_to_bin_indices = np.digitize(
        input_values, bin_cutoffs, right=False) - 1
    input_to_bin_indices[input_to_bin_indices < 0] = 0
    input_to_bin_indices[input_to_bin_indices > num_bins - 1] = num_bins - 1

    num_examples_by_bin = np.full(num_bins, -1, dtype=int)
    for j in range(num_bins):
        num_examples_by_bin[j] = np.sum(input_to_bin_indices == j)

    return input_to_bin_indices, num_examples_by_bin     


ku_hist = np.zeros((40,176))


for i in range(ku.shape[1]):
    
    ku_hist[:,i] = create_histogram(ku[:,i], 40,10,50)[1]   



dm_hist = np.zeros((25,176))


for i in range(d0.shape[1]):
    
    dm_hist[:,i] = create_histogram(d0[:,i],25,0.5,3)[1]
    

nw_hist = np.zeros((35,176))


for i in range(n0.shape[1]):
    
    nw_hist[:,i] = create_histogram(n0[:,i],35,1.5,5)[1]
    
    
    
ku_norm = (ku_hist - np.nanmin(ku_hist))/(np.nanmax(ku_hist) - np.nanmin(ku_hist))

    
ku_hist[ku_hist==0] = np.nan 
ku_hist[ku_hist>25] = np.nan    
ku_norm[ku_norm<0.001] = np.nan

x_bins = np.arange(10,50,1)
y_bins = np.arange(0,21.988032,0.124932)



fig,ax = plt.subplots(figsize=(14,14)) 
plt.pcolormesh(x_bins,y_bins,ku_hist.T, cmap = 'turbo')

xtick_label_size = 24
ytick_label_size = 26
tick_label_size = 30
title_size = 30



xlabels = np.arange(5,50,step = 5)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels)


plt.xlim(10,50)
plt.ylim(0,12)
#plt.clim(0,0.1)
plt.clim(0,15)

plt.ylabel('Altitude (km)', size = 26)
plt.xlabel('Bin range (dBZ)', size = 26)

ax.set_xticklabels(['10','15','20', '25', '30', '35', '40', '45', '50'], size = xtick_label_size)
plt.yticks(size =ytick_label_size )        
plt.title('8/27 1247 UTC Along-Track KuPR', size = title_size)


cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
#cbar.set_label(label = '[Frequency]',size = tick_label_size)
cbar.set_label(label = '[Frequency]',size = tick_label_size)


plt.show()   
 
#%%

#Now for Dm:


dm_hist[dm_hist==0] = np.nan
dm_hist[dm_hist>=50] = np.nan

x_bins = np.arange(0.5,3,0.1)
y_bins = np.arange(0,21.988032,0.124932)



fig,ax = plt.subplots(figsize=(14,14)) 
plt.pcolormesh(x_bins,y_bins,dm_hist.T, cmap = 'turbo')

xtick_label_size = 24
ytick_label_size = 26
tick_label_size = 30
title_size = 30



xlabels = np.arange(0,3.25,0.25)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels)


plt.xlim(0.25,3.25)
plt.ylim(0,12)
plt.clim(0,30)

plt.ylabel('Altitude (km)', size = 26)
plt.xlabel('Bin range (mm)', size = 26)

ax.set_xticklabels(['0.25','0.5','0.75','1.0','1.25','1.5','1.75','2.0','2.25','2.5','2.75','3.0','3.25'], size = xtick_label_size)
plt.yticks(size =ytick_label_size )        
plt.title(r'8/27 1247 UTC Along-Track $D_{M}$ ', size = title_size)


cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = '[Frequency]',size = tick_label_size)

plt.show()   


#%% 


#Now for log10Nw


nw_hist[nw_hist==0] = np.nan
nw_hist[nw_hist>=50] = np.nan

x_bins = np.arange(1.5,5.0,0.1)
y_bins = np.arange(0,21.988032,0.124932)



fig,ax = plt.subplots(figsize=(14,14)) 
plt.pcolormesh(x_bins,y_bins,nw_hist.T, cmap = 'turbo')

xtick_label_size = 24
ytick_label_size = 26
tick_label_size = 30
title_size = 30



xlabels = np.arange(0.5,5.75,0.5)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels)


plt.xlim(1,5.5)
plt.ylim(0,12)
plt.clim(0,40)

plt.ylabel('Altitude (km)', size = 26)
plt.xlabel(r'Bin range (mm $m^{-3}$)', size = 26)

ax.set_xticklabels(['1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0','5.5',''], size = xtick_label_size)
plt.yticks(size =ytick_label_size )        
plt.title(r'8/27 1247 UTC Along-Track $log_{10}(N_{W})$ ', size = title_size)


cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = tick_label_size)
cbar.set_label(label = '[Frequency]',size = tick_label_size)

plt.show() 


   
#%%    
#Take the transpose so it's row = bins, column = height
    
ku_hist_no_nan = ku_hist.copy() 
 

#NaN out zeros:

ku_hist[ku_hist==0] = np.nan

#%%

fig,ax = plt.subplots(figsize=(14,14))
ku[ku<0] = np.nan

ku_25 = np.nanpercentile(ku,25,axis = 0)
ku_50 = np.nanpercentile(ku,50,axis = 0)
ku_75 = np.nanpercentile(ku,75,axis = 0)
ku_95 = np.nanpercentile(ku,95,axis = 0)



x_bins = np.arange(10,60,1)
y_bins = np.arange(0,21.988032,0.124932)



fig,ax = plt.subplots(figsize=(14,14)) 

xtick_label_size = 20
ytick_label_size = 20
tick_label_size = 26
title_size = 28



xlabels = np.arange(10,65,5)
ylabels = np.arange(0,13.5,0.5)

plt.xticks(xlabels)
plt.yticks(ylabels)


plt.xlim(10,55)
plt.ylim(0,12)


plt.ylabel('Altitude (km)', size = 26)
plt.xlabel('dBZ', size = 26)

ax.set_xticklabels(['10','15','20', '25', '30', '35', '40', '45', '50','55','60'], size = xtick_label_size)
plt.yticks(size =ytick_label_size )        
plt.title('8/27 1247 UTC Along-Track KuPR Percentiles', size = title_size)


zero_deg_mean = np.nanmean(zero_deg_isotherm)

plt.plot(ku_25,y_bins, label = '25p', color = 'k')
plt.plot(ku_50,y_bins,label = '50p', color = 'r')
plt.plot(ku_75,y_bins,label = '75p', color = 'b')
plt.plot(ku_95,y_bins,label = '95p', color = 'g')
plt.axhline(y=zero_deg_mean,color = 'k', linestyle = '--')
plt.legend(fontsize = 'xx-large')

plt.show()

#%%


#Now index all points across all percentiles below the mean zero degree isotherm


warm_cloud_index = np.where((y_bins<4)&(y_bins>2))[0]


ku_25_warm = ku_25[warm_cloud_index]
ku_50_warm = ku_50[warm_cloud_index]
ku_75_warm = ku_75[warm_cloud_index]
ku_95_warm = ku_95[warm_cloud_index]

y_bins_warm = y_bins[warm_cloud_index]


#Now perform linear regression to compute slope across each percentile


ku25_slope = stats.linregress(ku_25_warm, y_bins_warm)
ku50_slope = stats.linregress(ku_50_warm, y_bins_warm)
ku75_slope = stats.linregress(ku_75_warm, y_bins_warm)
ku95_slope = stats.linregress(ku_95_warm, y_bins_warm)







#%%    


#Remove the values less than or equal to zero

n0[n0<=0] = np.nan
d0[d0<=0] = np.nan

    
    
#%% 
    
#Now we plot! #N-S (along track first)  


#Plot mean drop size first

dist = np.array([0,600])

x_loc = 500
y_loc = 14.5
label = '@NOAABrauer'

plt.figure(figsize=(10,10))

vmax = 3
vmin = 0

R_min = R_gpm.min()
R_max = R_gpm.max()

label_size = 20

cmin =0.; cmax = 3.; cint = 0.1; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, d0, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 8/27 1247 UTC $D_{M}$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(0,15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[mm]',size = label_size)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.clim(0,3)


#plt.text(x_loc, y_loc, label)


plt.show()    
    


####Number concentration (liquid water content)

plt.figure(figsize=(10,10))

vmax = 6
vmin = 1

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =1.; cmax = 6.; cint = 0.25; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, n0, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 8/27 1247 UTC $log_{10 }(N_{W})$ ', size = 20)
plt.xlim(dist[0],dist[1])
plt.ylim(0,15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = r'[$m^{-3}$ $mm^{-1}$]',size = label_size)
plt.clim(1,6)
#plt.text(x_loc, y_loc, label)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

plt.show()    
    

###And lastly Ku-band


plt.figure(figsize=(10,10))

vmax = 60
vmin =12

R_min = R_gpm.min()
R_max = R_gpm.max()



cmin =12.; cmax = 60.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='pyart_NWSRef',lut=nlevs)
pm = plt.pcolormesh(R_gpm/1000., y, ku, cmap='turbo',vmin=vmin,vmax=vmax)
pm2 = plt.plot(R_gpm/1000., zero_deg_isotherm, '--', color = 'k')
plt.xlabel('Along Track Distance (km)', size = 20)
plt.ylabel('Altitude (km)', size = 20)
plt.title(r'GPM Overpass 8/27 1247 UTC KuPR ', size = 20)
#plt.xlim(300,450)
plt.xlim(dist[0], dist[1])
plt.ylim(0,15)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = label_size)
cbar.set_label(label = '[dBZ]',size = label_size)
plt.clim(12,60)
plt.xticks(size = label_size)
plt.yticks(size = label_size)

#plt.text(x_loc, y_loc, label)
plt.show()   





#%%

#Determine median value of ind3 (intersection between longitudes)
cross_track_index = int(len(ind3)/2)+ind3[0]



#Let's plot a map with the center point location for our cross-sections

#Mask near surface reflectivity values

z = np.ma.masked_where(z<=12, z)



cmin = 12.; cmax = 70.; cint = 2.5; clevs = np.round(np.arange(cmin,cmax,cint),2)
nlevs = len(clevs) - 1; cmap = plt.get_cmap(name='turbo',lut=nlevs)
    
plt.figure(figsize=(10,10))
  
#xlim = np.array([-110,-75]); ylim = np.array([15,40])
xlim = np.array([-96,-89]); ylim = np.array([27,34])
 

  
m = Basemap(projection='cyl',lon_0=np.mean(xlim),lat_0=np.mean(ylim),llcrnrlat=ylim[0],urcrnrlat=ylim[1],llcrnrlon=xlim[0],urcrnrlon=xlim[1],resolution='i')
m.drawcoastlines(); m.drawstates(), m.drawcountries()
cs = m.contourf(lon,lat,z,clevs,cmap='turbo',extend='both')
cs2 = m.plot(lon[:,0]+0.03,lat[:,0]-0.03,'--k',zorder=4)
cs3 = m.plot(lon[:,-1]-0.03,lat[:,-1]+0.03,'--k',zorder=4)
m.plot(lon[ind3,ray],lat[ind3,ray],'-w',zorder=4,lw=0.25,alpha=0.75)


##Change this to modify star size (cross-section starting point)
m.plot(lon[ind3[0],ray],lat[ind3[0],ray],'*w',markersize = 30, zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
print(lon[ind3[0],ray])
print(lat[ind3[0],ray])


m.plot(lon[cross_track_index,:],lat[cross_track_index,:],'-w',zorder=4,lw=0.25,alpha=0.75)
m.plot(lon[cross_track_index,0],lat[cross_track_index,0],'*w', zorder=4,lw=0.25,alpha=0.75,markerfacecolor='k',markeredgewidth=0.25)
m.plot(lons[0:105],lats[0:105], 'k--', linewidth = 4)
m.drawcounties()
parallels = np.arange(15,25,step = 2)
m.drawparallels(parallels, labels = [True, False, False, False])

meridians = np.arange(-80,-70, step = 2)
m.drawmeridians(meridians, labels = [False, False, False, True])
cbar = m.colorbar(cs,size='2%')
cbar.ax.set_ylabel('[dBZ]',name='Calibri',size=18)
cticks = []
for i in clevs:
    cticks.append(int(i)) if i.is_integer() else cticks.append(i)
    cbar.set_ticks(clevs[::4])
    cbar.set_ticklabels(cticks[::4])
for i in cbar.ax.yaxis.get_ticklabels():
    i.set_family('Calibri')
    i.set_size(20)
   
plt.title('GPM Overpass 8/27 1247 UTC KuPR', size = 20)

markersize = 12


label = '1A'
x2star,y2star = m(-93.518746,30.224910)
m.plot(x2star,y2star,'or',markersize=markersize, color = 'r')
plt.text(x2star-0.35,y2star-0.15,label, weight = 'bold', size = 16)

label = '1B'
x2star,y2star = m(-93.359979,30.266642)
m.plot(x2star,y2star,'or',markersize=markersize, color = 'r')
plt.text(x2star+0.06,y2star+0.05,label, weight = 'bold', size = 16)

label = '2A'
x2star,y2star = m(-92.928744,30.245168)
m.plot(x2star,y2star,'or',markersize=markersize, color = 'r')
plt.text(x2star-0.15,y2star-0.25,label, weight = 'bold', size = 16)

label = '2B'
x2star,y2star = m(-92.739657,30.246693)
m.plot(x2star,y2star,'or',markersize=markersize, color = 'r')
plt.text(x2star+0.1,y2star-0.05,label, weight = 'bold', size = 16)

