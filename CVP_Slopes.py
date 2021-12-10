import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as dttime
import datetime
from scipy import stats
import matplotlib.dates as mdates







data = np.load('cvp_data.npy', allow_pickle = True).tolist()
ref = data['dbz'] #time x height shape: (80,240)
zdr = data['zdr']
height = data['height'] #shape 240


time = data['datetime']
datetimeformat = [dttime.num2date(date) for date in time]


#Remove first 11 elements of datetime objects in python


date_str = []


for i in range(len(datetimeformat)):
    
    date_time = datetimeformat[i]
    date = date_time.strftime('%H:%M:%S')
    
    
    date_str.append(date)
  


#Now let's calculate the slope
#Define warm cloud layer as 4.5 to 2 km

#Indices within the warm cloud layer we want to compute slopes for
height_index = np.where((height[0]<4500)&(height[0]>2000))[0]
height_warm = height[0][height_index]/1000

#Now set up an array of slopes for all times. Time has dimension of 80
#Loop through all times for the heights of interest
#First, will end up having an array of 80x49: time x depth
    
z_warm_layer = ref[:,height_index]  
zdr_warm_layer = zdr[:,height_index]  


slope_z = []
slope_zdr = []

for i in range(ref.shape[0]):
    
    ref_slope = stats.linregress(z_warm_layer[i,:],height_warm)[0]
    zdr_slope = stats.linregress(zdr_warm_layer[i,:],height_warm)[0]
    slope_z.append(ref_slope)
    slope_zdr.append(zdr_slope)


fontsize = 24

fig,ax = plt.subplots(figsize=(20,20))


fmt_hour = mdates.HourLocator(interval = 1)
ax.xaxis.set_major_locator(fmt_hour)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))




ax.plot(datetimeformat,slope_z, color = 'r', label = 'Slope of $Z_{H}$')
ax.plot(datetimeformat,slope_zdr, color = 'b', label = 'Slope of $Z_{DR}$')
#plt.hlines(0, date_str[0], date_str[79])
plt.xticks(size = fontsize)
plt.yticks(size = fontsize)
plt.xlabel('Time (UTC)', size = fontsize)
plt.ylabel('Slope', size = fontsize)
plt.legend(prop={'size': 16})
plt.title(r'4.5-2.0 km Slope of $Z_{H}$ and $Z_{DR}$', size = 30)
plt.show()




        



