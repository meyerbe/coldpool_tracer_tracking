import numpy as np
import numpy.ma as ma
import math
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# SETTINGS
back = 18   # serach 'back' timesteps back in time for highest updraft

# data input
data_in_precip_mask   = Dataset("../output/raincell/irt_tracks_mask.nc", 'r')
data_in_w = Dataset("/nbi/ac/conv1/henneb/Moseley2018/data/lind_p2K/lind_p2K.out.vol.w.nc", 'r')

# extracting the dimensions of the datafile
dim_x         = data_in_precip_mask.dimensions['x'].size
dim_y         = data_in_precip_mask.dimensions['y'].size
dim_tm        = data_in_precip_mask.dimensions['time'].size
dim_z         = data_in_w.dimensions['zm'].size
dim_tw        = data_in_w.dimensions['time'].size 
print 'time steps of mask' , dim_tm
print 'time steps of w' , dim_tw

# allocate fields
data_pmask        = np.zeros(shape=(dim_x,dim_y))
data_in_timemask  = np.zeros(shape=(dim_tm))
data_in_timew     = np.zeros(shape=(dim_tw))
data_in_level     = np.zeros(shape=(dim_z))

# get variables
data_in_timemask  = np.array(data_in_precip_mask.variables['time'][:]) #% of day
data_in_timew     = np.array(data_in_w.variables['time'][:]) #seconds
data_in_level     = np.array(data_in_w.variables['zm'][:])

# change units
data_in_timemask = np.round(data_in_timemask *24.*60./5.).astype(int)  # in timesteps
data_in_timew    = np.round(data_in_timew /(60.*5.)).astype(int)  # in timesteps
delta_ts = data_in_timew[0]-data_in_timemask[0]  #w data start delta_ts later than tracking (if <0 earlier)
#find level
level_index   = np.where( np.logical_and(data_in_level>=1000, data_in_level<=1200))[0] 
dim_z         = len(level_index)

# allocate again
data_w    = np.zeros(shape=(back+1,dim_x,dim_y,dim_z))

# get variables
starttime = data_in_timew[0] +delta_ts
data_w   = np.array(data_in_w.variables["w"][starttime-2-back:starttime-1,:,:,level_index[0]:level_index[dim_z-1]]) #level_index]) #[140,:,:,level_index])
print 'initial shape ' , data_w.shape
# loop through time
tmean = {}
tmax  = {}
dtmax = {}
dtmean = {}
tstart = {}
xstart = {}
ystart = {}
zstart = {}

ID_list = {}
ID_list[0] = 0
ID_list[-1] = -1

time_list = range(0,10,1) #range(0,min([dim_tm-1,dim_tw-1]),1)
for ti in time_list: 
  print ti, ' von ' , dim_tm
  data_pmask = np.array(data_in_precip_mask.variables["var1"][ti,:,:])
  current_CP_IDs  = np.unique(data_pmask).astype(int)
  # update wind data to current timestep 
  data_w = np.roll(data_w,-1,axis=0)
  data_w[back,:,:,:] = np.array(data_in_w.variables["w"][starttime+ti,:,:,level_index[0]:level_index[dim_z-1]]) 
  dw = data_w[:-1,:,:,:]-data_w[1:,:,:,:]
  print dw[1:2,0,0,0] 
#  data_w = ma.masked_where(data_w < 0, data_w)
  # data_w =np.where(data_w > 0, data_w, 0)
  for IDi in current_CP_IDs:
    if IDi not in ID_list: 
     indices = np.where(data_pmask == IDi)
     tmean[IDi] = np.mean(np.mean(data_w[:,indices[1],indices[2],:],axis=1),axis=1)
     dtmean[IDi] = np.mean(np.mean((dw[:,indices[1],indices[2],:]),axis=1),axis=1)

     tmax[IDi]  = np.max(np.max(data_w[:,indices[1],indices[2],:],axis=1),axis=1)
     dtmax[IDi] = np.max(np.max((dw[:,indices[1],indices[2],:]),axis=1),axis=1)

     # the timestep when precip is initiated
     start = np.argmax(np.argmax((dw[:,indices[1],indices[2],:]),axis=1),axis=1) 
     #tstart[IDi] = starttime+ti - start
     # the location where precip is initiated, where w is max at start time
     #xstart[IDi]  = np.where(data_w[start,:,:,:])[0]
     #ystart[IDi]  = np.where(data_w[start,:,:,:])[1]
     #zstart[IDi]  = np.where(data_w[start,:,:,:])[2]

     # update the list of rain tracks already used
     ID_list[IDi] = IDi
    
print('finished loop')
a = np.zeros(shape=(len(tmean),back+1))#,dim_z-1))
b = np.zeros(shape=(len(tmean),back+1))
d = np.zeros(shape=(len(dtmax),back))
e = np.zeros(shape=(len(dtmax),back))
t = np.zeros(shape=(len(dtmax),back+1))


c =0 
for key in tmean:
  a[c,:] = tmean[key]
  b[c,:] = tmax[key]
  d[c,:] = dtmax[key]
  e[c,:] = dtmean[key]
#  t[c,:] = tstart[key]
  c = c+1
#  if c>10:
#    return
#a = a[0:10,:]
#b = b[0:10,:]
#d = d[0:10,:]
#e = e[0:10,:]
#colors = t 
fig, ax = plt.subplots(2,2)
# make a little extra space between the subplots
fig.subplots_adjust(hspace=0.5)
ax[0,0].plot(np.argmax(a[0:10,:],axis=1),np.max(a[0:10,:],axis=1), 'o')
ax[0,0].plot(range(0,back+1,1), np.swapaxes(a[0:10,:],1,0))
ax[0,0].set_xlim(0, 18)
ax[0,0].set_xlabel('time before event')
ax[0,0].set_ylabel('average updraft in initial precip area')
ax[0,0].grid(True)

ax[1,0].plot(np.argmax(b[0:10,:],axis=1),np.max(b[0:10,:],axis=1), 'o')
ax[1,0].plot(range(0,back+1,1), np.swapaxes(b[0:10,:],1,0))
ax[1,0].set_xlabel('time before event')
ax[1,0].set_ylabel('max updraft in initial precip area')


ax[0,1].plot(np.argmax(e[0:10,:],axis=1),np.max(e[0:10,:],axis=1), 'o')
ax[0,1].plot(range(0,back,1), np.swapaxes(e[0:10,:],1,0))
ax[0,1].set_xlim(0, 17)
ax[0,1].set_xlabel('time before event')
ax[0,1].set_ylabel('avg updraft change')


ax[1,1].plot(np.argmax(d[0:10,:],axis=1),np.max(d[0:10,:],axis=1), 'o')
ax[1,1].plot(range(0,back,1), np.swapaxes(d[0:10,:],1,0))
ax[1,1].set_xlim(0, 17)
ax[1,1].set_xlabel('time before event')
ax[1,1].set_ylabel('max updraft change')


fig.savefig("test.png")
#plt.show()

fig2,ax2 = plt.subplots()
ax2.plot(np.argmax(e,axis=1),np.max(e,axis=1), 'o')

ax2.plot(np.argmax(a,axis=1),np.max(a,axis=1), 'o')
fig2.savefig("test2.png")
plt.show()





