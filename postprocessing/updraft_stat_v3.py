import numpy as np
import numpy.ma as ma
import os   # for environment 
import math
import matplotlib
import operator
import matplotlib.pyplot as plt
import datetime
from classes import RAIN, Marker
from operator import itemgetter 
from netCDF4 import Dataset   
from six.moves import cPickle   # to save class files
EXPID = 'test'

atime = datetime.datetime.now()
# SETTINGS
back = 18   # serach 'back' timesteps back in time for highest updraft

# data input
odir = '/nbi/ac/conv1/henneb/results/coldpool/'
print odir+'/'+EXPID+'/output/raincell/irt_tracks_mask.nc'
data_in_precip_mask   = Dataset(odir+'/'+EXPID+'/output/raincell/irt_tracks_mask.nc', 'r')
data_in_w = Dataset('/nbi/ac/conv1/henneb/Moseley2018/data/lind_p2K/lind_p2K.out.vol.w.nc', 'r')

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
#find level
level_index   = np.where( np.logical_and(data_in_level>=1000, data_in_level<=1200))[0] 
dim_z         = len(level_index)

# allocate again
data_w    = np.zeros(shape=(back+1,dim_x,dim_y,dim_z))

# get variables
starttime = data_in_timemask[0] - data_in_timew[0] #timestep of w at the time when tracking begins
data_w   = np.array(data_in_w.variables["w"][starttime-2-back:starttime-1,:,:,level_index[0]:level_index[dim_z-1]]) #level_index]) #[140,:,:,level_index])
print starttime

###############################################################
# find the initiation of precip 
###############################################################

# list of rain cell IDs
ID_list = {}
ID_list[0] = 0
ID_list[-1] = -1
# class object containing all information about rain cell (eg starttime and location)
a = []
btime = datetime.datetime.now()
print 'start main loop to find preicp initiation'
print 'took ',(btime-atime), 'for preparation'
for ti in range(0,min([dim_tm-1,dim_tw-1]),1): 
  print ti, ' von ' , dim_tm
  data_pmask = np.array(data_in_precip_mask.variables["var1"][ti,:,:])
  current_CP_IDs  = np.unique(data_pmask).astype(int)
  # update wind data to current timestep 
  data_w = np.roll(data_w,-1,axis=0)
  data_w[back,:,:,:] = np.array(data_in_w.variables["w"][starttime+ti,:,:,level_index[0]:level_index[dim_z-1]]) 
  dw = data_w[:-1,:,:,:]-data_w[1:,:,:,:]
  # go trough all precip Obj at starting current timestep
  for IDi in current_CP_IDs:
    if IDi not in ID_list: 
     indices = np.where(data_pmask == IDi)
     mask = data_pmask == IDi
     mask = mask[0,:,:]
     tmax   = np.max(np.max(np.max(data_w,axis=3)*mask,axis=1),axis=1)
     dtmax  = np.max(np.max(np.max(dw,axis=3)*mask,axis=1),axis=1)

     # the timestep when precip is initiated
     start = np.argmax(dtmax)  # relative to first precipitation time
     tstart = starttime+ti - start+data_in_timew[0]
     # the location where precip is initiated, where w is max at start time
     xstart  = np.where(np.max(data_w[start,:,:,:],axis=2)*mask == np.max(np.max(np.max(data_w[start,:,:,:],axis=2)*mask))) #[1]
     a.append(RAIN(int(start),int(tstart),int(xstart[0][0]),int(xstart[1][0]),int(IDi),indices[1],indices[2])) 
     # update the list of rain tracks already used
     ID_list[IDi] = IDi

ctime = datetime.datetime.now()
print 'took ',(ctime-btime), 'for loop to find updraft location'

###############################################################
# find the associated cold pools
############################################################
f = open(odir+EXPID+'/output/cp/coldpool_tracer_out.txt', 'r')
lines = f.readlines()
print ' run with value range'
#c = 0
b={}
for line in lines:
#    c += 1
    #line = f.next().strip()
    columns = line.split()
    tist = (int(columns[0]))
    xpos = (int(columns[7])) 
    ypos = (int(columns[6]))
    cCP  = (int(columns[3]))
    if cCP not in b.keys():
      b[cCP]=Marker(cCP)
    b[cCP].add_marker(tist,xpos,ypos)
#    print c, tist 
    for data in a:
      if tist == data.ts and xpos in data.xps and ypos in data.yps[np.where(data.xps == xpos)]:
         print 'halo', data.ts, data.x, data.y, data.ID, data.number   
         if cCP not in a[data.number].CPs:
           print 'yes'
           a[data.number].add_CP(cCP) 
           #a[data.number].noCPs = len(a[data.number].CPs)

dtime= datetime.datetime.now()
print 'took ',(dtime-ctime), 'for loop to find associated CPs'
f.close()

###############################################################
# get precip values for CELLS
#############################################################
f = open(odir+EXPID+'/output/raincell/irt_tracks_output_pure.txt', 'r')
lines = f.readlines()
rainID = 0
a.sort(key=operator.attrgetter('ID'))
allIDs = RAIN.IDs
print a[-1].ID
#for line in lines:
while rainID <  a[-1].ID-2: #  a..tolist().index(xpos)
  for line in lines:
    columns = line.split()
    rainID = (int(columns[0]))
    precip = (float(columns[5]))
    size = (float(columns[4]))
    if rainID in allIDs:
      a[allIDs.index(rainID)].add_PREC(precip,size)

#
#
f.close()
etime= datetime.datetime.now()

print 'took ',(etime-atime), 'for loop to find associated CPs'

######################################################
# save temorary data 
######################################################
os.mkdir('/nbi/ac/conv1/henneb/results/coldpool/'+EXPID+'/tempdata/')
for data in a:
 print data.CPs
 data.noCPs = len(data.CPs)
f = open('/data.save', 'wb')
cPickle.dump(a, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
f = open('/RAIN.save', 'wb')
cPickle.dump(RAIN, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
f = open('/Marker.save', 'wb')
cPickle.dump(b, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()



