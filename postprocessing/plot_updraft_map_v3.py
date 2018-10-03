from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import operator
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from classes import Marker, RAIN
from six.moves import cPickle
import os,sys
from matplotlib import cm
sys.path.append(os.path.abspath('.'))
#import colormaps as cmaps  
odir = '/nbi/ac/conv1/henneb/results/coldpool/test/'  
f = open(odir+'/tempdata/data.save', 'rb')
a = cPickle.load(f)
f.close()
f = open(odir+'/tempdata/Marker.save', 'rb')
b = cPickle.load(f)
f.close()
cmap = plt.get_cmap('hsv') #viridis')
fig, ax = plt.subplots(1)

a.sort(key=operator.attrgetter('ts'))
current = a[0].ts-1
xy_mask = np.zeros(shape=(1024,1024))
#xy_mask = np.nan
#xy_mask[xy_mask==0] = np.ma.masked
for data in a:
  if data.ts  > current:
    fig.savefig(odir+'/plots/test'+str(data.ts)+'.png')
    fig, ax = plt.subplots(1)
    ax.set_xlim(0, 320)
    ax.set_ylim(0, 320)
    ax.set_aspect(aspect='equal')

  for xi,yi in zip(data.xps,data.yps):
     xy_mask[yi,xi] = data.noCPs+1   
  ax.contourf(np.arange(0,1024 , 1),np.arange(0, 1024, 1),xy_mask,levels=[0.5,1,2,3,4,5])
  xy_mask[:,:] = 0
  current = data.ts

  for cpno in data.CPs:
    x = np.array(b[cpno].x)[np.where(np.array(b[cpno].ts) == data.ts)[0]] 
    y = np.array(b[cpno].y)[np.where(np.array(b[cpno].ts) == data.ts)[0]]
    ax.plot(x,y,'o', markersize=3,color=cmap(b[cpno].ID))        


