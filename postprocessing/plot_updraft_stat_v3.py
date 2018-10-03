import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from fncs/classes import RAIN, Marker
from six.moves import cPickle

f = open('data.save', 'rb')
a = cPickle.load(f)
f.close()
f = open('RAIN.save', 'rb')
RAIN = cPickle.load(f)
print RAIN.tss
for data in a:
  print data.noCPs
f.close()

allTimes = range(a[0].ts,a[-1].ts,1) #RAIN.tss

trange = range(a[0].ts,a[-1].ts,3) #range(30,50,3)
allP =np.zeros(len(trange))
someP=np.zeros([len(trange),4])
counter=np.zeros( len(trange))
somecounter =np.zeros([len(trange),4])
#prepare data 
for data in a:
#  if data.noCPs > 0:
  print data.ID, data.noCPs, data.CPs
  
#  'is ther precipitation? ', data.meanP
  for i in range(0,len(trange)-2,1):
    if trange[i] <= data.ts and data.ts < trange[i+1]:
      allP[i]    += np.sum(np.multiply(data.prec,data.size)) #data.vol #np.sum(data.prec)
      #allP[i]    += np.sum(np.multiply(data.prec,data.size)) #meanP #+= data.meanP
      counter[i] += 1
      someP[i,min([data.noCPs,3])] += np.sum(np.multiply(data.prec,data.size))#data.vol #np.sum(data.prec) #data.vol
      somecounter[i,min([data.noCPs,3])] += 1
allP = np.multiply(np.divide(allP,counter),3e-4) #change unit: mm/h * (200m)^2 -> 
someP= np.multiply(np.divide(someP,somecounter),3e-4)
fig, ax = plt.subplots()
ax.plot(trange, allP)
ax.plot(trange,someP[:,0], 'o')
ax.plot(trange,someP[:,1], 'o')
ax.plot(trange,someP[:,2], 'o')
ax.plot(trange,someP[:,3], 'o')

ax.set_xlim(135, 200)
ax.set_ylim(0, 80)
ax.set_ylabel('volume / '+'$\mathregular{10^3}$'+'$\mathregular{m^{-3}}$')

plt.show()
