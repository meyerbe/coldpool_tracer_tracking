import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import operator
from classes import CPlife,COL 
from six.moves import cPickle
from collections import OrderedDict
EXPID =  'test1plus4K'# 'lindp2K'#
odir = '/nbi/ac/conv1/henneb/results/coldpool/'
ngp = 320

f = open(odir+EXPID+'/tempdata/CPlife.save', 'rb')
a = cPickle.load(f)
f.close()

xxm = {}
xxn = {}
xxo = {}
xxp = {}
c = {}
for k in a.keys():
 if list(a[k].age[l] for l in sorted(a[k].age.keys())) != sorted(list(a[k].age[l] for l in sorted(a[k].age.keys())))  : 
    print a[k].age 
    print 'CP ', k
 else:
    for kk in a[k].age.keys():
      if a[k].age[kk] in xxm.keys():
        xxm[a[k].age[kk]] += a[k].noIT[kk]
        xxn[a[k].age[kk]] += a[k].noGP[kk]
        xxo[a[k].age[kk]] += a[k].noIGP[kk]
        xxp[a[k].age[kk]] += float(a[k].noIT[kk])/float(a[k].noT[kk])
        c[a[k].age[kk]] +=1
      else:
        xxm[a[k].age[kk]] = a[k].noIT[kk]
        xxn[a[k].age[kk]] = a[k].noGP[kk]
        xxo[a[k].age[kk]] = a[k].noIGP[kk]
        xxp[a[k].age[kk]] = a[k].noIT[kk]/a[k].noT[kk]
        c[a[k].age[kk]] = 1

for k in xxm.keys():
   print 'k', k
   print 'x', xxm[k]
   print 'c', c[k]
   xxm[k] = xxm[k]/c[k]
   xxn[k] = xxn[k]/c[k]
   xxo[k] = xxo[k]/c[k]
   xxp[k] = xxp[k]/c[k]

fig, ax = plt.subplots(4)
ax[0].plot(xxm.keys(),xxm.values())
ax[1].plot(xxn.keys(),xxn.values())
ax[2].plot(xxo.keys(),xxo.values())
ax[3].plot(xxo.keys(),xxp.values())

ax[3].set_xlabel('CP age/ 5 min')
ax[0].set_ylabel('# not collided tracer')
ax[1].set_ylabel('# gp occupied by CP edge')
ax[2].set_ylabel('# gp occupied by one CP edge')
ax[3].set_ylabel('% of not collided tracer')


#volume / '+'$\mathregular{10^3}$'+'$\mathregular{m^{-3}}$')
#for k in sorted(a.keys()):
#  xx = list(a[k].age.values())
#  yy = list(a[k].noIT.values())
#  ax[1].plot(xx,yy)
plt.show()

#######################################################################

#######################################################################

#f = open(odir+EXPID+'/tempdata/colli.save', 'rb')
#c = cPickle.load(f)
#f.close()
#
#for k in c.keys():
#  #print k 
# if len(k) > 1 :
#  print k, len(k)  #, ":",c[k].ntot
#  
#fig, ax = plt.subplots(1)
#for k in c.keys():
# if len(k) > 1 : # more than one CP  
#  yy = list(c[k].ntot.values())
#  xx = range(0,len(list(c[k].ntot.values())),1) #  list(a[k].noIT.values())
#  ax.plot(xx,yy)
#
#
##fig, ax = plt.subplots(4)
##ax[0].plot(range(1,len(c.ntot.keys),1),c.ntot.values())
#
##ax[0].set_xlabel('CP age/ 5 min')
##ax[0].set_ylabel('# not collided tracer')
#plt.show()






#if list(a[k] for k in a.keys().age[l] for l in sorted(a[k] for k in a.keys().age.keys())) != sorted(list(a[k] for k in k.keys().age[l] for l in sorted(a[k] for k in k.keys().age.keys())))  :
#  print 'CP ', k

