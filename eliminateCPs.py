# October 2018
# open cold pool life: 
#   get for every cp colliding cp and numer of common tracer gp, 
#   if 1/3 of CP tracer is sread over same gps as another tracer, terminate the tracers
#   therefore get tracer information by loading CP collisons
##########################################################
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import operator
from classes import CPlife,COL,CP_terminate 
from six.moves import cPickle
from collections import OrderedDict
EXPID =  'test1plus4K_circle'# 'lindp2K'#
odir =os.environ.get('results') + 'coldpool/'
ngp = 320

f = open(odir+EXPID+'/tempdata/CPlife.save', 'rb')
a = cPickle.load(f)
f.close()

f = open(odir+EXPID+'/tempdata/colli.save', 'rb')
CPL = cPickle.load(f)
f.close()

ne = 0
n = 0
cptermin = {}
perold = -1
for k in a.keys():
 for t,cp in a[k].combi.keys():
   print t,cp
   per = float(a[k].combi[t,cp])/float(a[k].noT[t]) 
   if per > 0.3:
     print 'greater 0.3'
     for kk in CPL.keys():
      if k in kk and cp in kk and t in CPL[kk].x.keys():
        print 'add', t,CPL[kk].x[t],CPL[kk].y[t]
        if not k in cptermin.keys():
         cptermin[k] = CP_terminate(k)
        cptermin[k].add(t,CPL[kk].x[t],CPL[kk].y[t])
     
f = open(odir+EXPID+'/output/cp/coldpool_tracer_out.txt', 'r')
lines = f.readlines()
fo = open(odir+EXPID+'/output/cp/coldpool_tracer_out_reduce.txt','wb')

for line in lines:
  columns = line.split()
  xpos = (int(columns[7]))
  ypos = (int(columns[6]))
  cp   = (int(columns[3]))
  t    = (int(columns[0]))
  if not cp in cptermin.keys():
    fo.write(line)
  else:  
#    print 'cp is in list', cp
    if not t in cptermin[k].x.keys():
      fo.write(line)
    else:
#     print 'at time', t
#     #if (xpos,ypos) in zip(cptermin[cp].x[t],cptermin[cp].y[t]): 
#     print xpos, cptermin[cp].x[t]
#     if xpos in cptermin[cp].x[t][0]:
#       print cptermin[cp].x[t][0]
#       print 'x found'
#     if ypos in cptermin[cp].y[t][0]:
#       print 'y found'
#      
     if xpos in cptermin[cp].x[t][0]: # and ypos in cptermin[cp].y[t][0]: # add where
       print cp, 'nearly stoped at' ,t, xpos, ypos
     else:
       if ypos in cptermin[cp].y[t][0]:
         print cp, ' stoped at' ,t, xpos, ypos
       else:
         fo.write(line)
 

   #if per <= perold: 
##   print 'combi', per, 'at ',t,a[k].age[t],':',k,'with',a[k].noT[t],'and',cp, 'with',a[k].combi[t,cp], 'combi', per 
#   print k, a[k].age[t]
#   n += 1
#   if  per > 0.3: ne+=1
  


print n, ne






#xxm = {}
#xxn = {}
#xxo = {}
#xxp = {}
#c = {}
#for k in a.keys():
# if list(a[k].age[l] for l in sorted(a[k].age.keys())) != sorted(list(a[k].age[l] for l in sorted(a[k].age.keys())))  : 
#    print a[k].age 
#    print 'CP ', k
# else:
#    for kk in a[k].age.keys():
#      if a[k].age[kk] in xxm.keys():
#        xxm[a[k].age[kk]] += a[k].noIT[kk]
#        xxn[a[k].age[kk]] += a[k].noGP[kk]
#        xxo[a[k].age[kk]] += a[k].noIGP[kk]
#        xxp[a[k].age[kk]] += float(a[k].noIT[kk])/float(a[k].noT[kk])
#        c[a[k].age[kk]] +=1
#      else:
#        xxm[a[k].age[kk]] = a[k].noIT[kk]
#        xxn[a[k].age[kk]] = a[k].noGP[kk]
#        xxo[a[k].age[kk]] = a[k].noIGP[kk]
#        xxp[a[k].age[kk]] = a[k].noIT[kk]/a[k].noT[kk]
#        c[a[k].age[kk]] = 1
#
#for k in xxm.keys():
#   print 'k', k
#   print 'x', xxm[k]
#   print 'c', c[k]
#   xxm[k] = xxm[k]/c[k]
#   xxn[k] = xxn[k]/c[k]
#   xxo[k] = xxo[k]/c[k]
#   xxp[k] = xxp[k]/c[k]
#
#fig, ax = plt.subplots(4)
#ax[0].plot(xxm.keys(),xxm.values())
#ax[1].plot(xxn.keys(),xxn.values())
#ax[2].plot(xxo.keys(),xxo.values())
#ax[3].plot(xxo.keys(),xxp.values())
#
#ax[3].set_xlabel('CP age/ 5 min')
#ax[0].set_ylabel('# not collided tracer')
#ax[1].set_ylabel('# gp occupied by CP edge')
#ax[2].set_ylabel('# gp occupied by one CP edge')
#ax[3].set_ylabel('% of not collided tracer')
#
#
##volume / '+'$\mathregular{10^3}$'+'$\mathregular{m^{-3}}$')
##for k in sorted(a.keys()):
##  xx = list(a[k].age.values())
##  yy = list(a[k].noIT.values())
##  ax[1].plot(xx,yy)
#plt.show()
#
########################################################################
#
########################################################################
#
##f = open(odir+EXPID+'/tempdata/colli.save', 'rb')
##c = cPickle.load(f)
##f.close()
##
##for k in c.keys():
##  #print k 
## if len(k) > 1 :
##  print k, len(k)  #, ":",c[k].ntot
##  
##fig, ax = plt.subplots(1)
##for k in c.keys():
## if len(k) > 1 : # more than one CP  
##  yy = list(c[k].ntot.values())
##  xx = range(0,len(list(c[k].ntot.values())),1) #  list(a[k].noIT.values())
##  ax.plot(xx,yy)
##
##
###fig, ax = plt.subplots(4)
###ax[0].plot(range(1,len(c.ntot.keys),1),c.ntot.values())
##
###ax[0].set_xlabel('CP age/ 5 min')
###ax[0].set_ylabel('# not collided tracer')
##plt.show()
#
#
#
#


#if list(a[k] for k in a.keys().age[l] for l in sorted(a[k] for k in a.keys().age.keys())) != sorted(list(a[k] for k in k.keys().age[l] for l in sorted(a[k] for k in k.keys().age.keys())))  :
#  print 'CP ', k

