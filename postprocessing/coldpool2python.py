###########################################################
# October 2018 OCH
# convert results from cold pool tracking into python data structure
#  * read ASCII output from cold pool tracking
#  * save no of CP, tracer etc in ghridpoint format CoPo 
#  * save collisions, CP information for combination of CPs which collide CPL 
#  * save CP information for single CP, including which are colliding with CP life
###########################################################
#
# TO DO
# * number of individual tracers (not colocated with tracer from other CP)
#   can increase when collided CPs terminate. How to handle this
import numpy as np
import os   # for environment 
import datetime
from classes import CP_map, COL, CPlife, COLDPOOL
import matplotlib.pyplot as plt
from six.moves import cPickle   # to save class files

EXPID = 'test1plus4K_circle'# 'lindp2K'#

atime = datetime.datetime.now()
odir = os.environ.get('results')+'/coldpool/'
f = open(odir+EXPID+'/output/cp/coldpool_tracer_out.txt', 'r')
lines = f.readlines()
nl = len(lines)
CoPo = {}
CPL = {}   #Cold Pool Collision
b={}
tist = 0
print "start first loop"
c = 0
for line in lines:
    c += 1
    print 'at', c, 'von', nl
    columns = line.split()
    tist = (int(columns[0]))
  #if tist < 166 and tist >160:
    age = (int(columns[1]))
    xpos = (int(columns[7]))
    ypos = (int(columns[6]))
    cCP  = (int(columns[3]))
    #print tist, xpos, ypos
    if (tist,xpos,ypos) not in CoPo.keys():
    #falsch: if (xpos,ypos,tist) != (k for k in CoPo.keys()) :
      #print 'make new at ',tist, xpos, ypos
      CoPo[tist,xpos,ypos]=CP_map() 
    CoPo[tist,xpos,ypos].add(cCP,age)  

btime = datetime.datetime.now()
print 'took ',(btime-atime), 'for loop '

################################################
# SAVE DATA
#####################

f = open(odir+EXPID+'/tempdata/CoPo.save', 'wb')
cPickle.dump(CoPo, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()

################################################
# Proceed
#####################

for (t,x,y) in CoPo.keys() :
    #print x,y,t     
    cCP = CoPo[t,x,y].CPs
    #print 'cps', list(cCP.values())
    # summerize collisions:
    if zip(*zip(cCP.values()))[0] not in CPL.keys():
        print list(cCP.values()), CPL.keys()
        #if this CP combination not yet exists, make it
        CPL[zip(*zip(list(cCP.values())))[0]]= COL()
    else:
        print x,y,t, 'list is added' 
    # add location of tracers 
    CPL[zip(*zip(list(cCP.values())))[0]].add(t,x,y,CoPo[t,x,y].nTrtot) 

    # summerize CP lifecycle
    for cp in list(cCP.values()):
      if cp not in b.keys():
        b[cp] = CPlife(cp,t)
      if len(list(cCP.values())) ==1 : 
        b[cp].add(CoPo[t,x,y].nTrCP[cp],CoPo[t,x,y].nTrCP[cp],1,t,CoPo[t,x,y].age[cp],x,y)  
      else:
        b[cp].add(CoPo[t,x,y].nTrCP[cp],0,0,t,CoPo[t,x,y].age[cp],x,y)
        for cps in list(cCP.values()):
         if not cps == cp:
          b[cp].add_others(t,cps,CoPo[t,x,y].nTrCP[cp],CoPo[t,x,y].nTrCP[cps] ) 


ctime = datetime.datetime.now()
print 'took ',(ctime-btime), 'for loop '

################################################
# SAVE DATA
#####################

f = open(odir+EXPID+'/tempdata/colli.save', 'wb')
cPickle.dump(CPL, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
f = open(odir+EXPID+'/tempdata/CPlife.save', 'wb')
cPickle.dump(b, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()


x = np.zeros(shape=(len(b) ,20  ))
y = np.zeros(shape=(len(b) ,20  )) #[]
c = 0
for k in b.keys():
 lala = b[k].age
 print lala 
 for kk in lala.keys():
  if lala[kk] < 20:
   print lala[kk]
   x[c,lala[kk]] = b[k].age[kk]  
   y[c,lala[kk]] = b[k].noT[kk]
 c += 1 
 

xx = list(b[54].age.values())
yy = list(b[54].noIT.values())
print xx
print yy

print c
fig, ax = plt.subplots()
for k in b.keys():
  xx = list(b[k].age.values())
  yy = list(b[k].noIT.values())
  ax.plot(xx,yy)

plt.show()











##    for c in cCP:
##      age = CoPo[t,x,y].age[c] 
##      if (c,age) not in b.keys():
##        b[c,age]=COLDPOOL(c,t,1)
##        
##      else:
##        b[c,age].add_marker(c,t,x,y)
##
##    ############
##    # find collision
##    ######
##    
##
##
##
##btime = datetime.datetime.now()
##print 'took ',(btime-atime), 'for loop '
##
##
#CoPo = {}
#CoPo[1,10,10] = 5
#CoPo[1,30,10] = 8
#CoPo[1,50,10] = 9
#CoPo[10,1,10] = 6
#a=1
#b=10
#c=10
#
#
#
#for k in CoPo.keys():
#  print 'k ' ,k
#if (a,b,c) != (k for k in CoPo.keys()) :
#    print 'yes' 
#if (a,b,c) in CoPo.keys(): 
#    print 'hier'
