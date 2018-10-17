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

EXPID = 'test1plus4K'# 'lindp2K'#

atime = datetime.datetime.now()

odir = os.environ.get('results')+'/coldpool/'
f = open(odir+EXPID+'/output/cp/coldpool_tracer_out.txt', 'r')
lines = f.readlines()

CPL = {}   #Cold Pool Collision
b={}

f = open(odir+EXPID+'/tempdata/CoPo.save', 'rb')
CoPo = cPickle.load(f)
f.close()

for (t,x,y) in CoPo.keys() :
    #print x,y,t     
    cCP = CoPo[t,x,y].CPs
    #print 'cps', list(cCP.values())
    # summerize collisions:
    if zip(*zip(cCP.values()))[0] not in CPL.keys():
        print list(cCP.values()) #, CPL.keys()
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


btime = datetime.datetime.now()
print 'took ',(btime-atime), 'for loop '

################################################
# SAVE DATA
#####################

f = open(odir+EXPID+'/tempdata/colli.save', 'wb')
cPickle.dump(CPL, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()
f = open(odir+EXPID+'/tempdata/CPlife.save', 'wb')
cPickle.dump(b, f, protocol=cPickle.HIGHEST_PROTOCOL)
f.close()


