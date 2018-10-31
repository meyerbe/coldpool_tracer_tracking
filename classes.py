import numpy as np

class CP_map():
    def __init__(self):
      # for one marker
      self.CPs = {}     # Colliding CPs
      self.age = {}
      self.nCPs = 0
      self.nTrtot = 0
      self.nTrCP = {} 

    def add(self,ID,age):
      if ID in self.CPs.keys():
        self.nTrCP[ID] += 1
      else:
        self.CPs[ID] = ID
        self.nTrCP[ID] = 1
        self.nCPs += 1
        self.age[ID] = age
      self.nTrtot += 1
########################################################
class COL():         # key are colliding CP combination
    def __init__(self):
      self.x = {}    # dict with timesteps as key giving a list 
      self.y = {}    # ...of all collision points for this CP combi
      self.ntot= {}  # ... of number of Tracers accumulating at collsion
    def add(self,t,x,y,ntot):
      if not t in self.x.keys(): 
        self.x[t]    = []
        self.y[t]    = []
        self.ntot[t] = 0
      self.x[t].append(x)
      self.y[t].append(y)
      self.ntot[t] += ntot

########################################################
class CPlife():
    def __init__(self,ID,t):
      self.ID      = ID     # ID
      self.start   = t
      self.age     = {}     # age of CP based on precip onset
                            # if CP result from merger tracer may have different ages
      self.age2nd  = {}     # CP age is than given by oldest precip event
      self.noT     = {}     # no of tracer
      self.noGP    = {}     # no of gp occupied from this CP
      self.noIT    = {}     # no of individual tracer (not colldided with other CP)
      self.noIGP   = {}     # no of individual grid point
      self.x       = {}
      self.y       = {}
      self.combi   = {} 
    def add(self,noT,noIT,noIGP,t,age,x,y):
      if t in self.noGP.keys():
        self.noGP[t]  += 1
        self.noT[t]   += noT
        self.noIT[t]  += noIT
        self.noIGP[t] += noIGP
        if age < self.age[t]:
          self.age2nd[t] = age
        else:
          self.age[t]= age        
        self.x[t].append(x) 
        self.y[t].append(y)
      else:
        self.noGP[t]  = 1
        self.noT[t]   = noT
        self.noIT[t]  = noIT
        self.noIGP[t] = noIGP
        self.age[t]   = age
        self.x[t] = []
        self.y[t] = []

    def add_others(self,t,cp,n,nother):
      if (t,cp) in self.combi.keys():
        self.combi[t,cp] += n 
      else: 
        self.combi[t,cp] = n 

      #self.age[t]   = age    # it can happen, that age of individual tracers of CP differe when precip events merge 
                              # and tracers for same Cp where set at different precip events
      #self.noIT    = noIT    # no of individual tracer (not colldided with other CP)
      #self.noIGP   = noIGP   # no of individual grid point
########################################################
class CP_terminate():
    def __init__(self,ID):
      self.ID = ID
      self.x = {}
      self.y = {}
    def add(self,t,x,y):
      if not t in self.x.keys():
        self.x[t] = []
        self.y[t] = []
      self.x[t].append(x)
      self.y[t].append(y)
      
    

class Pool():
    def __init__(self,ID):
      # for one marker
      self.ID = ID  # ID of CP marker belongs to 
      self.CPs = [] # Colliding CPs
      self.size = [] 
      self.ts = []
#      self.x = []   # 
#      self.y = []
    def add(self,t,r):
      self.ts.append(t)   #timestep 
      self.size.append(r) # size (average over all marker dists)
#      self.x.append(xp)
#      self.y.append(yp)

###########################################################################
class COLDPOOL():
    def __init__(self,ID,tstart,dur):
      self.ID      = ID      # ID 
      self.tstart  = tstart  # timestep when CP begins
      self.dur     = dur     # duration of CP

      self.ts = [] # timestep 
      self.x = []  # position of tracer at timestep
      self.y = []

      self.ColCPs  = {}      # colliding CPs
      self.locx    = {}      # location of collision
      self.locy    = {}      # location of collision
      self.Colt    = {}      # time of collision
      self.ColN    = {}      # number of tracer colliding
      self.ColSum  = {}      # summed numbr of tracer collided

    def add_marker(self,ID,t,xp,yp):
      self.ts.append(t)
      self.x.append(xp)
      self.y.append(yp)
 
    def add_coll(self,ID,tist,xp,yp):
      self.ColCPs[ID]=ID
#      self.ColCPs.append(ID)
      self.Colt[ID]   = tist #.append(tist)
      self.locx[ID]   = [xp]
      self.locy[ID]   = [yp]
      self.ColN[ID]   = 1 #.append(1)
      self.ColSum[ID] =1
    def add_loc(self,ID,xp,yp):
      self.locx[ID].append(xp)
    def add_t(self,ID,tist):
      self.ColN[ID] += 1

###########################################################################

class Marker():
    def __init__(self,ID):
      # for one marker
      self.ID = ID # ID of CP marker belongs to 
      self.ts = [] # timestep 
      self.x = []  # position of tracer at timestep
      self.y = []
    def add_marker(self,t,xp,yp):
      self.ts.append(t)
      self.x.append(xp)
      self.y.append(yp)
###########################################################################
class Marker2():
    def __init__(self,ID):
      # for one marker
      self.ID  = ID # ID of CP marker belongs to 
      self.ts  = [] # timestep 
      self.x   = [] # position of tracer at timestep
      self.y   = []
      self.age = []
      self.dist= []
      self.phi = []
    def add_marker(self,t,xp,yp,a,r,p):
      self.ts.append(t)
      self.x.append(xp)
      self.y.append(yp)
      self.age.append(a)
      self.dist.append(r)
      self.phi.append(p)
    def add_others(self,size):
      self.size=size

class RAIN():
   count = 0
   xs = []
   ys = []
   IDs = []
   tss = []
   xps = []
   yps = []

   def __init__(self,dt,ts,x,y,ID,xps,yps):
     self.number = RAIN.count
     self.ts = ts  #time when precip is initiated
     self.dt = dt  #time delay from initiation to onset
     self.x  = x   # location of strongest updraft
     self.y  = y
     self.ID = ID  # ID of precipitation object/track
     self.CPs = [] # create new emty list for involved CPs
     self.noTCP = [] # no of tracers from current CP found in area
     #self.CP1s = [] 
     self.CP2s = []  
     self.noTCP2 = []
     self.noCPs = len(self.CPs)
     self.noCP2s = len(self.CP2s)
     self.tCP2 = [] # when are CP-tracer pass area of precip initiation
     self.xps = xps # all gp of current precip
     self.yps = yps 
     self.prec = [] #precip intensity for all timesteps during the event
     self.size = [] 
     self.dur = len(self.prec)  #duration of precip event
     self.vol = np.sum(np.multiply(self.prec,self.size))
     self.meanP = 0
     self.firstP = 0

     RAIN.count += 1
     RAIN.IDs.append(ID)
     RAIN.xs.append(x)
     RAIN.ys.append(y)
     RAIN.tss.append(ts)
   def add_CP(self,CP):
     self.CPs.append(CP)
     self.noTCP.append(1) # counts number of tracer for this CP
   def add_TCP(self,CP):
     #self.noTCP[self.CPs==CP] += 1
     print self.CPs
     tin = self.CPs.index(CP)
     print self.noTCP[tin] 
     self.noTCP[tin] += 1
   def add_PREC(self,precip,sizes):
     self.prec.append(precip)
     self.size.append(sizes)
     self.firstP = self.prec[0]
     self.meanP = sum(self.prec) /max([self.dur,1])
   #def add_earlierCP(self,CP,t):
   #  self.CP1s.append(CP)
   def add_laterCP(self,CP,t):
     self.CP2s.append(CP)
     self.noCP2s = len(self.CP2s)
     self.tCP2.append(t)
     self.noTCP2.append(1)
   def add_TCP2(self,CP):
     self.noTCP2[self.CP2s==CP] += 1

######################################################
class RAIN2():
   count = 0
   xs = []
   ys = []
   IDs = []
   tss = []
   xps = []
   yps = []

   def __init__(self,dt,ts,x,y,ID,xps,yps):
     self.number = RAIN2.count
     self.ts = ts  #time when precip is initiated
     self.dt = dt  #time delay from initiation to onset
     self.x  = x   # location of strongest updraft
     self.y  = y
     self.ID = ID  # ID of precipitation object/track
     self.CPs = [] # create new emty list for involved CPs
     self.age = []
     self.noTCP = [] # no of tracers from current CP found in area
     self.noCPs = len(self.CPs)
     self.tCP = [] # when are CP-tracer pass area of precip initiation
     self.xps = xps # all gp of current precip
     self.yps = yps 
     self.prec = [] #precip intensity for all timesteps during the event
     self.size = [] 
     self.cogx = []
     self.cogy = []
     self.dur = len(self.prec)  #duration of precip event
     self.vol = np.sum(np.multiply(self.prec,self.size))
     self.meanP = 0
     self.firstP = 0
     RAIN2.count += 1
     RAIN2.IDs.append(ID)
     RAIN2.xs.append(x)
     RAIN2.ys.append(y)
     RAIN2.tss.append(ts)
   def add_PREC(self,precip,sizes,comx,comy):
     self.prec.append(precip)
     self.size.append(sizes)
     self.cogx.append(comx)
     self.cogy.append(comy)
     self.firstP = self.prec[0]
     self.meanP = sum(self.prec) /max([self.dur,1])
   #def add_earlierCP(self,CP,t):
   #  self.CP1s.append(CP)
   def add_CP(self,CP,t,a):
     self.CPs.append(CP)
     self.noCPs = len(self.CPs)
     self.tCP.append(t)
     self.noTCP.append(1)
     self.age.append(a)
   def add_TCP(self,CP):
     self.noTCP[self.CPs==CP] += 1

######################################################

