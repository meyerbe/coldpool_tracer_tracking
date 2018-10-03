import numpy as np

class Marker():
    def __init__(self,ID):
      self.ID = ID
      self.ts = []
      self.x = []
      self.y = []
    def add_marker(self,t,xp,yp):
      self.ts.append(t)
      self.x.append(xp)
      self.y.append(yp)


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
     self.ts = ts
     self.dt = dt 
     self.x  = x
     self.y  = y
     self.ID = ID
     self.CPs = [] # create new emty list for involved CPs
     self.noCPs = len(self.CPs)
     self.xps = xps
     self.yps = yps
     self.noCPs = len(self.CPs)
     self.prec = []
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
   def add_PREC(self,precip,sizes):
     self.prec.append(precip)
     self.size.append(sizes)
     self.firstP = self.prec[0]
     self.meanP = sum(self.prec) /max([self.dur,1])

