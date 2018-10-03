import numpy as np
a = np.zeros(shape=(112,200))

f = open('../output/cp/coldpool_tracer_out.txt', 'r')
firstline = f.readline()
columns0 = firstline.split()
timestep0 = columns0[0]
age = columns0[1]
print(type(a[1,1]))
print(type(columns0[11]))
a[1,1] = columns0[11] 
k=1
print(type(k))
for line in f:
    line = line.strip()
    columns = line.split()
    timestep = columns[0]
    if timestep == timestep0:
#      print(type(a))
#      print(type(columns[11]))
      a[k,1] = a[k,1] +1 # columns[11]
    else:
      k = k+1
print(a[1,1])
f.close()
