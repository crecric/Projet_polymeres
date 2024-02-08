import numpy as np
import math
import calcul
import visualisation

l=100
n=2*l
M=np.zeros((n+1,n+1,n+1))
M[math.floor(n/2),math.floor(n/2),math.floor(n/2)]=1
#print(M)
i=math.floor(n/2)
j=math.floor(n/2)
k=math.floor(n/2)
x=np.zeros(l)
y=np.zeros(l)
z=np.zeros(l)
for p in range(l):
   x[p]=i
   y[p]=j
   z[p]=k
   i,j,k=calcul.polymerize(M,i,j,k)
    
print(M)

x[-1]=i
y[-1]=j
z[-1]=k

visualisation.visu3D(x,y,z)