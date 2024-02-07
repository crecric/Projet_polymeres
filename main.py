import numpy as np
import math
import calcul
import visualisation

l=4
n=2*l
M=np.zeros((n+1,n+1))
M[math.floor(n/2),math.floor(n/2)]=1
print(M)
k=math.floor(n/2)
for i in range(l):
    k=calcul.polymerize(M(k))