import numpy as np
from calcul import LatticePolymer
import visualisation
import matplotlib as mpl
from calcul import MonteCarlo
from copy import copy
from numpy import linalg

# Params
N = 1000            # Number of monomers
beta_eps = -0  # beta*eps
n = 100              # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n,N)
mcgroup.rosenbluth(perm=False)
groupPos = np.array(mcgroup.history[0].pos)
#print(mcgroup.compute_re())
for i in range(1,n):
    groupPos = np.vstack((groupPos,mcgroup.history[i].pos))
#print(groupPos)
groupPos = np.array(groupPos).T


# Generating polymer
# polymer = LatticePolymer(N, constraint = "force", beta_eps=beta_eps)
# polymer.gen_walk()
# length= polymer.length()
# print(length)
# pos = np.array(polymer.pos).T



#visualisation.singPolyVisu3D(groupPos[0], groupPos[1], groupPos[2])

# Visualisation of a group of polymers
xmax = max(groupPos[0])
ymax = max(groupPos[1]) 
zmax = max(groupPos[2])
xmin = min(groupPos[0])
ymin = min(groupPos[1]) 
zmin = min(groupPos[2])

map=np.zeros((xmax-xmin+1,ymax-ymin+1,zmax-zmin+1))
for i in range(len(groupPos[0])):
    x = groupPos[0][i]
    y = groupPos[1][i]
    z = groupPos[2][i]
    map[x-xmin-1,y-ymin-1,z-zmin-1]+=1

unique = np.unique(groupPos, axis=1)
points = np.zeros(len(unique[0]))

for i in range(len(unique[0])):
    points[i] = map[unique[0][i]-xmin-1,unique[1][i]-ymin-1,unique[2][i]-zmin-1]

visualisation.polyCloud3D(unique, points)
