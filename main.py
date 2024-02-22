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


# Visualisation of a single polymer
#visualisation.singPolyVisu3D(pos[0], pos[1], pos[2])

# Visualisation of a group of polymers

visualisation.polyCloud3D(groupPos)
