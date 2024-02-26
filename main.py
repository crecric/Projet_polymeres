import numpy as np
from calcul import LatticePolymer, MonteCarlo
import visualisation
import matplotlib as mpl

# Params
N = 100         # Number of monomers
beta_eps = -0  # beta*eps
n = 10              # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n, N)
mcgroup.rosenbluth(perm=False, c_m=1)
groupPos = np.array(mcgroup.history['pos'][0])
for i in range(1,n):
     groupPos = np.vstack((groupPos,mcgroup.history['pos'][i]))
groupPos = np.array(groupPos).T
groupweight=np.array(mcgroup.history['weight'][0])
print(groupweight)
print(mcgroup.compute_re())

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
