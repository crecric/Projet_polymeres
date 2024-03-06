import numpy as np
from calcul import LatticePolymer, MonteCarlo
import visualisation
import matplotlib as mpl

# Params
N = 1000         # Number of monomers
beta_eps = -0  # beta*eps
n = 1000              # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n, N, beta_eps = beta_eps)
mcgroup.rosenbluth(perm=True, c_m=1)
print("Nombre de cul-de-sacs :",mcgroup.failed)
groupPos = np.array(mcgroup.history['pos'][0])
for i in range(1,n):
     groupPos = np.vstack((groupPos,mcgroup.history['pos'][i]))
groupPos = np.array(groupPos).T
groupweight=np.array(mcgroup.history['weight'])
#print(groupweight)
print("re :", mcgroup.compute_re())

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
