import numpy as np
from calcul import LatticePolymer, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt

# Params
N = 1000         # Number of monomers
beta_eps = -0  # beta*eps
n = 500              # Number of polymers
# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n, N, beta_eps = beta_eps)
mcgroup.rosenbluth(perm=True, c_m=0.2)
# groupPos = np.array(mcgroup.history['pos'][0])
# posi = [pos[:N] for pos in mcgroup.history['pos'] if pos.shape[0] >= N]
# print(len(posi), len(mcgroup.weights[N-1]))
# for i in range(1,n):
#      groupPos = np.vstack((groupPos,mcgroup.history['pos'][i]))
# groupPos = np.array(groupPos).T
# groupweight = np.array(mcgroup.history['weight'])
# print(groupweight)
def r(L):
     nu = 3/5
     return L**(2*nu)
ks = [N]
# x = np.linspace(80, 110, 100)
y = []
for k in ks:
     t = mcgroup.compute_re(k-1)
     y.append(t)
     print("re(%d):" % k, mcgroup.compute_re(k-1))
     print('r_theo(%d):' % k, r(k-1))

# plt.plot(ks, y, 'bo-')
# plt.plot(x, r(x), 'r-')
# plt.show()
# plt.savefig('r(L).jpg', dpi=300)
# print("Z :", mcgroup.Z)

# Generating polymer
# polymer = LatticePolymer(N, constraint = "force", beta_eps=beta_eps)
# polymer.gen_walk()
# length= polymer.length()
# print(length)
# pos = np.array(polymer.pos).T


# Visualisation of a single polymer
#visualisation.singPolyVisu3D(pos[0], pos[1], pos[2])

# Visualisation of a group of polymers
# visualisation.polyCloud3D(groupPos)
