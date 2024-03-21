import numpy as np
from calcul import LatticePolymer, MonteCarloFactory, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt

# Params
N = 500         # Number of monomers
beta_eps = -0  # beta*eps
n = 1000              # Number of polymers
# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarloFactory(n=n, N=N, beta_eps = beta_eps)
mcgroup.rosenbluth(perm=True, c_m=0.25, relaxation=5000000)
# print(mcgroup.history['origin'])
# print('Cloning in average every %f steps' % np.mean(mcgroup.c))
# print('Pruning in average every %f steps' % np.mean(mcgroup.k))



#      groupPos = np.vstack((groupPos,mcgroup.history['pos'][i]))
#groupPos = np.array(groupPos).T
groupweight = np.array(mcgroup.history['weight'])
print(groupweight)

def r(L):
     nu = 3/5
     return L**(2*nu)
ks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1]
ks = [int(N*k) for k in ks]
x = np.linspace(ks[0], ks[-1], 100)
y = []
for k in ks:
     t = mcgroup.compute_observable(LatticePolymer.length, k)
     y.append(t)
     print("re(%d):" % k, t)
     print('r_theo(%d):' % k, r(k-1))

plt.plot(ks, y, 'bo-')
plt.plot(x, r(x), 'r-')
plt.savefig('r(L).jpg', dpi=300)
plt.xlabel('Number of monomers')
plt.ylabel('R_e^2')
plt.show()
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
visualisation.polyCloud3D(mcgroup,N-50,n)
visualisation.polyCloud3D(mcgroup,N,n)