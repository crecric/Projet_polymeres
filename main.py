import numpy as np
from calcul import LatticePolymer, MonteCarloFactory, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt

# Params
N = 2000         # Number of monomers
beta_eps = -0  # beta*eps
n = 10000              # Number of polymers
poly_per_run = 100
runs = 50
c_m = 0.3
c_p = 3

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarloFactory(n=n, N=N, beta_eps = beta_eps)
mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                      save='%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (runs, N, poly_per_run, c_m, c_p))
# mcgroup.rosenbluth(perm=True, c_m=c_m, c_p=c_p, relaxation=100000)
# print(mcgroup.history['origin'])
# print('Cloning in average every %f steps' % np.mean(mcgroup.c))
# print('Pruning in average every %f steps' % np.mean(mcgroup.k))
# print('%f clones were generated' % mcgroup.n_c)
# print('%f polymers were pruned' % mcgroup.n_p)

# positions = [pos[:N] for i, pos in enumerate(mcgroup.history['pos']) if pos.shape[0] >= N and mcgroup.history['origin'][i] < N]
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
ks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1]
ks = [int(N*k) for k in ks]
x = np.linspace(ks[0], ks[-1], 100)
ts = []
es = []
# mcgroup = MonteCarloFactory(load='%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (runs, N, poly_per_run, c_m, c_p))
for k in ks:
     t = mcgroup.compute_observable(LatticePolymer.length, k)
     e = mcgroup.error(LatticePolymer.length, k)
     ts.append(t)
     es.append(e)
     print("re(%d) = %f +/- %f" % (k, t, e))
     print('r_theo(%d)=' % k, r(k-1))

plt.errorbar(ks, ts, yerr=es, fmt='bo-')
plt.plot(x, r(x), 'r-')
plt.savefig('Re_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (runs, N, poly_per_run, c_m, c_p), dpi=300)
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
# visualisation.polyCloud3D(groupPos)
