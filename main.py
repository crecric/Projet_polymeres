import numpy as np
from calcul import LatticePolymer, MonteCarlo
from visualisation import visu3D
import matplotlib as mpl

# Params
N = 100               # Number of monomers
beta_eps = 0          # beta*eps
n = 10                # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n=n, N=N, beta_eps=beta_eps)
mcgroup.rosenbluth(perm=True, c_m=1)
print(mcgroup.compute_re())

# Generating polymer
# polymer = LatticePolymer(N, constraint = "force",  beta_eps=beta_eps)
# polymer.gen_walk()
# pos1 = polymer.pos
# polymer.gen_walk(start=51)
# pos2 = polymer.pos
# print(pos1)
# length = polymer.gyration()
# print(length)
# pos = np.array(polymer.pos).T

# Global plotting parameters
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (6,5)
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['axes.linewidth'] = 3

# visu3D(pos[0], pos[1], pos[2])