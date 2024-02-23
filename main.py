import numpy as np
from calcul import LatticePolymer, MonteCarlo
from visualisation import visu3D
import matplotlib as mpl

# Params
N = 100            # Number of monomers
beta_eps = -10e-23  # beta*eps
n = 100             # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n)
mcgroup.rosenbluth(perm=True, c_m=1)
# print(mcgroup.compute_re())

# Generating polymer
# polymer = LatticePolymer(N, constraint = "force",  beta_eps=beta_eps)
# polymer.gen_walk()
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