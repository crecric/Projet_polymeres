import numpy as np
from calcul import LatticePolymer
from visualisation import visu3D
import matplotlib as mpl
from calcul import MonteCarlo
from copy import copy
from numpy import linalg

# Params
N = 1000            # Number of monomers
beta_eps = -10e-23  # beta*eps
n = 10              # Number of polymers

# Generating a group of polymers with Rosenbluth method
mcgroup = MonteCarlo(n)
mcgroup.rosenbluth(perm=False)
print(mcgroup.compute_re())

# Generating polymer
# polymer = LatticePolymer(N, constraint = "force", interacting = False,  beta_eps=beta_eps)
# polymer.genwalk()
# length= polymer.length()
# print(length)
# pos = np.array(polymer.pos).T

# Global plotting parameters
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (6,5)
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['axes.linewidth'] = 3

#visu3D(pos[0], pos[1], pos[2])