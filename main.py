import numpy as np
from calcul import LatticePolymer
from visualisation import visu3D
import matplotlib as mpl

# Params
n = 1000
beta_eps = -10e-23

# Generating polymer
polymer = LatticePolymer(n, constraint = "force", interacting = False,  beta_eps=beta_eps)
polymer.genwalk()
print(polymer.weight)
pos = np.array(polymer.pos).T

# Global plotting parameters
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (6,5)
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['axes.linewidth'] = 3

#visu3D(pos[0], pos[1], pos[2])