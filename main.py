import numpy as np
from calcul import LatticePolymer
from visualisation import visu3D
import matplotlib as mpl

# Grid and length params
n = 10000

# Generating polymer
polymer = LatticePolymer(n)
polymer.genwalk()
polymer.weight
pos = np.array(polymer.pos).T

# Global plotting parameters
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (6,5)
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['axes.linewidth'] = 3

visu3D(pos[0], pos[1], pos[2])