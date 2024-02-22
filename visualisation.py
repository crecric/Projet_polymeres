import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def singPolyVisu3D(x, y, z):

    # Global plotting parameters
    mpl.rcParams['font.size'] = 15
    mpl.rcParams['figure.figsize'] = (6,5)
    # mpl.rcParams['text.usetex'] = True
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['font.serif'] = ['Times']
    mpl.rcParams['axes.linewidth'] = 3

    # Plot
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(x, y, z, 'ro-', mfc='red', mec='k', ms=6.0, lw=1)

    # Marking the first and last point
    ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'bo', mfc='blue', mec='k', ms=6.0)

    ax.grid(True)
    plt.show()

def polyCloud3D(unique, points):
    ax = plt.figure().add_subplot(projection='3d')
    color = points**(1/6)
    size = points**(1)
    ax.scatter(unique[0], unique[1], unique[2], c=-color, cmap='inferno', s=size, alpha=0.5)
    
    ax.grid(True)
    plt.show()