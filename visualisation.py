import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from calcul import LatticePolymer
from calcul import MonteCarlo

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

def polyCloud3D(groupPos):

    # Creating the matrix to map the frequency of monomers situated on a certain position in the grid
    xmax = max(groupPos[0])
    ymax = max(groupPos[1]) 
    zmax = max(groupPos[2])
    xmin = min(groupPos[0])
    ymin = min(groupPos[1]) 
    zmin = min(groupPos[2])
    map=np.zeros((xmax-xmin+1,ymax-ymin+1,zmax-zmin+1))
    for i in range(len(groupPos[0])):
        x = groupPos[0][i]
        y = groupPos[1][i]
        z = groupPos[2][i]
        map[x-xmin-1,y-ymin-1,z-zmin-1]+=1

    # Taking every point where at least a polymer passed
    unique = np.unique(groupPos, axis=1)
    points = np.zeros(len(unique[0]))

    # Associating to each point its frequency value
    for i in range(len(unique[0])):
        points[i] = map[unique[0][i]-xmin-1,unique[1][i]-ymin-1,unique[2][i]-zmin-1]

    # Plotting the discrete heatmap
    ax = plt.figure().add_subplot(projection='3d')
    color = points**(1/6)
    size = points**(1/2)
    ax.scatter(unique[0], unique[1], unique[2], c=-color, cmap='inferno', s=size, alpha=0.5)
    
    ax.grid(True)
    plt.show()