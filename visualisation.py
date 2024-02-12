import numpy as np
import matplotlib.pyplot as plt

def visu3D(x, y, z):
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(x, y, z, 'ro-', mfc='red', mec='k', ms=6.0, lw=1)

    # Marking the first and last point
    ax.plot(x[[0,-1]], y[[0,-1]], z[[0,-1]], 'bo', mfc='blue', mec='k', ms=6.0)

    ax.grid(True)
    plt.show()