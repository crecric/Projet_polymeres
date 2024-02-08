import numpy as np
import matplotlib.pyplot as plt

def visu3D(x,y,z):
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(x, y, z)
    plt.show()
