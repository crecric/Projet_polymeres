import numpy as np
from random import choice

def random_step(N, coords):
    '''
    This function generates a random step (with periodic boundary conditions) at position pos in a 3-dimension n-length chain.
    N**3 is the number of accessible sites.
    '''
    x, y, z = coords
    neighbors = [[(x+1)%N, y, z], [(x-1)%N, y, z], [x, (y+1)%N, z], [x, (y-1)%N, z], [x, y, (z+1)%N], [x, y, (z-1)%N]]
    x, y, z = choice(neighbors)

    return x, y, z

def number_neighbors(N, coords, history):
    '''
    This function computes the number of occupied neighbors at a given site in a history of steps. 
    It is relevent to use in the case where a polymer chain can no longer be extended (the final step is surrounded by 6 neighbors)
    '''
    x, y, z = coords
    neighbors = [[(x+1)%N, y, z], [(x-1)%N, y, z], [x, (y+1)%N, z], [x, (y-1)%N, z], [x, y, (z+1)%N], [x, y, (z-1)%N]]
    c = 0
    for neighbor in neighbors:
        if neighbor in history:
            c += 1
    return c