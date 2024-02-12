import numpy as np
from random import choice

def random_step(N, pos):
    '''
    This function generates a random step (with periodic boundary conditions) at position pos in a 3-dimension n-length chain.
    N**3 is the number of accessible sites.
    '''
    x, y, z = pos
    x, y, z = choice([((x+1)%N, y, z), ((x-1)%N, y, z), (x, (y+1)%N, z), (x, (y-1)%N, z), (x, y, (z+1)%N), (x, y, (z-1)%N)])

    return x, y, z

