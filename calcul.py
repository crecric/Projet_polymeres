import numpy as np
from random import choice
from copy import copy
from numpy import linalg

class LatticePolymer:
    '''
    TODO: excluse stuck polymers 

    ### This should already be the case since we give it a weight of 0. ###
    '''
    def __init__(self, N=100, constraint='force', interacting=False,  **kwargs):
        '''
        Initializes a polymer chain that is to be simulated by a self-avoiding random walk. 
        Creates a grid and initializes the weight according to the choice between interacting and not interacting random walk.
        '''
        self.interacting = interacting
        #self.set_weight()
        self.N = N
        self.constraint = constraint
        if self.constraint not in ['force', 'length']:
            raise NotImplementedError('Please select constraint in ["force", "length"].')
        self.beta_eps = kwargs.get("beta_eps",None)
        # Creating a sufficiently big grid
        grid = 2*self.N + 1
        self.weight = 1
        if self.interacting:
            # This initialization is necessary because the first iteration will wrongly count -1 occupied neighboring sites,
            # so the effect has to be inversed.
            self.weight = np.exp(-self.beta_eps)

    def genwalk(self, compute_weight=True):
        '''
        Generates a chain of random steps to simulate the polymer. It starts at the center of the grid. It creates
        '''

        # Positioning the initial monomer in the center of the grid
        self.pos = [[self.N, self.N, self.N]]
        
        # Looping on the walk
        for step in range(self.N):
            if compute_weight:
                
                
                self.update_weight()
                #print(self.weight)
            if self.number_neighbors() == 0:
                # Stoping the walk when it reaches a closed-loop of neighbors
                break
            x, y, z = self.random_step()
            
            while [x, y, z] in self.pos:
                # Generating new step if the step is already present in the history of steps
                x, y, z = self.random_step()

            self.pos.append([x,y,z])

        self.pos=np.array(self.pos)

    def number_neighbors(self):
        '''
        This function computes the number of free neighboring sites at a given site in a history of steps. 
        It is relevant to use in the case where a polymer chain can no longer be extended (the final step is surrounded by 6 neighbors), 
        or to calculate the weight.
        '''
        x, y, z = self.pos[-1]
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        c = 0
        for neighbor in neighbors:
            if neighbor not in self.pos:
                c += 1
        return c
    
    def number_pairs(self):
        '''
        This function computes the number of occupied neighbors at a given site in a history of steps. 
        It is relevent to use in the case where a polymer chain can no longer be extended (the final step is surrounded by 6 neighbors),
        or to calculate the weight.
        '''
        x, y, z = self.pos[-1]
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        # Counting the number of occupied neighbors
        c = -1  
        # Since we will certainly count the occupied neighbor from which we just moved, we start the count at -1. 
        # An adjustment was added in the __init__ function to account for the difference in treatment for the very first iteration,
        # where there is no occupied neighbor that we shouldn't count.
        for neighbor in neighbors:
            if neighbor in self.pos:
                c += 1

        return c
    
    def random_step(self):
        '''
        This function generates a random step at position pos in a 3-dimension n-length chain.
        N**3 is the number of accessible sites.
        '''
        x, y, z = self.pos[-1]
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        x, y, z = choice(neighbors)

        return x, y, z
    
    def update_weight(self):
        '''
        Updates weight according to the chosen random walk pattern
        '''
        if not self.interacting:
            self.weight *= self.number_neighbors()

        else: 
            self.weight *= self.number_neighbors()*np.exp(-self.beta_eps*self.number_pairs())


    def length(self):
        '''
        Computes the square length of a polymer (squared norm between beginning and end of said polymer)
        '''
        return np.linalg.norm(self.pos[-1]-self.pos[0],2)
    
    
    # def perm():
        
class MonteCarlo(LatticePolymer):
    '''
    Generates collection of polymers
    Returns thermodynamic observables
    '''
    def __init__(self, n=10, N=100, constraint='force', interacting=False,  **kwargs):
        self.n = n
        LatticePolymer.__init__(self)
        self.history = np.empty(shape=self.n, dtype=MonteCarlo)

    def rosenbluth(self, perm=False):
        '''
        Fills the history with multiple polymers simulated by a random walk with the normal Rosenbluth method.
        '''
        for trial in range(self.n):
            poly = copy(self)

            poly.genwalk()
            self.history[trial]=poly
    
    def compute_re(self):
        '''
        Computes the weighted average squared norm of a group of polymers
        '''
        
        return (np.sum([self.history[trial].weight*self.history[trial].length() for trial in range(self.n)]))/ \
                np.sum([self.history[trial].weight for trial in range(self.n)])
    
# polymer.sample_re(rosenbluth='perm')
# class Observables(LatticePolymer):