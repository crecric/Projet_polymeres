import numpy as np
from random import choice, uniform
from copy import deepcopy as copy
from tqdm import tqdm
from colorama import Fore, Style
from bisect import insort

class LatticePolymer:
    def __init__(self, N, constraint='force', beta_eps=0):
        '''
        Initializes a polymer chain that is to be simulated by a self-avoiding 
        random walk. 
        Initializes the weight according to the choice between an interacting and not interacting random walk.
        Parameters
        ----------
        N : int
            Polymer length
        constraint : string
            Thermo-mecanical constraint applied to the chain. Only 'force' implemented for the moment.
        beta_eps : float
            strength of interacting energy compared to inverse temperature. 
            If beta_eps = 0 there is no closest-non-paired neighbor interaction.
        '''
        # Setting Interaction 
        self.beta_eps = beta_eps
        if self.beta_eps == 0:
            self.interacting = False
        else:
            self.interacting = True
        
        self.N = N
        self.constraint = constraint
        if self.constraint not in ['force', 'length']:
            raise NotImplementedError('Please select constraint in ["force", "length"].')

    def gen_walk(self, start=1, perm=False, c_m=1):
        '''
        Generates a chain of random steps to simulate the polymer. It starts at the center of the grid.
        Parameters
        ----------
        start : int
            Number of monomers already present in the polymer
        perm : bool
            Tells if pruning/enriching is applied during sampling
        c_m : float
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        '''
        # Positioning the initial monomer
        if start == 1:
            self.pos = [[0, 0, 0]]
            self.weight = 0
            if self.interacting:
                # The first iteration will wrongly count -1 occupied neighboring sites (cf number_pairs)
                # The weight has to be balanced in accordance.
                self.weight = np.log10(np.exp(-self.beta_eps))
        
        self.status = 'survive'
        self.heatup = 0
        # Looping on the walk
        try:
            for step in tqdm(range(start, self.N)):
                self.update_weight()
                self.heatup+=1
                try:
                    # Finding the max weight
                    self.weights[step].append(self.weight)
                    w_max = np.array(self.weights[step]).max()
                    y = [x-w_max for x in self.weights[step]]

                    zfactor = np.sum(np.power(10, y)) 
                    trials = len(self.weights[step])
                    self.Z[step] = np.log10(1/(trials)) + np.log10(zfactor) + w_max
        
                except OverflowError:
                    print('%sOVERFLOWERROR passed%s' % (Fore.RED, Style.RESET_ALL))
                    pass
    
                if perm and self.heatup >= self.N//50 and step <= int(0.9*self.N):
                    # Pruning/enriching
                    self.control_weight(step, c_m)

                # Stoping the walk when it reaches a closed-loop of neighbors
                if self.number_neighbors() == 0:
                    self.failed += 1
                    self.status = 'killed'
                    break

                # Generating a new direction
                x, y, z = self.random_step()
                while [x, y, z] in self.pos:
                    # Generating new step if the step is already present in the history of steps
                    x, y, z = self.random_step()

                self.pos.append([x,y,z])

        # If control_weight prematuraly kills a polymer
        except BreakException:
            pass
        
        self.pos = np.array(self.pos)

    def number_neighbors(self):
        '''
        This function computes the number of free neighboring sites of the last visited site in a history of steps. 
        It is relevant to use in the case where a polymer chain can no longer be extended 
        (the final step is surrounded by 6 neighbors) and to calculate the weight.
        '''
        neighbors = self.neighborhood(self.pos[-1])
        c = 0
        for neighbor in neighbors:
            if neighbor not in self.pos:
                c += 1
        return c
    
    def random_step(self):
        '''
        This function generates a random step starting from last visited site.
        '''
        x, y, z = choice(self.neighborhood(self.pos[-1]))
        return x, y, z
    
    def control_weight(self, step, c_m):
        '''
        This function applies the pruning/enriching algorithm to the Rosenbluth sampling.
        '''
        # PERM parameters
        c_p = 10*c_m
        # Current estimator of partition function
        W_m = np.log10(c_m)+self.Z[step]
        W_p = np.log10(c_p)+self.Z[step]

        # Pruning
        if self.weight < W_m:
            if uniform(0, 1) < 0.5:
                print('%sPolymer killed!%s' % (Fore.RED, Style.RESET_ALL))
                self.status = 'killed'
                self.failed += 1
                raise BreakException()
            else:
                print('%sPolymer survived!%s' % (Fore.GREEN, Style.RESET_ALL))
                self.weight += np.log10(2)
        elif self.weight > W_p:
            self.weight -= np.log10(2)
            self.clones.append(self.checkpoint())
            print('%sPolymer has been cloned!%s' % (Fore.CYAN, Style.RESET_ALL))
            self.heatup = 0
    
    def checkpoint(self):
        '''
        This function saves the key properties of a polymer at any step of its growth.
        '''
        return {'weight': copy(self.weight), 'pos': copy(self.pos)}
    
    def reset(self, weight, pos):
        self.weight = weight
        self.pos = pos

    def update_weight(self, step):
        '''
        Updates weight according to the chosen random walk pattern.
        '''
        numb_neigh = self.number_neighbors()
        if numb_neigh !=0:
            if not self.interacting:
                self.weight += np.log10(numb_neigh)
            else: 
                numb_pairs = 5 - numb_neigh
                self.weight += np.log10(numb_neigh*np.exp(-self.beta_eps*numb_pairs))
        else:
            self.weight = 0

        # Calculation of Z (for Monte Carlo purposes)
        try:
            # Finding the max weight
            self.weights[step].append(self.weight)
            w_max = np.array(self.weights[step]).max()
            y = [x-w_max for x in self.weights[step]]

            zfactor = np.sum(np.power(10, y)) 
            trials = len(self.weights[step])
            self.Z[step] = np.log10(1/(trials)) + np.log10(zfactor) + w_max

        except OverflowError:
            print('%sOVERFLOWERROR passed%s' % (Fore.RED, Style.RESET_ALL))
            pass
        
    @staticmethod
    def length(pos):
        '''
        Computes the squared start-to-end length of a set of vectors.
        '''
        return np.linalg.norm(pos[-1]-pos[0], 2)**2
    
    @staticmethod
    def gyration(pos):
        '''
        Computes the radius of gyration (how "curled up" the polymer is) of a set of vectors.
        '''
        rCM = np.average(pos, axis=0)
        return np.average(np.linalg.norm(pos - rCM, ord=2, axis=1)**2)
    
    @staticmethod
    def neighborhood(r):
        '''
        Checks the neighboring sites of a given vector r.
        '''
        x, y, z = r
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        return neighbors


class MonteCarlo(LatticePolymer):
    '''
    Generates collection of polymers.
    Returns thermodynamic observables.
    '''
    def __init__(self, n=1, N=100, constraint='force', beta_eps=0):
        '''
        Parameters
        ----------
        n : int
          Number of monte carlo steps (number of generated polymers)
        '''
        self.failed = 0
        self.n = n
        LatticePolymer.__init__(self, N, constraint, beta_eps)
        self.weights = [[] for _ in range(self.N)]
        self.Z = np.zeros(shape=self.N)
        self.history = {'weight': [], 'pos': []}
        self.clones = []
        # self.clones will eventually take the form [clone0_properties, clone1_properties, ...]

    def rosenbluth(self, perm=False, **kwargs):
        '''
        Fills the history with the polymers simulated by a random walk with the normal Rosenbluth method.
        Parameters
        ----------
        perm : bool
            Tells if pruning/enriching is applied during sampling
        c_m : float
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        '''

        self.perm = perm
        c_m = kwargs.get('c_m', 1)   # lower threshold
        start = 3                    # pruning/enriching is only applied after some trials
        self.trial = 0

        while self.trial < self.n:
            print('Simulating polymer %d:' % self.trial)
            if self.trial < start or not self.perm:
                self.gen_walk(perm = False)
                
            else:
                # Cheking if a clone has been generated for this trial
                if self.clones: # and self.history[trial-1].status == 'killed':
                    clone = self.clones[-1]
                    m = len(clone['pos'])                           # Number of monomers already present in present polymer
                    self.reset(clone['weight'], clone['pos'])
                    self.gen_walk(m, perm = True, c_m = c_m)        # Processing polymer growth on top of the clone

                    self.clones.remove(clone)

                # Else generating polymer from scratch
                else:    
                    self.gen_walk(perm=True, c_m=c_m)

            # if self.status == 'survive':
            self.history['weight'].append(self.weight)
            self.history['pos'].append(self.pos) 
            self.trial += 1

    def compute_re(self, N):
        '''
        Computes the observable squared norm of a polymer at a given number of monomers.
        '''
        logweights = self.weights[N]
        trials = len(self.weights[N]) 
        positions = [pos[:N] for pos in self.history['pos'] if pos.shape[0] >= N]
        lengths = [self.length(pos) for pos in positions]
        weights = np.power(10, [w-np.log10(trials)-self.Z[N] for w in logweights])

        return np.average(lengths, weights=weights)

class BreakException(Exception):
    pass