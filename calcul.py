import numpy as np
from random import choice, uniform
from copy import deepcopy as copy
from tqdm import tqdm
from sklearn.preprocessing import normalize

class LatticePolymer:
    def __init__(self, N=100, constraint='force', beta_eps=0):
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
        
        self.weight = 1
        if self.interacting:
            # The first iteration will wrongly count -1 occupied neighboring sites (cf number_pairs)
            # The weight has to be balanced in accordance.
            self.weight = np.exp(-self.beta_eps)

        self.weights = np.zeros(shape=(self.N))
        self.weights[0] = 1

        # :TODO: there is probably a more clever way to keep track of the sequential weights
        # probably just dividing the global weight by some configurational weight

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
        self.pos = [[0, 0, 0]]

        # Looping on the walk
        try:
            self.status = 'done'
            for step in tqdm(range(start, self.N)):
                self.update_weight()
                # :TODO: The following if condition is dumb
                if perm and step >= 3:
                    # Pruning/enriching
                    self.control_weight(step, c_m)
                self.weights[step] = self.weight

                # Stoping the walk when it reaches a closed-loop of neighbors
                if self.number_neighbors() == 0:
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
            self.status = 'killed'
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
    
    def number_pairs(self):
        '''
        This function computes the number of occupied neighbors at a given site in a history of steps. 
        It is relevent to use in the case where a polymer chain can no longer be extended (the final step is surrounded by 6 neighbors),
        or to calculate the weight.
        '''
        neighbors = self.neighborhood(self.pos[-1])
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
        W_m = c_m*self.Z[step]
        W_p = c_p*self.Z[step]

        # Pruning
        if self.weight < W_m:
            if uniform(0, 1) < 0.5:
                raise BreakException()
            else:
                self.weight *= 2

        elif self.weight > W_p:
            self.weight /= 2
            clone = copy(self)
            self.clones.append(clone)
            # print('Polymer has been cloned! (step %d)' % step)

    def update_weight(self):
        '''
        Updates weight according to the chosen random walk pattern.
        '''
        if not self.interacting:
            self.weight *= self.number_neighbors()
        else: 
            self.weight *= self.number_neighbors()*np.exp(-self.beta_eps*self.number_pairs())

    def length(self):
        '''
        Computes the squared length of a polymer (squared norm between beginning and end of said polymer).
        '''
        return np.linalg.norm(self.pos[-1]-self.pos[0], 2)**2
    
    def gyration(self):
        '''
        Computes the radius of gyration (how "curled up" the polymer is).
        '''
        N = self.pos.shape[0]                       # this in case the polymer is stuck before reaching self.N monomers
        rCM = np.sum(self.pos, axis=0)/N
        return np.sum(np.linalg.norm(self.pos - rCM, ord=2, axis=1)**2)/N
    
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
    def __init__(self, n=10, N=100, constraint='force', beta_eps=0):
        '''
        Parameters
        ----------
        n : int
          Number of monte carlo steps (number of generated polymers)
        '''
        self.n = n
        LatticePolymer.__init__(self, N, constraint, beta_eps)
        self.history = np.empty(shape=self.n, dtype=MonteCarlo) # np.full(shape=self.n, fill_value=np.nan)   # history of MC steps
        # We fill the history with nans because later we will have to compute averages on parts of the history
        # (cf estimate_Z)
        self.Z = np.empty(shape=(self.N))

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
        c_m = kwargs.get('c_m', 0.1) # lower threshold
        start = 1                    # pruning/enriching is only applied after some trials
        trial = 0
        self.clones = []
        while trial < self.n:
            print('Simulating polymer %d:' % trial)

            if trial < start or not perm:
                poly = copy(self)
                poly.gen_walk(perm=False)
                self.history[trial] = poly
            else:
                # :TODO: this is bad
                self.estimate_Z(trial)

                # Cheking if a clone has been generated for this trial
                if self.clones: # and self.history[trial-1].status == 'killed':
                    clone = self.clones[-1]
                    m = len(clone.pos)                 # number of monomers already present in present polymer
                    clone.gen_walk(m, perm, c_m)       # Processing polymer growth on top of the clone
                    self.history[trial] = clone        

                # Else generating polymer from scratch
                else:    
                    poly = copy(self)
                    poly.gen_walk(perm=perm, c_m=c_m)
                    self.history[trial] = poly

            trial += 1

    def compute_re(self):
        '''
        Computes the weighted average squared norm of a group of polymers
        '''
        # Weighted averaged length
        return np.average([self.history[trial].length() for trial in range(self.n)], \
                          weights=[self.history[trial].weight for trial in range(self.n)])
    
    def estimate_Z(self, trials):
        '''
        This function estimates a partition function for sized-L polymers (all possible Ls) with a specific number of trials.
        We only use it once.
        '''
        W = np.array([[self.history[trial].weights[L] for trial in range(trials)] for L in range(self.N)]) 
        self.Z = np.average(W, axis=1)

    # def update_Z(self, trials):
    #     self.Z = (1/trials)*((trials-1)*self.Z + np.array(self.history[trials].weights))

class BreakException(Exception):
    pass