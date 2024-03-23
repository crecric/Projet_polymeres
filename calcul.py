import numpy as np
from random import choice, choices, uniform
from copy import deepcopy as copy
from tqdm import tqdm
from colorama import Fore, Style
import pickle 
from itertools import zip_longest

class LatticePolymer:
    def __init__(self, N, boltzmann_energy=1, boltzmann_force=1):
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
        # Setting Interaction and constraints
        self.q = boltzmann_energy
        self.b = boltzmann_force
        if self.q == 1:
            self.interacting = False
        else:
            self.interacting = True
        if self.b != 1:
            self.forced = True
        else:
            self.forced = False

        self.N = N

        self.n_c = 0
        self.n_p = 0

    def gen_walk(self, start=1, perm=False, c_m=0.2, c_p=2):
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
            # if self.interacting:
            #     # The first iteration will wrongly count -1 occupied neighboring sites (cf number_pairs)
            #     # The weight has to be balanced in accordance.
            #     self.weight = np.log10(self.q)
        
        self.heatup = 0
        heatup_thres = max(10, self.N // 500)
        # if self.tours.count(self.tour) >= self.relaxation and self.origin <= int(0.2*self.N):
        #     print('%sRelaxing PERM...%s' % (Fore.YELLOW, Style.RESET_ALL))
        self.c0 = len(self.clones)
        self.k_ = 0
        self.c_ = 0

        # Looping on the walk
        try:
            for step in tqdm(range(start, self.N)):
                # Stoping the walk when it reaches a closed-loop of neighbors
                if self.number_neighbors() == 0:
                    # del self.weights[step][-1]
                    break
                
                self.update_weight(step)
            
                self.heatup+=1
                # if start == 1 and not self.clones:
                #     heatup_thres = self.start_heatup
                # else:
                #     heatup_thres = max(5, self.N//500)
                if perm: # and self.heatup >= heatup_thres:
                        # input('%sAdding FORCING to polymer.%s' % (Fore.YELLOW, Style.RESET_ALL))
                        # c_m /= 10
                    # Pruning/enriching
                    self.control_weight(step, c_m, c_p)
                    
                if len(self.clones) == self.c0:
                    self.c_ += 1
                self.k_ += 1
                if len(self.clones) == self.c0:
                    self.c_ += 1
                self.k_ += 1
                    
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
        if not self.forced:
            x, y, z = choice(self.neighborhood(self.pos[-1]))
        else:
            weights = [np.sqrt(self.b), 1/np.sqrt(self.b), 1, 1, 1, 1]
            
            s = np.sum(weights)
            weights = [w/s for w in weights]
            x, y, z = choices(self.neighborhood(self.pos[-1]), weights=weights)[0]
        return x, y, z
    
    def control_weight(self, step, c_m, c_p):
        '''
        This function applies the pruning/enriching algorithm to the Rosenbluth sampling.
        '''
        # Current estimator of partition function
        W_m = np.log10(c_m)+self.Z[step]
        W_p = np.log10(c_p)+self.Z[step]

        # if self.tours.count(self.tour) >= self.relaxation and self.origin <= int(0.2*self.N):
        #     W_m = 0
        #     W_p = 10**(10000)

        self.heatup=0 
        # Pruning
        if self.weight < W_m:
            self.n_p += 1
            self.n_p += 1
            # print(W_m)
            # print(self.weights[step][-1])
            # print(self.weights[step][-2])
            self.k.append(self.k_)
            self.k_ = 0
            # print(self.zfactor)
            self.k.append(self.k_)
            self.k_ = 0
            # print(self.zfactor)
            #print(np.power(10, self.y)) 
            if uniform(0, 1) < 0.5:
                print('%sPolymer has been KILLED!%s' % (Fore.RED, Style.RESET_ALL))
                #del self.weights[step][-1]
                raise BreakException()
            else:
                print('%sPolymer has SURVIVED!%s' % (Fore.GREEN, Style.RESET_ALL))
                self.weight += np.log10(2)
                # self.weights[step][-1] += np.log10(2)
                # for s in range(1, step+1):
                #     self.weights[s][-1] += np.log10(2)

        elif self.weight > W_p and step != self.N-1: # and step <= int(0.95*self.N):
            self.n_c += 1
            self.n_c += 1
            self.weight -= np.log10(2)
            # self.weights[step][-1] -= np.log10(2)
            # for s in range(1, step+1):
            #     # self.weights[s] = [w-np.log10(2) for w in self.weights[s]]
            #     self.weights[s].append(self.weights[s][-1])
            self.clones.append(self.checkpoint())
            print('%sPolymer has been CLONED!%s' % (Fore.MAGENTA, Style.RESET_ALL))
            self.c.append(self.c_)
            self.c0 = len(self.clones)
            self.c.append(self.c_)
            self.c0 = len(self.clones)
    
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
        # Generating a new direction
        x, y, z = self.random_step()
        while [x, y, z] in self.pos:
            # Generating new step if the step is already present in the history of steps
            x, y, z = self.random_step()
        self.pos.append([x,y,z])

        if not self.forced:
            if not self.interacting:
                self.weight += np.log10(numb_neigh)
            else: 
                numb_pairs = 5 - self.number_neighbors()
                self.weight += np.log10(numb_neigh*self.q**(numb_pairs))
        else:
            numb_pairs = 5 - self.number_neighbors()
            delta_x = self.pos[-1][0] - self.pos[-2][0]
            if delta_x == 0:
                self.weight += np.log10(self.q**(numb_pairs))
            elif delta_x == 1:
                self.weight += np.log10(self.q**(numb_pairs)*np.sqrt(self.b))
            else:
                self.weight += np.log10(self.q**(numb_pairs)/np.sqrt(self.b))

        # Calculation of Z (for Monte Carlo purposes)
        # try:
        # Finding the max weight
        self.weights[step].append(self.weight)
        w_max = np.array(self.weights[step]).max()
        self.y = [x-w_max for x in self.weights[step]]

        self.zfactor = np.sum(np.power(10, self.y))
        trials = len(self.weights[step])
        self.Z[step] = np.log10(1/(trials)) + np.log10(self.zfactor) + w_max

        # except OverflowError:
        #     print('%sOVERFLOWERROR passed%s' % (Fore.RED, Style.RESET_ALL))
        #     pass

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
    def extension(pos):
        '''
        Computes the extension of a set of vector in the x direction.
        '''
        return pos[-1][0] - pos[0][0]
    
    @staticmethod
    def neighborhood(r):
        '''
        Checks the neighboring sites of a given vector r.
        '''
        x, y, z = r
        neighbors = [[(x+1), y, z], [(x-1), y, z], [x, (y+1), z], [x, (y-1), z], [x, y, (z+1)], [x, y, (z-1)]]
        return neighbors

def MonteCarloFactory(load=None, *args, **kwargs):
    '''
    This function serves as a factory to load or initiate a new MonteCarlo instance.
    '''
    if load:
        with open(load, "rb") as f:
            return pickle.load(f)
    else:
        return MonteCarlo(*args, **kwargs)
    
class MonteCarlo(LatticePolymer):
    '''
    Generates collection of polymers.
    Returns thermodynamic observables.
    '''
    def __init__(self, n=1, N=100, boltzmann_energy=1, boltzmann_force=1):
        '''
        Parameters
        ----------
        n : int
          Number of monte carlo steps (number of generated polymers)
        '''
        self.n = n
        LatticePolymer.__init__(self, N, boltzmann_energy, boltzmann_force)
        self.weights = [[] for _ in range(self.N)]
        self.Z = np.zeros(shape=self.N)
        self.history = {'weight': [], 'pos': [], 'origin': []}

    def rosenbluth(self, perm=False, save=None, **kwargs):
        '''
        Fills the history with the polymers simulated by a random walk with the normal Rosenbluth method.
        Parameters
        ----------
        perm : bool
            Tells if pruning/enriching is applied during sampling
        c_m : float
            Pruning strength (as in lower threshold = c_m * current estimator of Z)
        '''
        # Run variables
        self.clones = []
        # self.clones will eventually take the form [clone0_properties, clone1_properties, ...]

        # Run globals
        self.perm = perm
        c_m = kwargs.get('c_m', 0.2)     # lower threshold
        c_p = kwargs.get('c_p', 10*c_m)
        # self.relaxation = kwargs.get('relaxation', max(250, self.n//5))
        start = 50 # self.N//100                       # pruning/enriching is only applied after some trials

        # Run iterators
        self.trial = 0                  
        self.desired_trials = 0
        # self.start_heatup = max(5, self.N//500)
        self.cloning_freeze = 0
        self.tour = 0
        self.tours = []
        self.c = []
        self.k = []
        self.c = []
        self.k = []

        while self.desired_trials < self.n:
            print('Simulating Polymer %d / Trial %d / Tour %d' % (self.desired_trials, self.trial, self.tour))
            self.origin = 0
           # if self.cloning_freeze >= 750:
           #     break
            if (self.trial < start and self.r==1) or not self.perm:
                self.gen_walk(perm = False)
            else:
                # Cheking if a clone has been generated for this trial
                if self.clones:
                    clone = self.clones[-1]
                    m = len(clone['pos'])                                # Number of monomers already present in present polymer
                    self.origin = m
                    self.reset(clone['weight'], clone['pos'])
                    self.gen_walk(m, perm=True, c_m=c_m, c_p=c_p)        # Processing polymer growth on top of the clone
                    self.clones.remove(clone)
                    self.cloning_freeze = 0
                # Else generating polymer from scratch
                else:    
                    self.gen_walk(perm=True, c_m=c_m, c_p=c_p)
                    self.cloning_freeze += 1

            self.history['weight'].append(self.weight)
            self.history['pos'].append(self.pos) 
            self.history['origin'].append(self.origin)
            # if self.cloning_freeze >= 40:
            #     self.start_heatup += self.N//500
            #     input('%sPress SPACE to increase PERM cooldown to %d steps%s' % (Fore.CYAN, self.start_heatup, Style.RESET_ALL))

            if not self.clones:
                self.tour += 1
            self.tours.append(self.tour)
            self.trial += 1
            if self.pos.shape[0] == self.N:
                self.desired_trials += 1
        
        # Saving
        if save != None:
            self.save(save)

    def multiple_PERM(self, runs=50, poly_per_run=100, c_m=(1/3), c_p=3, save=None):
        '''
        This function runs a collection of PERM simulations to achieve tour-decorrelation.
        '''
        self.n = poly_per_run
        # Variables
        # self.weights_stack = []
        # self.history_stack = []

        # Iterators
        self.r = 1
        self.all_tours = 0
        self.all_trials = 0
        self.all_desired_trials = 0

        # Sequential decorrelated runs
        while self.r <= runs:
            print('%sStarting run %d/%d%s' % (Fore.YELLOW, self.r, runs, Style.RESET_ALL))
            self.rosenbluth(perm=True, c_m=c_m, c_p=c_p, save=None)
            # self.weights_stack.append([list(i) for i in zip_longest(*self.weights)])
            # self.history_stack.append(self.history)

            self.all_tours += self.tour
            self.all_desired_trials += self.desired_trials
            self.all_trials += self.trial
            self.r += 1
            print('%sSimulated so far %d polymers / %d trials / %d tours%s' % (Fore.YELLOW, self.all_desired_trials, \
                                                                               self.all_trials, self.all_tours, Style.RESET_ALL))
            
        # Rebuilding stacked MC attributes

        # all_origins = [x['origin'] for x in self.history_stack]
        # all_positions = [x['pos'] for x in self.history_stack]

        # self.history['origin'] = sum(all_origins, [])
        # self.history['pos'] = sum(all_positions, [])
        # # self.weights = list(np.concatenate(self.weights_stack).T)
        # transposed_weights = sum(self.weights_stack, [])
        # self.weights = [list(i) for i in zip_longest(*transposed_weights)]

        # # Computing whole Z
        # self.Z = np.zeros(shape=self.N)

        # # Finding the max weight
        # for step in range(1, self.N):
        #     logweights = [x for x in self.weights[step] if x is not None]
        #     w_max = np.array(logweights).max()
        #     self.y = [x-w_max for x in logweights]

        #     self.zfactor = np.sum(np.power(10, self.y))
        #     trials = len(self.weights[step])
        #     self.Z[step] = np.log10(1/(trials)) + np.log10(self.zfactor) + w_max

        # Saving    
        if save != None:
            self.save(save)

    def compute_observable(self, obs, N):
        '''
        Computes an observable average at a given number of monomers.
        '''
        logweights = [x for x in self.weights[N-1] if x is not None]
        trials = len(logweights)
        positions = [pos[:N] for i, pos in enumerate(self.history['pos']) if pos.shape[0] >= N and self.history['origin'][i] < N]

        observables = [obs(pos) for pos in positions]
        weights = np.power(10, [w-np.log10(trials)-self.Z[N-1] for w in logweights])
        
        return np.average(observables, weights=weights)

    def error(self, obs, N):
        '''
        Estimate the Monte Carlo error using approximate analytical formula and boostrapping
        '''
        # Approximate analytical formula
        # err_a = np.sqrt((n/(n-1))*(np.sum([self.history[trial].weight**2*(self.history[trial].length()**2- np.sum([self.history[j].length() for j in range(self.n)]))**2) for trial in range(self.n)]/(w**2))

        logweights = [x for x in self.weights[N-1] if x is not None]
        trials = len(logweights)
        positions = [pos[:N] for i, pos in enumerate(self.history['pos']) if pos.shape[0] >= N and self.history['origin'][i] < N]

        observables = np.array([obs(pos) for pos in positions])
        weights = np.power(10, [w-np.log10(trials)-self.Z[N-1] for w in logweights])

        # Boostrapping
        num_bootstrap_samples = trials//10
    
         # Table of zero
        bootstrap_obs = np.zeros(num_bootstrap_samples)
    
        for i in range(num_bootstrap_samples):
            bootstrap_indices = np.random.choice(trials, size=trials, replace=True)
            bootstrap_obs[i] = np.average(observables[bootstrap_indices], weights=weights[bootstrap_indices])

         # standard error 
        error = np.std(bootstrap_obs)
        return error
    
    def save(self, file):
        with open(file, "wb") as f:
            pickle.dump(self, f)

class BreakException(Exception):
    pass