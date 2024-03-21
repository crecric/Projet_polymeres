import numpy as np
from random import choice, uniform
from copy import deepcopy as copy
from tqdm import tqdm
from colorama import Fore, Style
import pickle 

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
            if self.interacting:
                # The first iteration will wrongly count -1 occupied neighboring sites (cf number_pairs)
                # The weight has to be balanced in accordance.
                self.weight = np.log10(np.exp(-self.beta_eps))
        
        self.heatup = 0
        heatup_thres = self.N // 60
        if self.tours.count(self.tour) >= self.relaxation:
            print('%sRelaxing PERM...%s' % (Fore.YELLOW, Style.RESET_ALL))
        # self.c0 = len(self.clones)
        # self.k_ = 0
        # self.c_ = 0

        # Looping on the walk
        try:
            for step in tqdm(range(start, self.N)):
                # Stoping the walk when it reaches a closed-loop of neighbors
                if self.number_neighbors() == 0:
                    # del self.weights[step][-1]
                    break

                self.update_weight(step)
                # Generating a new direction
                x, y, z = self.random_step()
                while [x, y, z] in self.pos:
                    # Generating new step if the step is already present in the history of steps
                    x, y, z = self.random_step()

                self.pos.append([x,y,z])
            
                self.heatup+=1
                # if start == 1 and not self.clones:
                #     heatup_thres = self.start_heatup
                # else:
                #     heatup_thres = max(5, self.N//500)
                if perm :#and self.heatup >= heatup_thres:
                        # input('%sAdding FORCING to polymer.%s' % (Fore.YELLOW, Style.RESET_ALL))
                        # c_m /= 10
                    # Pruning/enriching
                    self.control_weight(step, c_m, c_p)
                    
                # if len(self.clones) == self.c0:
                #     self.c_ += 1
                # self.k_ += 1
                    
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
    
    def control_weight(self, step, c_m, c_p):
        '''
        This function applies the pruning/enriching algorithm to the Rosenbluth sampling.
        '''
        # Current estimator of partition function
        W_m = np.log10(c_m)+self.Z[step]
        W_p = np.log10(c_p)+self.Z[step]

        if self.tours.count(self.tour) >= self.relaxation:
            W_m = 0
            W_p = 10**(10000)
        # if self.weight<self.Z[step]:
        #     print('Z=',self.Z[step])
        #     print("W_M=", W_m)
        #     print("Weight=", self.weight)
        self.heatup=0 
        # Pruning
        if self.weight < W_m:
            # print(W_m)
            #print(self.weights[step])
            # print(self.weights[step][-2])
            # self.k.append(self.k_)
            # self.k_ = 0
            
            
            
            
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
            self.weight -= np.log10(2)
            # self.weights[step][-1] -= np.log10(2)
            # for s in range(1, step+1):
            #     # self.weights[s] = [w-np.log10(2) for w in self.weights[s]]
            #     self.weights[s].append(self.weights[s][-1])
            self.clones.append(self.checkpoint())
            print('%sPolymer has been CLONED!%s' % (Fore.MAGENTA, Style.RESET_ALL))
            # self.c.append(self.c_)
            # self.c0 = len(self.clones)
    
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
            self.y = [x-w_max for x in self.weights[step]]
            

            self.zfactor = sum(np.power(10, self.y))
            trials = len(self.weights[step])
            self.Z[step] = np.log10(self.zfactor) + w_max + np.log10(1/(trials))

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
    def __init__(self, n=1, N=100, constraint='force', beta_eps=0, load=''):
        '''
        Parameters
        ----------
        n : int
          Number of monte carlo steps (number of generated polymers)
        '''
        self.n = n
        LatticePolymer.__init__(self, N, constraint, beta_eps)
        self.weights = [[] for _ in range(self.N)]
        self.Z = np.zeros(shape=self.N)
        self.history = {'weight': [], 'pos': [], 'origin': []}
        self.clones = []
        # self.clones will eventually take the form [clone0_properties, clone1_properties, ...]

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
        self.perm = perm
        c_m = kwargs.get('c_m', 0.2)    # lower threshold
        c_p = kwargs.get('c_p', 20*c_m)
        self.relaxation = kwargs.get('relaxation', max(250, self.n//5))
        start = 10                       # pruning/enriching is only applied after some trials
        self.trial = 0                  
        self.desired_trials = 0
        # self.start_heatup = max(5, self.N//500)
        # self.cloning_freeze = 0
        self.tour = 0
        self.tours = []
        # self.c = []
        # self.k = []

        while self.desired_trials < self.n:
            print('Simulating Polymer %d / Trial %d / Tour %d' % (self.desired_trials, self.trial, self.tour))
            self.origin = 0
            if self.trial < start or not self.perm:
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
                    # self.cloning_freeze = 0
                # Else generating polymer from scratch
                else:    
                    self.gen_walk(perm=True, c_m=c_m, c_p=c_p)
                    # self.cloning_freeze += 1

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

    def compute_observable(self, obs, N):
        '''
        Computes an observable average at a given number of monomers.
        '''
        logweights = self.weights[N-1]
        trials = len(logweights)
        positions = [pos[:N] for i, pos in enumerate(self.history['pos']) if pos.shape[0] >= N and self.history['origin'][i] < N]

        observables = [obs(pos) for pos in positions]
        weights = np.power(10, [w-np.log10(trials)-self.Z[N-1] for w in logweights])
        print(trials)
        print(len(positions))
        return np.average(observables, weights=weights)

    def save(self, file):
        with open(file, "wb") as f:
            pickle.dump(self, f)

class BreakException(Exception):
    pass