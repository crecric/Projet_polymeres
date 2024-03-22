import numpy as np
import sys
from calcul import LatticePolymer, MonteCarloFactory, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# Plot params
mpl.rcParams['figure.figsize'] = (8,7)
mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.linewidth'] = 3
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']

def make_axis(axe):
    axe.xaxis.set_minor_locator(AutoMinorLocator())
    axe.yaxis.set_minor_locator(AutoMinorLocator())
    axe.xaxis.set_tick_params(which='major', direction='in', size=10, width=2, top='on', bottom='on')
    axe.xaxis.set_tick_params(which='minor', direction='in', size=6, width=1, top='on', bottom='on')
    axe.yaxis.set_tick_params(which='major', direction='in', size=10, width=2, right='on', left='on')
    axe.yaxis.set_tick_params(which='minor', direction='in', size=6, width=1, right='on', left='on')
    pass

# Run params
arg = sys.argv[1]
if arg not in ['sarw', 'isarw', 'bisaw']:
    raise NotImplementedError("Please provide a run type in ['sarw', 'isarw', 'bisaw']")

# Params
N = 2000               # Number of monomers
n = 10000              # Number of polymers
poly_per_run = 100
runs = 100
c_m = 0.3
c_p = 3

if arg == 'sarw':
    # Generating a group of polymers with Rosenbluth method
    mcgroup = MonteCarloFactory(n=n, N=N)
    mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                        save='%s_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (arg, runs, N, poly_per_run, c_m, c_p))
    
    def r(L):
        nu = 3/5
        return L**(2*nu)
    
    ks = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1]
    ks = [int(N*k) for k in ks]
    x = np.linspace(ks[0], ks[-1], 100)
    ts = []
    es = []
    # mcgroup = MonteCarloFactory(load='%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (runs, N, poly_per_run, c_m, c_p))
    for k in ks:
        t = mcgroup.compute_observable(LatticePolymer.length, k)
        e = mcgroup.error(LatticePolymer.length, k)
        ts.append(t)
        es.append(e)
        print("re(%d) = %f +/- %f" % (k, t, e))
        print('r_theo(%d)=' % k, r(k-1))

    plt.errorbar(ks, ts, yerr=es, fmt='ro', mfc='red', mec='k', ms=6.0, capsize=3, lw=0.9, label='PERM')
    plt.plot(x, r(x), 'b-', lw=2, label='Approximate law')
    ax = plt.gca()
    make_axis(ax)
    plt.savefig('Re_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (runs, N, poly_per_run, c_m, c_p), dpi=300)
    plt.legend(loc=2, frameon=False)
    plt.ylabel(r'$\langle r_e^2\rangle$')
    plt.xlabel(r'$N$')
    plt.show()