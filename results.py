import numpy as np
import sys
from calcul import LatticePolymer, MonteCarloFactory, MonteCarlo
import visualisation
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import time

start_time = time.time()

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
if arg not in ['sarw', 'isaw', 'bisaw']:
    raise NotImplementedError("Please provide a run type in ['sarw', 'isarw', 'bisaw']")

# Params
N = 2000               # Number of monomers
n = 100              # Number of polymers
poly_per_run = 500
runs = 3
c_m = 0.8
c_p = 15.0
log = True

if arg == 'sarw':
    # Generating a group of polymers with Rosenbluth method
    mcgroup = MonteCarloFactory(n=n, N=N)
    mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                        save='%s_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (arg, runs, N, poly_per_run, c_m, c_p))
    
    def r(L):
        nu = 0.586
        return L**(2*nu)
    
    if log ==True:
        ks=np.logspace(1, np.log10(N), num=15)  
        ks = [int(k) for k in ks]
    else :
        ks= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1]
        ks = [int(N*k) for k in ks]

    x = np.linspace(ks[0], ks[-1], 100)
    ts = []
    es = []
    #mcgroup = MonteCarloFactory(load='sarw_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % (runs, N, poly_per_run, c_m, c_p))
    for k in ks:
        t = mcgroup.compute_observable(LatticePolymer.length, k)
        e = mcgroup.error(LatticePolymer.length, k)
        if log ==True :
            ts.append(t/r(k-1)) 
            es.append(e/r(k-1))
        else :
            ts.append(t)
            es.append(e)
        print("re(%d) = %f +/- %f" % (k, t, e))
        print('r_theo(%d)=' % k, r(k-1))

    plt.errorbar(ks, ts, yerr=es, fmt='ro', mfc='red', mec='k', ms=6.0, capsize=3, lw=0.9, label='PERM')
    if not log : plt.plot(x, r(x), 'b-', lw=2, label='Approximate law')
    if log==True : plt.xscale('log')
    ax = plt.gca()
    make_axis(ax)
    plt.legend(loc=2, frameon=False)
    plt.ylabel(r'$\langle r_e^2\rangle$')
    plt.xlabel(r'$N$')
    plt.savefig('Re_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (runs, N, poly_per_run, c_m, c_p), dpi=300)

    end_time = time.time()
    plt.show()



elif arg == 'isaw':
    energy = [1.3, 1.305, 1.3087, 1.310, 1.315]
    ks = np.logspace(1, np.log10(N), num=15)
    ks = [int(k) for k in ks]
    fmts = ['r+-', 'gs-', 'bx-', 'm^-', 'cD-']
    for i, en in enumerate(energy):
        mcgroup = MonteCarloFactory(n=n, N=N, boltzmann_energy=en)
        mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                            save='%s_%.3fq_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, en, runs, N, poly_per_run, c_m, c_p))
        
        ts = np.empty(shape=(len(energy), len(ks)))
        es = np.empty(shape=(len(energy), len(ks)))
        mcgroup = MonteCarloFactory(load='%s_%.3fq_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, en, runs, N, poly_per_run, c_m, c_p))
        for j, k in enumerate(ks):
            t = mcgroup.compute_observable(LatticePolymer.length, k)
            e = mcgroup.error(LatticePolymer.length, k)
            ts[i,j] = (t/k)
            es[i,j] = (e/k)

        plt.plot(ks, ts[i, :], fmts[i], mfc='none', ms=8.0, lw=0.5, label='q=%.4f' % en)
        
    ax = plt.gca()
    make_axis(ax)
    plt.xscale('log')
    plt.legend(loc=2)
    plt.ylabel(r'$\frac{1}{N}\langle r_e^2\rangle$')
    plt.xlabel(r'$N$')
    plt.savefig('Re_%s_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (arg, runs, N, poly_per_run, c_m, c_p), dpi=300)

    end_time = time.time()
    plt.show()

else:
    energy = 1.5
    force = [1, 1.45, 1.55, 1.6, 1.65, 1.70]
    ks = np.linspace(2, N, num=50)
    ks = [int(k) for k in ks]
    fmts = ['r-', 'g.', 'c--', 'm-', 'y.', 'b-']
    for i, f in enumerate(force):
        mcgroup = MonteCarloFactory(n=n, N=N, boltzmann_energy=energy, boltzmann_force=f)
        mcgroup.multiple_PERM(runs=runs, poly_per_run=poly_per_run, c_m=c_m, c_p=c_p, \
                            save='%s_%.2fb_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, f, runs, N, poly_per_run, c_m, c_p))
        
        ts = np.empty(shape=(len(force), len(ks)))
        es = np.empty(shape=(len(force), len(ks)))
        mcgroup = MonteCarloFactory(load='%s_%.2fb_%druns_%dmonom_%dpoly_%.2f_%.2f.pkl' % \
                                (arg, f, runs, N, poly_per_run, c_m, c_p))
        for j, k in enumerate(ks):
           
            t = mcgroup.compute_observable(LatticePolymer.extension, k)
            e = mcgroup.error(LatticePolymer.extension, k)
            ts[i,j] = t
            es[i,j] = e

        plt.plot(ks, ts[i, :], fmts[i], lw=0.5, label='b=%.3f' % f)
        
    
    ax = plt.gca()
    make_axis(ax)
    plt.legend(loc=1)
    plt.ylabel(r'$\langle x\rangle$')
    plt.xlabel(r'$N$')
    plt.savefig('X_%s_%druns_%dmonom_%dpoly_%.2f_%.2f.jpg' % (arg, runs, N, poly_per_run, c_m, c_p), dpi=300)
    end_time = time.time()
    plt.show()



#Heatmap
#visualisation.polyCloud3D(mcgroup,N,n)

#elapsed_time = end_time - start_time

#print("Elapsed time:", elapsed_time, "seconds")