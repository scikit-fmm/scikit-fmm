import numpy as np
from skfmm import distance
import time
import pylab as pl

def p(m): pl.matshow(m); pl.colorbar(); pl.show()
def c(m): pl.contour(m); pl.colorbar(); pl.show()
def cf(m): pl.contourf(m); pl.colorbar(); pl.show()


def test0():
    "circular level"
    N = 1500
    X, Y = np.meshgrid(np.linspace(-1,1,N),np.linspace(-1,1,N))
    r = 0.5
    dx = [2./(N-1),2./(N-1)]
    phi = (X)**2 + (Y)**2-r**2
    phi = np.ones_like(phi)
    phi[0][0]=-1
    t0=time.time()
    d = distance(phi,dx)
    t1=time.time()
    print "benchmark time", t1-t0
    return d

import cProfile
cProfile.run('test0()', 'data')

import pstats
p=pstats.Stats('data')
p.strip_dirs().sort_stats('cumulative').print_stats(10)


