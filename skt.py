import numpy as np
import pylab as plt

from skfmm import distance

"""

"""

for N in [21]:#, 41, 81, 161]:

    X, Y = np.meshgrid(np.linspace(-1,1,N), np.linspace(-1,1,N))
    dx = 2.0/(N-1)
    r     = 0.5
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2


    # plt.matshow(phi)
    # plt.colorbar()
    # plt.show()

    #d1 = distance(phi,dx=dx,order=1)
    d2 = distance(phi,dx=dx,initorder=2)
    #plt.matshow(d2)
    #plt.colorbar()
    #plt.show()
    #import IPython; IPython.embed()

    #print "order 1", N, d1[0,0] - (np.sqrt(2.0)-0.5)
    print "order 2", N, d2[0,0] - (np.sqrt(2.0)-0.5)

# from  master 78a243d
# order 2 21 0.00721008286517
# order 2 41 0.000794988804239
# order 2 81 -0.000549263254452
# order 2 161 -7.52168527809e-05

# this seems to choke for small values of phi?
# from bicubic-init fcf3315
# order 2 21 -0.00517334060839
# order 2 41 -0.000233720147591
# order 2 81 -0.0389212541035
# order 2 161 -0.0150101932032
