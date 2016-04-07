import numpy as np
import pylab as plt

from skfmm import distance

"""

"""

for N in [21, 41, 81, 161]:

    X, Y = np.meshgrid(np.linspace(-1,1,N), np.linspace(-1,1,N))
    dx = 2.0/(N-1)
    r     = 0.5
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    exact = np.sqrt(X**2+Y**2) - r

    d1 = distance(phi,dx=dx)
    d2 = distance(phi,dx=dx,initorder=2)


    plt.contour(X,Y,d1,[0.4],colors=["red"])
    plt.contour(X,Y,d2,[0.4],colors=["black"])
    plt.contour(X,Y,exact,[0.4],colors=["green"])
    plt.gca().set_aspect(1)
    plt.show()

    #import IPython; IPython.embed()


    print "linear L2 Norm", abs((d1 - exact)**2).sum()
    print "bi-cubic L2 Norm", abs((d2 - exact)**2).sum()
    print "linear L infinity", abs(d1 - exact).max()
    print "bi-cubic L infinity", abs(d2 - exact).max()
    print "central point error order 1", N, d1[0,0] - (np.sqrt(2.0)-0.5)
    print "central point error order 2", N, d2[0,0] - (np.sqrt(2.0)-0.5)


    print
    plt.subplot(121)
    plt.imshow(exact-d1)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(exact-d2)
    plt.colorbar()
    plt.show()
