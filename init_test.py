import numpy as np
import pylab as plt

class testinit(object):
    def __init__(self, phi, h):
        self.phi = phi
        self.h = h
        assert len(phi.shape)==2

        self.find_frozen()

    def find_frozen(self):
        gr = np.gradient(phi)
        aborders = np.zeros_like(phi,dtype=bool)
        x, y = phi.shape
        for i in range(x-2):
            for j in range(y-2):
                ii=i+1
                jj=j+1
                for k in [-1,1]:
                    if phi[ii,jj] * phi[ii+k,jj] < 0: aborders[ii,jj] = True
                    elif phi[ii,jj] * phi[ii,jj+k] < 0: aborders[ii,jj] = True
        self.aborders = aborders

    def ini_frozen(self):
        for i in range(x-2):
            for j in range(y-2):
                ii=i+1
                jj=j+1
                if aborders(ii, jj):
                    pass


    def find_distance(self):
        pass

if __name__ == '__main__':
    N = 150
    x, h = np.linspace(-2.5,2.5,N,retstep=True)
    extent = (x[0],x[-1],x[0],x[-1])
    X,Y = np.meshgrid(x,x)

    # solve the eikonal equation: |grad phi| = 1
    # given elliptic zero contour 0.5*x^2 + y^2 = 1

    phi = 0.5*X**2 + Y**2 - 1.

    plt.matshow(phi)
    plt.show()

    a=testinit(phi, h)

    plt.matshow(a.aborders)
    plt.show()
