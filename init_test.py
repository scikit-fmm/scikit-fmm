import numpy as np
import pylab as plt
from sys import float_info
from scipy.optimize import fsolve, fmin_l_bfgs_b, minimize

ainv = np.matrix([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                 [-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                 [0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0],
                 [0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0],
                 [-3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0],
                 [0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0],
                 [9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1],
                 [-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1],
                 [2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0],
                 [0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0],
                 [-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1],
                 [4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1]], dtype=np.double)

class bc_interp(object):
    def __init__(self, a):
        self.a = a.reshape((4,4),order="f")
    def __call__(self, x, y):
        a=self.a
        return  a[0,0] + a[0,1]*y + a[0,2]*y**2 + a[0,3]*y**3 + \
                a[1,0]*x + a[1,1]*x*y + a[1,2]*x*y**2 + a[1,3]*x*y**3 + \
                a[2,0]*x**2 + a[2,1]*x**2*y + a[2,2]*x**2*y**2 + \
                a[2,3]*x**2*y**3 + a[3,0]*x**3 + a[3,1]*x**3*y + \
                a[3,2]*x**3*y**2 + a[3,3]*x**3*y**3

class bc_interp_eq2(bc_interp):
    def __init__(self, interp, b0, b1):
        self.a = interp.a
        self.b0, self.b1 = b0, b1
    def __call__(self, x, y):
        a=self.a
        b0, b1 = self.b0, self.b1
        return -(b0 - x)* \
            (a[0,1] + 2*a[0,2]*y + 3*a[0,3]*y**2 + \
                           a[1,1]*x + 2*a[1,2]*x*y + \
                           3*a[1,3]*x*y**2 + a[2,1]*x**2 +\
                           2*a[2,2]*x**2*y + 3*a[2,3]*x**2*y**2 +\
                           a[3,1]*x**3 + 2*a[3,2]*x**3*y +\
                           3*a[3,3]*x**3*y**2) + (b1 - y)* \
            (a[1,0] + a[1,1]*y + a[1,2]*y**2 + a[1,3]*y**3 +\
             2*a[2,0]*x + 2*a[2,1]*x*y +\
             2*a[2,2]*x*y**2 + 2*a[2,3]*x*y**3 +\
             3*a[3,0]*x**2 + 3*a[3,1]*x**2*y +\
             3*a[3,2]*x**2*y**2 +\
             3*a[3,3]*x**2*y**3)

class testinit(object):
    def __init__(self, phi, h, X, Y):
        self.phi = phi
        self.d = np.ones_like(phi) * float_info.max
        self.d1 = np.ones_like(phi) * float_info.max
        self.h = h
        assert len(phi.shape)==2
        self.X, self.Y = X, Y

        gx, gy = np.gradient(phi)
        self.gx, self.gy = gx/2./h, gy/2./h

        xgr = np.zeros_like(phi)
        # xgr[1:-1,1:-1] = (phi[2:,2:] - phi[2:,1:-1] - phi[1:-1,2:] \
        #                   + 2*phi[1:-1,1:-1] - phi[:-2,1:-1] - phi[1:-1,:-2] \
        #                   + phi[:-2,:-2])/(2 * h**2)

        xgr[1:-1,1:-1] = (phi[2:,2:] - phi[2:,:-2] - phi[:-2,2:] + \
                          phi[:-2,:-2])/(4*h**2)
        # need to also define this for the outside edges.
        self.xgr = xgr
        self.find_frozen()
        self.ini_frozen()

    def find_frozen(self):
        aborders = np.zeros_like(phi,dtype=bool)
        border_cells = np.zeros_like(phi, dtype=bool)[:-1,:-1]
        x, y = phi.shape
        for i in range(x-2):
            for j in range(y-2):
                ii=i+1
                jj=j+1
                for k in [-1,1]:
                    if phi[ii,jj] * phi[ii+k,jj] < 0:
                        aborders[ii,jj] = True
                        border_cells[ii][jj] = True
                        border_cells[ii][jj-1] = True
                    elif phi[ii,jj] * phi[ii,jj+k] < 0:
                        aborders[ii,jj] = True
                        border_cells[ii][jj] = True
                        border_cells[ii-1][jj] = True

        self.aborders = aborders
        self.border_cells = border_cells

    def ini_frozen(self):
        x, y = self.border_cells.shape
        for i in range(x):
            for j in range(y):
                if self.border_cells[i, j]:
                    self.process_cell(i,j)


    def process_cell(self,i,j):
        # find interpolation values for this cell
        #

        coords = ((i, i+1, i, i+1), (j, j, j+1, j+1))
        #print coords
        X = np.hstack((self.phi[coords], self.gx[coords],
                       self.gy[coords], self.xgr[coords]))
        #print i,j, self.phi[i,j]
        #print X
        a =  ainv *np.matrix(X).T
        #print a
        interp = bc_interp(a)

        #print
        np.testing.assert_allclose(self.phi[i,j],     interp(0,0))
        np.testing.assert_allclose(self.phi[i+1,j],   interp(1,0))
        np.testing.assert_allclose(self.phi[i,j+1],   interp(0,1))
        np.testing.assert_allclose(self.phi[i+1,j+1], interp(1,1))

        # print interp(0,0)
        # print interp(1,0)
        # print interp(0,1)
        # print interp(1,1)

        # print "phi at cell center", interp(0.5,0.5)

        self.process_point(i,j,0,0,interp)
        self.process_point(i+1,j,1,0,interp)
        self.process_point(i,j+1,0,1,interp)
        self.process_point(i+1,j+1,1,1,interp)

    def process_point(self, i, j, ii, jj, interp):
        # ii and jj are b here in the dimensionless reference cell
        #print ii,jj
        eq2 = bc_interp_eq2(interp, ii, jj)

        def eqns(p):
            c0, c1 = p
            return interp(c0, c1)**2 + eq2(c0, c1)**2

        def eqns2(p):
            c0, c1 = p
            return (interp(c0, c1), eq2(c0, c1))


        sx,sy = fsolve(eqns2, (0.5,0.5))
        #print sx,sy, interp(sx,sy), eq2(sx,sy)

        # res = minimize(eqns, (0.5,0.5), bounds=((0,1),(0,1)),
        #                method="TNC", tol=1e-18)
        # assert res.success
        # rx, ry = res.x
        # print rx, ry, interp(rx, ry), eq2(sx,sy)
        # assert 0 <= rx <=1
        # assert 0 <= ry <=1

        if 0 <= sx <= 1 and 0 <= sy <= 1:
            dist = np.sqrt((sx-ii)**2 + (sy-jj)**2)

            if self.d[i,j] > dist:
                #print "setting d",i,j,dist
                self.d[i,j]=dist


if __name__ == '__main__':
    N = 150
    x, h = np.linspace(-2.5,2.5,N,retstep=True)
    extent = (x[0],x[-1],x[0],x[-1])
    X,Y = np.meshgrid(x,x)

    # solve the eikonal equation: |grad phi| = 1
    # given elliptic zero contour 0.5*x^2 + y^2 = 1

    phi = 0.5*X**2 + Y**2 - 1.
    #plt.matshow(phi)

    a=testinit(phi, 1, X, Y)

    # plt.matshow(a.aborders)
    # plt.matshow(a.border_cells)
    # plt.matshow(a.xgr[1:-1,1:-1])

    arr=np.copy(a.d)
    mask = arr==float_info.max
    bb=np.ma.MaskedArray(arr,mask)

    # make sure that all point adjacent to the zero contour got a value
    assert np.logical_and(a.aborders==True, np.logical_not(mask)).sum() == \
        (a.aborders==True).sum()

    plt.matshow(bb)
    plt.colorbar()

    plt.show()
