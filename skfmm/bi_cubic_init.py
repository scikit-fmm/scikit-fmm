import numpy as np
import pylab as plt
from sys import float_info
from scipy.optimize import fsolve

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



class BiCubicInit(object):
    def __init__(self, phi, h):
        self.phi = phi
        self.d = np.ones_like(phi) * float_info.max
        self.h = h
        self.pdict = {}
        assert len(phi.shape)==2

        gx, gy = np.gradient(phi, edge_order=2)
        self.gx, self.gy = gx/h, gy/h # is this correct for h != 1?

        xgr = np.zeros_like(phi)

        denom = (4*h**2)
        xgr[1:-1,1:-1] = (phi[2:,2:] - phi[2:,:-2] - phi[:-2,2:] + phi[:-2,:-2])
        # there are 8 special cases, 4 corners and 4 sides
        # The cross derivative stencil looks like this:
        #   C---.---A
        #   |       |     phi_xy = (phi_A - phi_B - phi_C + phi_D) / (4*h**2)
        #   .   .   .
        #   |       |
        #   D---.---B
        #                A            B           C           D
        nx, ny = phi.shape
        xgr[   0,   0] = -(phi[ 2, 0] - phi[ 2, 2] - phi[ 0, 0] + phi[ 0, 2])
        xgr[nx-1,ny-1] = -(phi[-1,-3] - phi[-1,-1] - phi[-3,-3] + phi[-3,-1])
        xgr[   0,ny-1] =  (phi[ 0,-3] - phi[ 0,-1] - phi[ 2,-3] + phi[ 2,-1])
        xgr[nx-1,   0] = -(phi[-1, 0] - phi[-1, 2] - phi[-3, 0] + phi[-3, 2])
        xgr[0,1:-1] = -(phi[ 2, :-2] - phi[ 2, 2:] - phi[0, :-2] + phi[0, 2:])
        xgr[nx-1,1:-1] =  (phi[ -3, :-2] - phi[ -3, 2:] - phi[-1, :-2] + phi[-1, 2:])
        xgr[1:-1, 0] =   -(phi[ 2:,  0] - phi[2:, 2] - phi[:-2,   0] + phi[:-2, 2])
        xgr[1:-1,ny-1] = -(phi[ 2:, -3] - phi[2:, -1] - phi[:-2, -3] + phi[:-2, -1])
        xgr /= denom

        self.xgr = xgr
        self.find_frozen()
        self.ini_frozen()

        mask = self.d==float_info.max
        if not np.logical_and(self.aborders==True, np.logical_not(mask)).sum() == \
            (self.aborders==True).sum():
            raise RuntimeError("scikit-fmm (skfmm) error, some initially \
            frozen points were not initilized correctly.")


    def find_frozen(self):
        phi = self.phi
        aborders = np.zeros_like(phi,dtype=bool)
        aborders[phi==0.0] = True
        self.d[phi==0.0] = 0.0
        border_cells = np.zeros_like(phi, dtype=bool)[:-1,:-1]
        nx, ny = phi.shape
        for i in range(nx):
            for j in range(ny):
                for k in [-1,1]: # each direction
                    # logic to not read off the edge of the arrays
                    # skip this test if:
                    ## i==0 and k == -1
                    ## i==nx-1 and k==1
                    ok_first =  not ((i==0 and k==-1) or (i==(nx-1) and k==1))
                    ok_second =  not ((j==0 and k==-1) or (j==(ny-1) and k==1))
                    if ok_first and phi[i, j] * phi[i+k,j] < 0:
                        aborders[i,j] = True
                        border_cells[i][j] = True
                        if j > 0:
                            border_cells[i][j-1] = True
                    elif ok_second and phi[i,j] * phi[i,j+k] < 0:
                        aborders[i,j] = True
                        border_cells[i][j] = True
                        if j > 0:
                            border_cells[i][j-1] = True
                        if i > 0:
                            border_cells[i-1][j] = True

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
        #if not self.aborders[i,j]: return
        coords = ((i, i+1, i, i+1), (j, j, j+1, j+1))
        X = np.hstack((self.phi[coords], self.gx[coords],
                       self.gy[coords], self.xgr[coords]))
        a =  ainv *np.matrix(X).T
        interp = bc_interp(a)

        np.testing.assert_almost_equal(self.phi[i,j],     interp(0,0))
        np.testing.assert_almost_equal(self.phi[i+1,j],   interp(1,0))
        np.testing.assert_almost_equal(self.phi[i,j+1],   interp(0,1))
        np.testing.assert_almost_equal(self.phi[i+1,j+1], interp(1,1))

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
        nx, ny = self.phi.shape
        if not self.aborders[i,j]: return
        if abs(self.phi[i,j]) < float_info.epsilon:
            self.d[i,j] = 0.0
            return

        eq2 = bc_interp_eq2(interp, ii, jj)

        def eqns(p):
            c0, c1 = p
            # try scaling these?
            return (interp(c0, c1), eq2(c0, c1))

        fprime = None
        # maybe try using the bilnear approximation as a starting value?
        # or try specifying the Jacobian of the functions
        # http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
        sol, info, ier, mesg = fsolve(eqns, (0.5,0.5),
                                      full_output=1, fprime=fprime)
        sx,sy = sol

        if ier==1:
            if 0 <= sx <= 1 and 0 <= sy <= 1:
                #print i+sx,j+sy
                #print "test", eqns(sol)
                dist = np.sqrt((sx-ii)**2 + (sy-jj)**2)
                if self.d[i,j] > dist:
                    self.d[i,j]=dist
                    self.pdict[(i,j)] = (sx-ii,sy-jj)
            else:
                # for boundary points we need to accept points outside [0,1]
                # this could/should be done more cleverly/safely
                # which direction to let sx or sy exceed the interval should be known?
                if i==0 or j==0 or i==nx-1 or j==ny-1:
                    dist = np.sqrt((sx-ii)**2 + (sy-jj)**2)
                    if self.d[i,j] > dist:
                        self.d[i,j]=dist
                        self.pdict[(i,j)] = (sx-ii,sy-jj)

                # print "-"*20
                # print sx, sy
                # print "point {} cell {}".format((i, j), (ii, jj))
                # dist = np.sqrt((sx-ii)**2 + (sy-jj)**2)
                # print "distance {}".format(dist)
                # print eqns(sol)
                # print
        else:
            #pass
            print ier, mesg
            print sx,sy, interp(sx,sy), eq2(sx,sy)


if __name__ == '__main__':
    N = 100
    x, h = np.linspace(-2.5,2.5,N,retstep=True)
    extent = (x[0],x[-1],x[0],x[-1])
    X,Y = np.meshgrid(x,x)

    # solve the eikonal equation: |grad phi| = 1
    # given elliptic zero contour 0.5*x^2 + y^2 = 1

    phi = X**2 + Y**2 - 1.
    plt.matshow(phi)
    plt.colorbar()
    plt.show()

    a=BiCubicInit(phi, 1)

    plt.matshow(a.aborders)
    plt.colorbar()
    plt.show()

    #plt.matshow(a.border_cells)
    #plt.matshow(a.xgr[1:-1,1:-1])

    arr=np.copy(a.d)
    mask = arr==float_info.max
    bb=np.ma.MaskedArray(arr,mask)
    bb[phi<0] *= -1

    # make sure that all point adjacent to the zero contour got a value
    assert np.logical_and(a.aborders==True, np.logical_not(mask)).sum() == \
        (a.aborders==True).sum()

    plt.matshow(bb)
    plt.colorbar()
    plt.show()
