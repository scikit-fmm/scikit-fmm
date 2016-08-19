# this is a test that we are calculating the cross derivatives correctly

import numpy as np
import pylab as plt

x, h= np.linspace(0,5,100,retstep=True)

X,Y = np.meshgrid(x,x)

phi = (X-1.5)**2*(Y-2.1)**3
# from sympy.abc import x, y
# from sympy import diff
# expr = (x-1.5)**2 * (y-2.1)**3
# print diff(diff(expr,x),y)
dphidxdx = 3*(2*X - 3.0)*(Y - 2.1)**2 # exact cross derivative

# calculate dxdy

xgr = np.zeros_like(phi)

denom = (4*h**2)
xgr[1:-1,1:-1] = (phi[2:,2:] - phi[2:,:-2] - phi[:-2,2:] + phi[:-2,:-2])

# there are 8 special cases, 4 corners and 4 sides
#plt.matshow(phi); plt.colorbar(); plt.show()


nx, ny = phi.shape
# The cross derivative stencil looks like this:
#   C---.---A
#   |       |     phi_xy = (phi_A - phi_B - phi_C + phi_D) / (4*h**2)
#   .   .   .
#   |       |
#   D---.---B
#                A            B           C           D
xgr[   0,   0] = -(phi[ 2, 0] - phi[ 2, 2] - phi[ 0, 0] + phi[ 0, 2])
xgr[nx-1,ny-1] = -(phi[-1,-3] - phi[-1,-1] - phi[-3,-3] + phi[-3,-1])
xgr[   0,ny-1] =  (phi[ 0,-3] - phi[ 0,-1] - phi[ 2,-3] + phi[ 2,-1])
xgr[nx-1,   0] = -(phi[-1, 0] - phi[-1, 2] - phi[-3, 0] + phi[-3, 2])

# left side
xgr[0,1:-1] = -(phi[ 2, :-2] - phi[ 2, 2:] - phi[0, :-2] + phi[0, 2:])
# right side
xgr[nx-1,1:-1] =  (phi[ -3, :-2] - phi[ -3, 2:] - phi[-1, :-2] + phi[-1, 2:])
# top row
xgr[1:-1, 0] =   -(phi[ 2:,  0] - phi[2:, 2] - phi[:-2,   0] + phi[:-2, 2])
# bottom row
xgr[1:-1,ny-1] = -(phi[ 2:, -3] - phi[2:, -1] - phi[:-2, -3] + phi[:-2, -1])


xgr /= denom

print xgr
plt.matshow(xgr.T); plt.colorbar(); plt.show()
plt.matshow(dphidxdx.T); plt.colorbar(); plt.show()

eps = (dphidxdx-xgr)/(dphidxdx.max()-dphidxdx.min())
plt.matshow(eps.T); plt.colorbar(); plt.show()
# close enough?
