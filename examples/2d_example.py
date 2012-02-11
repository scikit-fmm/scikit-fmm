import numpy as np
import pylab as pl
import skfmm

X, Y = np.meshgrid(np.linspace(-1,1,200), np.linspace(-1,1,200))
phi = -1*np.ones_like(X)

phi[X>-0.5] = 1
phi[np.logical_and(np.abs(Y)<0.25, X>-0.75)] = 1

pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.title('Boundary location: the zero contour of phi')
pl.savefig('2d_phi.png')
pl.show()

d = skfmm.distance(phi, dx=1e-2)
pl.title('Distance from the boundary')
pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.contour(X, Y, d, 15)
pl.colorbar()
pl.savefig('2d_phi_distance.png')
pl.show()

speed = np.ones_like(X)
speed[Y>0] = 1.5
t = skfmm.travel_time(phi, speed, dx=1e-2)

pl.title('Travel time from the boundary')
pl.contour(X, Y, phi,[0], linewidths=(3), colors='black')
pl.contour(X, Y, t, 15)
pl.colorbar()
pl.savefig('2d_phi_travel_time.png')
pl.show()

mask = np.logical_and(abs(X)<0.1, abs(Y)<0.5)
phi  = np.ma.MaskedArray(phi, mask)
t    = skfmm.travel_time(phi, speed, dx=1e-2)
pl.title('Travel time from the boundary with an obstacle')
pl.contour(X, Y, phi, [0], linewidths=(3), colors='black')
pl.contour(X, Y, phi.mask, [0], linewidths=(3), colors='red')
pl.contour(X, Y, t, 15)
pl.colorbar()
pl.savefig('2d_phi_travel_time_mask.png')
pl.show()
