import numpy as np
import pylab as pl
import skfmm

X, Y = np.meshgrid(np.linspace(-1,1,501), np.linspace(-1,1,501))
phi = np.ones_like(X)
phi[(X+0.8)**2+(Y+0.8)**2<0.01] = -1
speed = 1+X**2+Y**2

pl.subplot(221)
pl.title("Zero-contour of phi")
pl.contour(X, Y, phi, [0], colors='black', linewidths=(3))
pl.gca().set_aspect(1)

pl.subplot(222)
pl.title("Distance")
pl.contour(X, Y, phi, [0], colors='black', linewidths=(3))
pl.contour(X, Y, skfmm.distance(phi, dx=2.0/500), 15)
pl.colorbar()
pl.gca().set_aspect(1)

pl.subplot(223)
pl.title("Distance with PBC")
pl.contour(X, Y, phi, [0], colors='black', linewidths=(3))
pl.contour(X, Y, skfmm.distance(phi, dx=2.0/500, periodic=True), 15)
pl.colorbar()
pl.gca().set_aspect(1)

pl.subplot(224)
pl.title("Travel time with PBC in y-direction")
pl.contour(X, Y, phi, [0], colors='black', linewidths=(3))
pl.contour(X, Y, skfmm.travel_time(phi, speed, dx=2.0/500, periodic=(1,0)), 15)
pl.colorbar()
pl.gca().set_aspect(1)

pl.show()
