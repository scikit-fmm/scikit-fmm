import numpy as np
from skfmm import extension_velocities

N     = 150
X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
r     = 1.75
dx    = 2.0 / (N - 1)
phi   = (X) ** 2 + (Y+1.85) ** 2 - r ** 2
speed = X+1.25
d, f_ext = extension_velocities(phi, speed, dx)

import pylab as pl

pl.subplot(131)
pl.title("Zero-contour of phi")
pl.contour(X, Y, phi,[0], colors='black', linewidths=(3))
pl.gca().set_aspect(1)

pl.subplot(132)
pl.title("Interface velocity")
pl.contour(X, Y, phi,[0], colors='black', linewidths=(3))
pl.contourf(X, Y, speed)
pl.gca().set_aspect(1)

pl.subplot(133)
pl.title("Extension velocities")
pl.contour(X, Y, phi,[0], colors='black', linewidths=(3))
pl.contourf(X, Y, f_ext)
pl.gca().set_aspect(1)
#pl.colorbar()

pl.show()

