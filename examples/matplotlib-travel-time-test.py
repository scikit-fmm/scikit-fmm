import skfmm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.figure()
X, Y = np.meshgrid(np.linspace(-1,1,501), np.linspace(-1,1,501))
phi = (X)**2+(Y)**2
speeds = 1+X**4

# TODO speeds as array of speed arrays? or list

drivers = phi * 0

print(drivers)

plt.subplot(221)
plt.title("Zero-contour of phi")
plt.contour(X, Y, phi, [1e-6], colors='black', linewidths=(3))
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(222)
plt.title("Travel time")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=2.0/500), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(223)
plt.title("Travel time with x- \nand y- directions periodic")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=2.0/500, periodic=True), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(224)
plt.title("Travel time with y- \ndirection periodic ")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=2.0/500, periodic=(1,0)), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.savefig("testgraph.png")
