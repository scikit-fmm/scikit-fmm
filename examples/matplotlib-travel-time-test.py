import skfmm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# resolution of grid for plots:
taps = 500
x_width = 2.0
y_width = 2.0

plt.figure()
X, Y = np.meshgrid(np.linspace(-0.5 * x_width, +0.5 * x_width, taps + 1), 
                   np.linspace(-0.5 * y_width, +0.5 * y_width, taps + 1))
phi = (X)**2+(Y)**2
drivers = {1: [0.5, 0]} # a dictionary with n entries
speeds = [1+X**4, 3+X**4] # a list of 2^n speed functions

print(drivers)

plt.subplot(221)
plt.title("Zero-contour of phi")
plt.contour(X, Y, phi, [1e-6], colors='black', linewidths=(3))
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(222)
plt.title("Travel time")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=x_width/taps), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(223)
plt.title("Travel time with x- \nand y- directions periodic")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=x_width/taps, periodic=True), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.subplot(224)
plt.title("Travel time with y- \ndirection periodic ")
plt.contour(X, Y, phi, [0], colors='black', linewidths=(3))
plt.contour(X, Y, skfmm.travel_time_genes(phi, drivers, speeds, dx=x_width/taps, periodic=(1,0)), 15)
plt.gca().set_aspect(1)
plt.xticks([]); plt.yticks([])

plt.savefig("testgraph.png")
