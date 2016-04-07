import numpy as np
import pylab as plt

from skfmm import distance

x, h= np.linspace(0,5,100,retstep=True)

X,Y = np.meshgrid(x,x)

xc,yc = -0.5e-1, -0.5e-1
#xc,yc = 0,0
phi = (X - xc)**2 + (Y-yc)**2 - 3**2
exact = np.sqrt((X - xc)**2 + (Y-yc)**2) - 3

d1 = distance(phi,initorder=1)
d2 = distance(phi,initorder=2)

ng1 = np.linalg.norm(np.gradient(d1, edge_order=2),axis=0)
ng2 = np.linalg.norm(np.gradient(d2, edge_order=2),axis=0)


plt.subplot(211)
plt.title("Distance function gradient magnitude\nFirst order init")
plt.imshow(ng1)
plt.colorbar()
plt.subplot(212)
plt.title("Distance function gradient magnitude\nSecond order init")
plt.imshow(ng2)
plt.colorbar()
plt.show()
