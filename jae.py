import numpy as np
import matplotlib.pyplot as plt
import skfmm

# grid info
N = 1000
x, h = np.linspace(-2.5,2.5,N,retstep=True)
extent = (x[0],x[-1],x[0],x[-1])
X,Y = np.meshgrid(x,x)

# solve the eikonal equation: |grad phi| = 1
# given elliptic zero contour 0.5*x^2 + y^2 = 1
phi = np.ones(X.shape)
phi[0.5*X**2 + Y**2 < 1.] = -1.
phi = skfmm.distance(phi, h, order=2)

# print max/mean of Q = ||grad phi| - 1 | in a band around the zero-contour
# gradient is approximated by central difference
normgrad = np.zeros_like(phi)
normgrad[1:-1,1:-1] = np.sqrt( (phi[1:-1,2:]-phi[1:-1,:-2])**2 +
                               (phi[2:,1:-1]-phi[:-2,1:-1])**2 )/(2.*h)
band =  (0.5*X**2 + Y**2 < 1.1)*~(0.5*X**2 + Y**2 < .9)
Q = np.abs(normgrad[band]-1)
print('evaluating\tQ = ||grad phi| - 1|\tnear zero-contour:')
print 'max(Q) =', np.max(Q) #this gives a value around 0.1450 independent of N
print 'mean(Q) =', np.mean(Q)
print 'grid spaceing =', h
print 'grid shape =', phi.shape

# laplacian approximated by 9pt stencil
# I need this laplacian to approximate the curvature of the zero-countour
lapl = np.zeros_like(phi)
lapl[1:-1,1:-1] = (phi[ :-2, :-2] + phi[2:  ,2:  ] + phi[2:  , :-2] +
                   phi[ :-2,2:  ] + phi[1:-1, :-2] + phi[1:-1,2:  ] +
                   phi[ :-2,1:-1] + phi[2:  ,1:-1] - 8.*phi[1:-1,1:-1])/h**2

# contour plot of phi and grad^2 phi
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.imshow(phi,extent=extent,origin='lower',alpha=0.5)
cs1 = ax.contour(X,Y,phi,[-.5, -.25, 0., .25, .5, .75], alpha=1, colors='k')
cs2 = ax.contour(X,Y,lapl, alpha=0.5)
ax.clabel(cs1, inline=1, fontsize=10)
plt.savefig('example.png')
plt.show()
