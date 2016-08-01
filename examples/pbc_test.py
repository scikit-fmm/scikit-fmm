import numpy as np
import pylab as pl
import skfmm

X, Y = np.meshgrid(np.linspace(-2,2,1000), np.linspace(-2,2,1000))

# Test 0: periodic=(0,0), periodic=False, and periodic option absent should give same result
phi = -1*np.ones_like(X); phi[X**2+(Y-0.9)**2<0.5] = 1.0
speed = np.ones_like(X); speed[(X-0.9)**2+Y**2<1.0] = 2.0
if np.allclose(skfmm.distance(phi),skfmm.distance(phi,periodic=False)) and np.allclose(skfmm.distance(phi),skfmm.distance(phi,periodic=(0,0))):
    print "Test 0/distance passed"
if np.allclose(skfmm.travel_time(phi,speed),skfmm.travel_time(phi,speed,periodic=False)) and np.allclose(skfmm.travel_time(phi,speed),skfmm.travel_time(phi,speed,periodic=(0,0))):
    print "Test 0/travel_time passed"

# Test 1: If phi and speed symmetric, periodic=True/False should give the same result
phi = -1*np.ones_like(X); phi[X**2+Y**2<0.5] = 1.0
speed = np.ones_like(X); speed[X**2+Y**2<1.0] = 2.0
if np.allclose(skfmm.distance(phi),skfmm.distance(phi,periodic=True)) and np.allclose(skfmm.distance(phi),skfmm.distance(phi,periodic=(1,1))):
    print "Test 1/distance passed"
if np.allclose(skfmm.travel_time(phi,speed),skfmm.travel_time(phi,speed,periodic=True)) and np.allclose(skfmm.travel_time(phi,speed),skfmm.travel_time(phi,speed,periodic=(1,1))):
    print "Test 1/travel_time passed"

# Test 2: With periodic boundary conditions, rolling the numpy array should commute with distance and travel_time functions
phi = -1*np.ones_like(X); phi[X**2+(Y-0.9)**2<0.5] = 1.0
speed = np.ones_like(X); speed[(X-0.9)**2+Y**2<1.0] = 2.0
if np.allclose(skfmm.distance(np.roll(phi,137,axis=0),periodic=True),np.roll(skfmm.distance(phi,periodic=True),137,axis=0)):
    print "Test 2/distance passed"
if np.allclose(skfmm.travel_time(np.roll(phi,-77,axis=1),np.roll(speed,-77,axis=1),periodic=True),np.roll(skfmm.travel_time(phi,speed,periodic=True),-77,axis=1)):
    print "Test 2/travel_time passed"

# Test 3: Visual inspection - also serving as example

X, Y = np.meshgrid(np.linspace(-1,1,200), np.linspace(-1,1,200))
phi = -1*np.ones_like(X)

phi[X>-0.5] = 1
phi[np.logical_and(np.abs(Y)<0.25, X>-0.75)] = 1
speed = np.ones_like(X)
speed[Y>0] = 1.5
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

t    = skfmm.travel_time(phi, speed, dx=1e-2, periodic=(0,1))
pl.title('Travel time from the boundary with an obstacle - periodic in x-direction')
pl.contour(X, Y, phi, [0], linewidths=(3), colors='black')
pl.contour(X, Y, phi.mask, [0], linewidths=(3), colors='red')
pl.contour(X, Y, t, 15)
pl.colorbar()
pl.savefig('2d_phi_travel_time_mask_periodic01.png')
pl.show()

t    = skfmm.travel_time(phi, speed, dx=1e-2, periodic=(1,0))
pl.title('Travel time from the boundary with an obstacle - periodic in y-direction')
pl.contour(X, Y, phi, [0], linewidths=(3), colors='black')
pl.contour(X, Y, phi.mask, [0], linewidths=(3), colors='red')
pl.contour(X, Y, t, 15)
pl.colorbar()
pl.savefig('2d_phi_travel_time_mask_periodic10.png')
pl.show()

t    = skfmm.travel_time(phi, speed, dx=1e-2, periodic=True)
pl.title('Travel time from the boundary with an obstacle - periodic in both directions')
pl.contour(X, Y, phi, [0], linewidths=(3), colors='black')
pl.contour(X, Y, phi.mask, [0], linewidths=(3), colors='red')
pl.contour(X, Y, t, 15)
pl.colorbar()
pl.savefig('2d_phi_travel_time_mask_periodic11.png')
pl.show()
