.. scikit-fmm documentation master file, created by
   sphinx-quickstart on Wed Feb  8 06:45:28 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

scikit-fmm documentation
========================

:py:mod:`scikit-fmm` is an extension module which implements the fast
marching method.

The fast marching method is used to model the evolution of boundaries
and interfaces in a variety of application areas.

More specifically, the fast marching method is a numerical technique for
finding approximate solutions to boundary value problems of the
Eikonal equation,

.. math::
   F(x) | \nabla T(x) | = 1

where :math:`F>0` and the unknown :math:`T(x)` is interpreted as the
travel time from an initial boundary to a point x. The initial
location of the boundary is defined by the zero contour (or zero
level-set) of a scalar function.

In this document the scalar function containing the initial interface
location is referred to as phi. The scalar function phi can be thought to
exist in a dimension higher than the boundary of interest and only the
zero contour of the function is physically meaningful. The boundary
grows outward in the local normal direction at a speed given by
:math:`F(x)`.

:py:mod:`scikit-fmm` is a simple module which provides two functions:
:py:func:`skfmm.distance` and :py:func:`skfmm.travel_time`. The import
name of :py:mod:`scikit-fmm` is :py:mod:`skfmm`.


Examples
========
First, a simple one-dimensional example::

    >>> import skfmm
    >>> phi = [-1,-1,-1,1,1,1]
    >>> skfmm.distance(phi)
    array([-2.5, -1.5, -0.5,  0.5,  1.5,  2.5])

Here the zero contour of phi is between elements 2 and 3. The return
value of :py:func:`skfmm.distance` gives the signed distance from zero
contour. No grid spacing is given, so it is taken as 1. To specify a
spacing use the dx argument::

    >>> skfmm.distance(phi, dx=0.25)
    array([-0.625, -0.375, -0.125,  0.125,  0.375,  0.625])

A two-dimensional example:

.. image:: 2d_phi.png

The boundary is specified as the zero contour of a scalar function phi:
::

 >>> import numpy as np
 >>> import pylab as pl
 >>> X, Y = np.meshgrid(np.linspace(-1,1,200), np.linspace(-1,1,200))
 >>> phi = -1*np.ones_like(X)
 >>> phi[X>-0.5] = 1
 >>> phi[np.logical_and(np.abs(Y)<0.25, X>-0.75)] = 1

.. image:: 2d_phi_distance.png

::

 >>> d = skfmm.distance(phi, dx=1e-2)

:py:mod:`scikit-fmm` can also calculate travel times from an interface
given an array containing the interface propogation speed at each
point. Using the same initial interface position as above we set the
interface propagation speed to be 1.5 times greater in the upper half
of the domain.

.. image:: 2d_phi_travel_time.png

::

 >>> speed = np.ones_like(X)
 >>> speed[Y>0] = 1.5
 >>> t = skfmm.travel_time(phi, speed, dx=1e-2)

both :py:func:`skfmm.travel_time` and :py:func:`skfmm.distance`
support masked arrays for input. This allows an obstacle to be introduced.

.. image:: 2d_phi_travel_time_mask.png

::

 >>> mask = np.logical_and(abs(X)<0.1, abs(Y)<0.5)
 >>> phi  = np.ma.MaskedArray(phi, mask)
 >>> t    = skfmm.travel_time(phi, speed, dx=1e-2)

The full example is in examples/2d_example.py.
:doc:`examples`



Limitations:
============
:py:mod:`scikit-fmm` only works for regular Cartesian grids, but grid cells may
have a different (uniform) length in each dimension.

:py:mod:`scikit-fmm` does not compute level set extension or solve the
time-dependent level set equation.


Function Reference
==================

.. autofunction:: skfmm.distance

.. autofunction:: skfmm.travel_time


.. toctree::
   :maxdepth: 2

   examples


