.. py:module:: skfmm

scikit-fmm documentation
========================

:py:obj:`scikit-fmm` is a python extension module which implements the
fast marching method.

The fast marching method is used to model the evolution of boundaries
and interfaces in a variety of application areas.

More specifically, the fast marching method is a numerical technique for
finding approximate solutions to boundary value problems of the
Eikonal equation,

.. math::
   F(x) | \nabla T(x) | = 1

Typically, such a problem describes the evolution of a closed curve as
a function of time :math:`T` with speed :math:`F(x)>0` in the normal
direction at a point x on the curve. The speed function is specified,
and the time at which the contour crosses a point x is obtained by
solving the equation. The initial location of the boundary is defined
by the zero contour (or zero level-set) of a scalar function.

In this document the scalar function containing the initial interface
location is referred to as phi. The scalar function phi can be thought to
exist in a dimension higher than the boundary of interest and only the
zero contour of the function is physically meaningful. The boundary
grows outward in the local normal direction at a speed given by
:math:`F(x)`.

:py:obj:`scikit-fmm` is a simple module which provides the functions:
:py:func:`distance`, :py:func:`travel_time` and
:py:func:`extension_velocities`. The import name of
:py:obj:`scikit-fmm` is :py:mod:`skfmm`.


Examples
========
First, a simple example::

    >>> import skfmm
    >>> import numpy as np
    >>> phi = np.ones((3, 3))
    >>> phi[1, 1] = -1
    >>> skfmm.distance(phi)
    array([[ 1.20710678,  0.5       ,  1.20710678],
           [ 0.5       , -0.35355339,  0.5       ],
           [ 1.20710678,  0.5       ,  1.20710678]])

Here the zero contour of phi is around the (1, 1) point. The return
value of :py:func:`distance` gives the signed distance from zero
contour. No grid spacing is given, so it is taken as 1. To specify a
spacing use the optional dx argument::

    >>> skfmm.distance(phi, dx=0.25)
    array([[ 0.3017767 ,  0.125     ,  0.3017767 ],
           [ 0.125     , -0.08838835,  0.125     ],
           [ 0.3017767 ,  0.125     ,  0.3017767 ]])

A more detailed example:

.. image:: 2d_phi.png

The boundary is specified as the zero contour of a scalar function phi:

::

 >>> import numpy as np
 >>> import pylab as pl
 >>> X, Y = np.meshgrid(np.linspace(-1,1,200), np.linspace(-1,1,200))
 >>> phi = -1 * np.ones_like(X)
 >>> phi[X > -0.5] = 1
 >>> phi[np.logical_and(np.abs(Y) < 0.25, X > -0.75)] = 1
 >>> d = skfmm.distance(phi, dx=1e-2)

.. image:: 2d_phi_distance.png

:py:obj:`scikit-fmm` can also calculate travel times from an interface
given an array containing the interface propogation speed at each
point. Using the same initial interface position as above we set the
interface propagation speed to be 1.5 times greater in the upper half
of the domain.

::

 >>> speed = np.ones_like(X)
 >>> speed[Y > 0] = 1.5
 >>> t = skfmm.travel_time(phi, speed, dx=1e-2)

.. image:: 2d_phi_travel_time.png

Consider an obstacle within which the speed is vanishing. In principle this may
lead to singularities (division by zero) in the update algorithm. Therefore,
both :py:func:`travel_time` and :py:func:`distance` support masked arrays for
input. This allows an obstacle to be introduced. (Note that if the speed in a cell
is less than machine precision, a cell is masked internally to prevent division by 0.)

::

 >>> mask = np.logical_and(abs(X) < 0.1, abs(Y) < 0.5)
 >>> phi  = np.ma.MaskedArray(phi, mask)
 >>> t    = skfmm.travel_time(phi, speed, dx=1e-2)

.. image:: 2d_phi_travel_time_mask.png

The distance function, travel time or extension velocities can be
limited to with in a narrow band around the zero contour by specifying
the `narrow` keyword.

::

 >>> phi = -1 * np.ones_like(X)
 >>> phi[X > -0.5] = 1
 >>> phi[np.logical_and(np.abs(Y) < 0.25, X > -0.75)] = 1
 >>> d = skfmm.distance(phi, dx=1e-2, narrow=0.3)

.. image:: 2d_phi_distance_narrow.png

The full example is in examples/2d_example.py.
:doc:`examples`



An example of using periodic boundary conditions.

::

 >>> X, Y = np.meshgrid(np.linspace(-1,1,501), np.linspace(-1,1,501))
 >>> phi = (X+0.8)**2+(Y+0.8)**2 - 0.01
 >>> speed = 1+X**2+Y**2
 >>> skfmm.distance(phi, dx=2.0/500)
 >>> skfmm.distance(phi, dx=2.0/500, periodic=True)
 >>> skfmm.travel_time(phi, speed, dx=2.0/500, periodic=(1,0))

.. image:: periodic.png

The full example is in examples/boundaryconditions_example.py
:doc:`examples`

An example of using :py:obj:`scikit-fmm` to compute extension velocities.

::

 >>> N     = 150
 >>> X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
 >>> r     = 1.75
 >>> dx    = 2.0 / (N - 1)
 >>> phi   = (X) ** 2 + (Y+1.85) ** 2 - r ** 2
 >>> speed = X + 1.25
 >>> d, f_ext = extension_velocities(phi, speed, dx)

.. image:: extension_velocity.png

The full example is in examples/extension_velocities_example.py.
:doc:`examples`

.. toctree::
   :maxdepth: 2

   examples
   testing

Limitations:
============
:py:obj:`scikit-fmm` only works for regular Cartesian grids, but grid cells may
have a different (uniform) length in each dimension.

Function Reference
==================

.. autofunction:: distance

.. autofunction:: travel_time

.. autofunction:: extension_velocities

.. autoclass:: heap
   :members:

Testing
=======

To run all the tests use

   $ python -c "import skfmm; skfmm.test()"

See the full :doc:`testing`.
