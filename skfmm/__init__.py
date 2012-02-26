"""
scikit-fmm is an extension module which implements the fast marching method.

The fast marching method is used to model the evolution of boundaries
and interfaces in a variety of application areas. More specifically,
the fast marching method is a numerical technique for finding
approximate solutions to boundary value problems of the Eikonal
equation:

F(x) | grad T(x) | = 1.

Typically, such a problem describes the evolution of a closed curve as
a function of time T with speed F(x)>0 in the normal direction at a
point x on the curve. The speed function is specified, and the time at
which the contour crosses a point x is obtained by solving the
equation.

scikit-fmm is a simple module which provides two functions:
distance(phi, dx=1.0) and travel_time(phi, speed, dx=1.0).

The functions calculate the signed distance and travel time to an
interface described by the zero contour of the input array phi.

>>> import skfmm
>>> import numpy as np
>>> phi = np.ones((3, 3))
>>> phi[1, 1] = -1
>>> skfmm.distance(phi)
array([[ 1.20710678,  0.5       ,  1.20710678],
       [ 0.5       , -0.35355339,  0.5       ],
       [ 1.20710678,  0.5       ,  1.20710678]])

>>> skfmm.travel_time(phi, speed = 3.0 * np.ones_like(phi))
array([[ 0.40236893,  0.16666667,  0.40236893],
       [ 0.16666667,  0.11785113,  0.16666667],
       [ 0.40236893,  0.16666667,  0.40236893]])

The input array can be of 1, 2, 3 or higher dimensions and can be a
masked array.

Documentation: http://packages.python.org/scikit-fmm

PyPI: http://pypi.python.org/pypi/scikit-fmm

Source Code: https://github.com/scikit-fmm/scikit-fmm

Requirements: Numpy and a C/C++ compiler (gcc/MinGW)

Bugs, questions, patches, feature requests:
  Please use the scikit-fmm Github issue tracker.
  https://github.com/scikit-fmm/scikit-fmm/issues

Installing:
 $ python setup.py install

Testing (requires nose):
  $ python tests/test_fmm.py

Building documentation (required sphinx and numpydoc):
  $ make html

Version History:

0.0.1: February 13 2012
  Initial release

0.0.2: February 26th 2012
  Including tests and docs in source distribution. Minor changes to
  documentation.

:Copyright: Copyright 2012 The scikit-fmm team.
:License: BSD-style license. See LICENSE.txt in the scipy source directory.
"""
__version__ = "0.0.2"

from pfmm import distance, travel_time
