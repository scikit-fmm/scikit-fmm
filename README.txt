scikit-fmm is an extension module which implements the fast marching method.

The fast marching method is used to model the evolution of boundaries
and interfaces in a variety of application areas.

scikit-fmm is a simple module which provides two functions:
distance(phi, dx=1.0) and travel_time(phi, speed, dx=1.0).

The functions calculate the signed distance and travel time to an
interface described by the zero contour of the input array phi.

>>> import skfmm
>>> phi = [-1,-1,-1,1,1,1]
>>> skfmm.distance(phi)
array([-2.5, -1.5, -0.5,  0.5,  1.5,  2.5])

>>> skfmm.travel_time(phi, [2,2,2,2,2,2])
array([ 1.25,  0.75,  0.25,  0.25,  0.75,  1.25])

The input array can be of 1, 2, 3 or higher dimensions and can be a
masked array.

Documentation: http://packages.python.org/scikit-fmm

PyPI: http://pypi.python.org/pypi/scikit-fmm

Source Code: https://github.com/scikit-fmm/

Requirements: Numpy and a C/C++ compiler (gcc/MinGW)

Installing:
 $ python setup.py install

Testing (requires nose):
  $ python tests/test_fmm.py

Building documentation (require sphinx):
  $ make html
