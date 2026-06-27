![Linux build](https://github.com/scikit-fmm/scikit-fmm/actions/workflows/build.yaml/badge.svg) ![Windows build](https://github.com/scikit-fmm/scikit-fmm/actions/workflows/windows_build.yaml/badge.svg) [![PyPI version](https://badge.fury.io/py/scikit-fmm.svg)](http://pypi.python.org/pypi/scikit-fmm)[![Documentation Status](https://readthedocs.org/projects/scikit-fmm/badge/?version=latest)](https://scikit-fmm.readthedocs.io/en/latest/?badge=latest) ![Contributors](https://img.shields.io/github/contributors/scikit-fmm/scikit-fmm.svg)<a href="https://pepy.tech/project/scikit-fmm"><img alt="Downloads" src="https://pepy.tech/badge/scikit-fmm"></a> [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE.txt)

# *Warning* this is an experimental branch of scikit-fmm

This branch uses bicubic interpolation to initialize the initially
frozen cells (adjacent to the zero-levelset). This is a work in progress.

# scikit-fmm: the fast marching method for Python


`scikit-fmm` is a Python extension module which implements the fast marching method.
The fast marching method is used to model the evolution of boundaries
and interfaces in a variety of application areas. More specifically,
the fast marching method is a numerical technique for finding
approximate solutions to boundary value problems of the Eikonal
equation:

F(x) | grad T(x) | = 1

Typically, such a problem describes the evolution of a closed curve as
a function of time T with speed F(x)>0 in the normal direction at a
point x on the curve. The speed function is specified, and the time at
which the contour crosses a point x is obtained by solving the
equation.

### Installing
* Via pip: `pip install scikit-fmm`
  * scikit-fmm is a compiled c extension module. Windows binary wheels
  are provided, Linux pip  install requires the
  compilers to be installed.
* Debian: https://tracker.debian.org/pkg/scikit-fmm
* Ubuntu 25.10 and later: `sudo apt update && sudo apt install python3-scikit-fmm`
* Anaconda packages: https://anaconda.org/channels/conda-forge/packages/scikit-fmm/overview

### Quick Start
scikit-fmm is a simple module which provides functions to calculate
the signed distance and travel time to an interface described by the
zero contour of the input array phi.

```python
import skfmm
import numpy as np
phi = np.ones((3, 3))
phi[1, 1] = -1
skfmm.distance(phi)
```

  ```python
   array([[ 1.20710678,  0.5       ,  1.20710678],
          [ 0.5       , -0.35355339,  0.5       ],
          [ 1.20710678,  0.5       ,  1.20710678]])
  ```
---
```python
skfmm.travel_time(phi, speed = 3.0 * np.ones_like(phi))
```

   ```python
   array([[ 0.40236893,  0.16666667,  0.40236893],
          [ 0.16666667,  0.11785113,  0.16666667],
          [ 0.40236893,  0.16666667,  0.40236893]])
   ```
---

The input array can be of 1, 2, 3 or higher dimensions and can be a
masked array. A function is provided to compute extension velocities.

### Documentation
* http://scikit-fmm.readthedocs.org/en/master/

### PyPI
* http://pypi.python.org/pypi/scikit-fmm


### Bugs, questions, patches, feature requests, discussion & cetera
* GitHub issues: https://github.com/scikit-fmm/scikit-fmm/issues

### Building and installing from Source
Building requires the pypa/build module
(https://github.com/pypa/build) and a C/C++ compiler.
```
pip install build
python -m build
pip install .
```
For more on development and making releases see: [DEV_NOTES.md](DEV_NOTES.md)

### Running Tests
* `python -c "import skfmm; skfmm.test(True)"` (Do not run the tests from the source directory.)
* Tests are doctests in `skfmm/__init__.py`

### Building documentation
* Requires sphinx and numpydoc
* `make html`

### Publications using scikit-fmm
see [PUBLICATIONS.md](PUBLICATIONS.md)

### Version History:
see [CHANGELOG.md](CHANGELOG.md)


Copyright 2025 The scikit-fmm team.

BSD-style license. See LICENSE.txt in the source directory.
