![Linux build](https://github.com/scikit-fmm/scikit-fmm/actions/workflows/build.yaml/badge.svg) ![Windows build](https://github.com/scikit-fmm/scikit-fmm/actions/workflows/windows_build.yaml/badge.svg) [![PyPI version](https://badge.fury.io/py/scikit-fmm.svg)](http://pypi.python.org/pypi/scikit-fmm)[![Documentation Status](https://readthedocs.org/projects/scikit-fmm/badge/?version=latest)](https://scikit-fmm.readthedocs.io/en/latest/?badge=latest) ![Contributors](https://img.shields.io/github/contributors/scikit-fmm/scikit-fmm.svg)<a href="https://pepy.tech/project/scikit-fmm"><img alt="Downloads" src="https://pepy.tech/badge/scikit-fmm"></a>

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

### Requirements
* NumPy
* Building requires the pypa/build module
  (https://github.com/pypa/build) and a C/C++ compiler
### Bugs, questions, patches, feature requests, discussion & cetera
* Open a GitHub pull request or a GitHub issue
* Email list: http://groups.google.com/group/scikit-fmm
  * Send an email to scikit-fmm+subscribe@googlegroups.com to subscribe.

### Installing
* Via pip: `pip install scikit-fmm`
* Anaconda linux-64 and linux-ppc64le packages:
  * `conda install scikit-fmm`
* Ubuntu PPA
  * https://launchpad.net/~nvidia-digits/+archive/ubuntu/dev
* Debian
  * https://tracker.debian.org/pkg/scikit-fmm

### Building and installing from Source
* `pip install build`
* `python -m build`
* `pip install .`


### Running Tests
* `python -c "import skfmm; skfmm.test(True)"` (Do not run the tests from the source directory.)
* Tests are doctests in `skfmm/__init__.py`

### Building documentation
* Requires sphinx and numpydoc
* `make html`

### Publications using scikit-fmm

* Akinola, I., J Varley, B. Chen, and P.K. Allen
  (2018) "Workspace Aware Online Grasp Planning" arXiv:1806.11402v1
  [cs.RO] 29 Jun 2018 https://arxiv.org/pdf/1806.11402.pdf

* Bortolussi, V., B. Figliuzzi, F. Willot, M.
  Faessel, M. Jeandin (2018) "Morphological modeling of cold spray
  coatings" Image Anal Stereol 2018;37:145-158 doi:10.5566/ias.1894
  https://hal.archives-ouvertes.fr/hal-01837906/document

* Chalmers, S., C.D. Saunter, J.M. Girkin and J.G. McCarron (2016)
  "Age decreases mitochondrial motility and increases mitochondrial
  size in vascular smooth muscle." Journal of Physiology, 594.15 pp
  4283–4295.

* Diogo Brandão Amorim (2014) "Efficient path planning of a mobile
  robot on rough terrain" Master's Thesis, Department of Aerospace
  Engineering, University of Lisbon.

* Giometto, A., D.R. Nelson, and A.W. Murray (2018)
  "Physical interactions reduce the power of natural selection in
  growing yeast colonies", PNAS November 6, 2018 115 (45) 11448-11453;
  published ahead of print October 23, 2018
  https://doi.org/10.1073/pnas.1809587115

* Joshua A. Taillon, Christopher Pellegrinelli, Yilin Huang, Eric D.
  Wachsman, and Lourdes G. Salamanca-Riba (2014) "Three Dimensional
  Microstructural Characterization of Cathode Degradation in SOFCs
  Using Focused Ion Beam and SEM" ECS Trans. 2014 61(1): 109-120;
  https://www.joshuataillon.com/pdfs/2015-08-06%20jtaillon%203D%20SOFC%20cathode%20degradation.pdf

* Marshak, C., I. Yanovsky, and L. Vese (2017) "Energy
  Minimization for Cirrus and Cumulus Cloud Separation in Atmospheric
  Images" IGARSS 2018 - 2018 IEEE International Geoscience and Remote
  Sensing Symposium DOI: 10.1109/IGARSS.2018.8517940
  ftp://ftp.math.ucla.edu/pub/camreport/cam17-68.pdf

* Moon, K. R., V. Delouille, J.J. Li, R. De Visscher, F. Watson and
  A.O. Hero III (2016) "Image patch analysis of sunspots and active
  regions." J. Space Weather Space Clim., 6, A3, DOI:
  10.1051/swsc/2015043.

* Tao, M., J. Solomon and A. Butscher (2016) "Near-Isometric Level Set
  Tracking." in Eurographics Symposium on Geometry Processing 2016
  Eds: M. Ovsjanikov and D. Panozzo. Volume 35 (2016), Number 5

* Thibaut, R., Laloy, E., Hermans, T., 2021. A new framework for
  experimental design using Bayesian Evidential Learning: The
  case of wellhead protection area. J. Hydrol. 603, 126903.
  https://doi.org/10.1016/j.jhydrol.2021.126903

* Vargiu, Antioco, M. Marrocu, L. Massidda (2015) "Implementazione e
  valutazione su un caso reale del servizio di Cloud Computing per la
  simulazione di incendi boschivi in Sardegna" (Implementation and
  evaluation on a real case of Cloud computing service for simulation
  of Forest fires in Sardinia). Sardinia Department of Energy and
  Environment. CRS4 PIA 2010 D5.4.

* Wronkiewicz, M. (2018) "Mapping buildings with help from machine
  learning" Medium article, June 29th 2018
  https://medium.com/devseed/mapping-buildings-with-help-from-machine-learning-f8d8d221214a

* Makki, K., Ben Salem, D., Ben Amor, B. (2021) "Toward the Assessment
  of Intrinsic Geometry of Implicit Brain MRI Manifolds" IEEE Access,
  volume 9, pages  131054 - 131071 (September 2021)
  DOI: 10.1109/ACCESS.2021.3113611
  https://ieeexplore.ieee.org/abstract/document/9540688

### Version History:
* 0.0.1: February 13 2012
  * Initial release
* 0.0.2: February 26th 2012
  * Including tests and docs in source distribution. Minor changes to
    documentation.
* 0.0.3: August 4th 2012
  * Extension velocities.
  * Fixes for 64 bit platforms.
  * Optional keyword argument for point update order.
  * Bug reports and patches from three contributors.
* 0.0.4: October 15th 2012
   * Contributions from Daniel Wheeler:
     * Bug fixes in extension velocity.
     * Many additional tests and migration to doctest format.
     * Additional optional input to extension_velocities() for FiPy compatibly.
* 0.0.5: May 12th 2014
   * Fix for building with MSVC (Jan Margeta).
   * Corrected second-order point update.
* 0.0.6: February 20th 2015
   * Documentation clarification (Geordie McBain).
   * Python 3 port (Eugene Prilepin).
   * Python wrapper for binary min-heap.
   * Freeze equidistant narrow-band points simultaneously.
* 0.0.7: October 21st 2015
   * Bug fix to upwind finite difference approximation for negative
     phi from Lester Hedges.
* 0.0.8: March 9th 2016
   * Narrow band capability: an optional "narrow" keyword argument
     limits the extent of the marching algorithm (Adrian Butscher).
* 0.0.9: August 5th 2016
   * Periodic boundaries: an optional "periodic" keyword argument
     enables periodic boundaries in one or more directions (Wolfram Moebius).
* 2019.1.30 January 30th 2019
   * Abrupt change to version numbering scheme.
   * Bugfix in setup.py to allow installing via pip with numpy (ManifoldFR).
   * Handle C++ exceptions during fast marching (Jens Glaser).
   * Accept a zero discriminant in travel time point update.
* 2021.1.20 January 20th 2021
   * Fix divide by zero bugs in travel_time and extension_velocities
   * Contributions from Murray Cutforth, f-fanni, and okonrad
* 2021.1.21 January 21st 2021
   * Minor C++ change (removed the auto keyword) to fix the compile on TravisCI.
* 2021.2.2 February 2nd 2021
   * Add a pyproject.toml file to specify numpy as a build
     requirement, this is needed to build with new version of pip
     (David Parsson).
* 2021.7.8 July 8th 2021
   * Add a pyproject.toml file to the MANIFEST.in file to fix the
     numpy build dependency (David Parsson). Fix numpy deprecation
     warnings and improve source code formatting (Robin Thibaut).
* 2021.9.23 September 23rd 2021
   * Make the pyproject.toml file specify the oldest supported
     numpy as a build requirement, to allow using wheels with any
     numpy version.
     (David Parsson).
* 2021.10.29 October 29th 2021
   * Fix for point update discriminant exactly equal to zero
   * Fall back calculation for point update when discriminant becomes negative
   * (Joshua Gehre)
* 2022.02.02 February 2nd 2022
   * Fixes for Python 3.10 compatibility
   * (Amin Sadeghi, Xylar Asay-Davis, David Parsson)
* 2022.03.26 March 26th 2022
   * Following the breaking changes in setuptools v61.0.0 it is
     suggested to set py_modules to disable auto-discovery behavior.
   * (Daniel Ammar)
* 2022.08.15 August 15th 2022
   * Following the breaking changes in setuptools v65 pin setuptools
     to v64
   * (DorSSS)
* 2022.04.02 April 2nd 2023
   * Build fixes for Python 3.11 (update pheap cython wrapper)
   * No solver changes
* 2024.05.29 May 29th 2024
   * Update build system to use meson
   * Python 3.12 support
   * No solver changes
* 2025.01.29 January 29th 2025
   * NumPy 2.0 support on Linux and Windows
   * Update to Github workflows to use v4 api
   * Windows wheels
   * Support for Python 3.13
* 2025.01.29 June 23rd 2025
   * Fixes for deprecated NumPy API and Python 3.11
   * Github Workflow fixes



Copyright 2025 The scikit-fmm team.

BSD-style license. See LICENSE.txt in the source directory.
