# Version history of scikit-fmm

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
* 2025.6.23 June 23rd 2025
   * Fixes for deprecated NumPy API and Python 3.11
   * Github Workflow fixes
