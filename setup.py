#!/usr/bin/env python
import os
import sys
import setuptools

long_description = """scikit-fmm is a Python extension module which implements the fast
marching method.

- Signed distance functions
- Travel time transforms (solutions to the Eikonal equation)
- Extension velocities

https://github.com/scikit-fmm/scikit-fmm
"""

DISTNAME         = "scikit-fmm"
DESCRIPTION      = "An extension module implementing the fast marching method"
MAINTAINER       = "Jason Furtney"
MAINTAINER_EMAIL = "jkfurtney@gmail.com"
VERSION          = "2021.1.21"
URL              = 'https://github.com/scikit-fmm/scikit-fmm'
LICENSE          = 'BSD'
KEYWORDS         = "fast marching method, Eikonal equation, interface, boundary"

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True)

    config.add_subpackage('skfmm')

    return config


def parse_setuppy_commands():
    """Check the commands and respond appropriately

    Inspired heavily from SciPy's setup script : https://github.com/scipy/scipy/blob/master/setup.py
    """
    if len(sys.argv) < 2:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands= ['--help-commands',
                        'egg_info',
                        '--version',
                        'clean',
                        'install_egg_info',
                        'rotate'
                   ]

    for command in info_commands:
        if command in sys.argv[1:]:
            return False

    good_commands = ('develop', 'sdist', 'build', 'build_ext', 'build_py',
                     'build_clib', 'build_scripts', 'bdist_wheel', 'bdist_rpm',
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg',
                     'build_sphinx')

    for command in good_commands:
        if command in sys.argv[1:]:
            return True

    # The following commands are supported, but we need to show more
    # useful messages to the user
    if 'install' in sys.argv[1:]:
        print("""
            Note: if you need to uninstall you should `pip install scikit-fmm` instead of using `setup.py install`
            """)
        return True


def setup_package():

    metadata = dict(
        name             = DISTNAME,
        version          = VERSION,
        maintainer       = MAINTAINER,
        maintainer_email = MAINTAINER_EMAIL,
        description      = DESCRIPTION,
        url              = URL,
        license          = LICENSE,
        keywords         = KEYWORDS,
        long_description = long_description,
        configuration    = configuration,
        install_requires = ['numpy >= 1.0.2'],
        classifiers      = ["Development Status :: 5 - Production/Stable",
                            "License :: OSI Approved :: BSD License",
                            "Operating System :: OS Independent",
                            "Topic :: Scientific/Engineering",
                            "Intended Audience :: Science/Research",
                            "Programming Language :: C++",
                            "Programming Language :: Python :: 2",
                            "Programming Language :: Python :: 3"]
    )

    if "--force" in sys.argv:
        run_build = True
        sys.argv.remove('--force')
    else:
        # Raise errors for unsupported commands, improve help output, etc.
        run_build = parse_setuppy_commands()

    from setuptools import setup

    if run_build:
        from numpy.distutils.core import setup

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
