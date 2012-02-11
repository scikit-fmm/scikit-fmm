from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('skfmm', parent_package, top_path)
    config.add_extension("cfmm",
                          sources=["fmm.cpp","heap.cpp", "fast_marching.cpp"],
                          include_dirs=['.'])
    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
