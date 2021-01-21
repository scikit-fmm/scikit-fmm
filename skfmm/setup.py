from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    config = Configuration('skfmm', parent_package, top_path)
    config.add_extension("cfmm",
                          sources=["fmm.cpp",
                                   "heap.cpp",
                                   "base_marcher.cpp",
                                   "distance_marcher.cpp",
                                   "travel_time_marcher.cpp",
                                   "extension_velocity_marcher.cpp"],
                         include_dirs=['.'])
    config.add_extension("pheap", sources=["pheap.cpp", "heap.cpp"])
    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
