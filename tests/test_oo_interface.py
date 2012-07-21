import numpy as np
from skfmm import distance_variable

# test non-spherical initial shape
# test variable speed


def test():
    """ simple test two point test"""
    phi = np.array([-1,1])
    dv = distance_variable(phi)
    dv.calculate()
    np.testing.assert_array_equal(dv.distance,
                                  [-0.5, 0.5])
