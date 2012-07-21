import numpy as np
from skfmm import distance_variable

# test non-spherical initial shape
# test variable speed


def test():
    """ OO interface test"""
    phi = np.array([-1,1], dtype="double")
    dv = distance_variable(phi, [1])
    dv.calculate()
    np.testing.assert_array_equal(dv.distance,
                                  [-0.5, 0.5])
