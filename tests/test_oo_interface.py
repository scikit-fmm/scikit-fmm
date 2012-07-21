import numpy as np
from skfmm import distance_variable, distance

# test non-spherical initial shape
# test variable speed


def test():
    """ OO interface test"""
    phi = np.array([-1,1], dtype="double")
    dv = distance_variable(phi, [1])
    dv.calculate()
    np.testing.assert_array_equal(dv.distance,
                                  [-0.5, 0.5])
def test_o2_distance_dw__():
    """ OO more detailed test """
    phi = np.array(((-1, -1, 1, 1),
                   (-1, -1, 1, 1),
                   (1, 1, 1, 1),
                   (1, 1, 1, 1)), dtype='double')
    dv = distance_variable(phi, [1.0, 1.0])
    dv.calculate()
    dw_o2 = [[-1.30473785,  -0.5,          0.5,         1.49923009],
             [-0.5,        -0.35355339,    0.5,         1.45118446],
             [ 0.5,         0.5,         0.97140452,  1.76215286],
             [ 1.49923009,  1.45118446,  1.76215286,  2.33721352]]
    np.testing.assert_allclose(dv.distance, dw_o2)

phi = np.ones((400,400), dtype=np.float64)
phi[0,0] = -1.0
dv = distance_variable(phi, [1.0, 1.0])


