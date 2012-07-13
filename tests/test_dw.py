import numpy as np
from skfmm import distance, travel_time

def test_o1_distance_dw():
    phi = np.array(((-1, -1, 1, 1),
                   (-1, -1, 1, 1),
                   (1, 1, 1, 1),
                   (1, 1, 1, 1)))
    o1 = distance(phi, order=1, self_test=True)
    dw_o1 =   [[-1.20710678, -0.5,         0.5,         1.5],
               [-0.5,        -0.35355339,  0.5,         1.5],
               [ 0.5,         0.5,         1.20710678,  2.04532893],
               [ 1.5,         1.5,         2.04532893,  2.75243571]]
    np.testing.assert_allclose(o1, dw_o1)


def test_o1_traveltime_dw():
    phi = np.array(((-1, -1, 1, 1),
                   (-1, -1, 1, 1),
                   (1, 1, 1, 1),
                   (1, 1, 1, 1)))
    o1 = travel_time(phi, np.ones_like(phi), order=1, self_test=True)
    dw_o1 =   [[-1.20710678, -0.5,         0.5,         1.5],
               [-0.5,        -0.35355339,  0.5,         1.5],
               [ 0.5,         0.5,         1.20710678,  2.04532893],
               [ 1.5,         1.5,         2.04532893,  2.75243571]]
    np.testing.assert_allclose(o1, np.abs(dw_o1))


def test_o2_distance_dw():
    phi = np.array(((-1, -1, 1, 1),
                   (-1, -1, 1, 1),
                   (1, 1, 1, 1),
                   (1, 1, 1, 1)))
    o2 = distance(phi, self_test=True)
    dw_o2 = [[-1.30473785,  -0.5,          0.5,         1.49923009],
             [-0.5,        -0.35355339,    0.5,         1.45118446],
             [ 0.5,         0.5,         0.97140452,  1.76215286],
             [ 1.49923009,  1.45118446,  1.76215286,  2.33721352]]
    np.testing.assert_allclose(o2, dw_o2)


def test_o2_travel_time_dw():
    phi = np.array(((-1, -1, 1, 1),
                   (-1, -1, 1, 1),
                   (1, 1, 1, 1),
                   (1, 1, 1, 1)))
    o2 = travel_time(phi, np.ones_like(phi), self_test=True)
    dw_o2 = [[-1.30473785,  -0.5,          0.5,         1.49923009],
             [-0.5,        -0.35355339,    0.5,         1.45118446],
             [ 0.5,         0.5,         0.97140452,  1.76215286],
             [ 1.49923009,  1.45118446,  1.76215286,  2.33721352]]
    np.testing.assert_allclose(o2, np.abs(dw_o2))


def test_order_error():
    try:
        distance([-1,1], order=0)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "order must be 1 or 2"
    else:
        raise Exception("no exception thrown")

    try:
        distance([-1,1], order=3)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "order must be 1 or 2"
    else:
        raise Exception("no exception thrown")


