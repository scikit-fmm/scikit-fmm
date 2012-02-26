import numpy as np
from skfmm import distance, travel_time

# test non-spherical initial shape
# test variable speed


def test():
    """ simple test two point test"""
    np.testing.assert_array_equal(distance([-1, 1]),
                                  [-0.5, 0.5])

    np.testing.assert_allclose(distance([-1, -1, -1, 1, 1, 1]),
                               [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])

    np.testing.assert_allclose(distance([1, 1, 1, -1, -1, -1]),
                               [2.5, 1.5, 0.5, -0.5, -1.5, -2.5])


def test_simple_case():
    """ more simple testing -- three point test """
    np.testing.assert_array_equal(distance([-1, 0, 1]),           [-1, 0, 1])
    np.testing.assert_array_equal(distance([-1, 0, 1], dx=[2]),   [-2, 0, 2])
    np.testing.assert_array_equal(distance([-1, 0, 1], dx=2),     [-2, 0, 2])
    np.testing.assert_array_equal(distance([-1, 0, 1], dx=2.0),   [-2, 0, 2])

    np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1]),
                                  [1, 0, 1])
    np.testing.assert_array_equal(travel_time([-1, 0, 1], [1, 1, 1]),
                                  [1, 0, 1])
    np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1], dx=2),
                                  [2, 0, 2])
    np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1], dx=[2]),
                                  [2, 0, 2])
    np.testing.assert_array_equal(travel_time([1, 0, -1], [1, 1, 1], dx=2.0),
                                  [2, 0, 2])


def test_traveltime_0():
    np.testing.assert_allclose(travel_time([0, 1, 1, 1, 1], [2, 2, 2, 2, 2]),
                               [0, 0.5, 1.0, 1.5, 2.0])
    np.testing.assert_array_equal(travel_time([1, 0, -1], [2, 2, 2]),
                                  [0.5, 0, 0.5])


def test_traveltime_1():
    phi   = [1, 1, 1, -1, -1, -1]
    t     = travel_time(phi, np.ones_like(phi))
    exact = [2.5, 1.5, 0.5, 0.5, 1.5, 2.5]
    np.testing.assert_allclose(t, exact)


def test_traveltime_2():
    phi   = [-1, -1, -1, 1, 1, 1]
    t     = travel_time(phi, np.ones_like(phi))
    exact = [2.5, 1.5, 0.5, 0.5, 1.5, 2.5]
    np.testing.assert_allclose(t, exact)


def test_corner_case0():
    """ test corner case: zero input """
    np.testing.assert_array_equal(distance([0, 0]), [0, 0])
    np.testing.assert_array_equal(travel_time([0, 0], [1, 1]), [0, 0])


def test_zero_dx():
    try:
        distance([1, 0, 1, 1], 0)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "dx must be greater than zero"
    else:
        raise Exception("no exception thrown")


def test_dx_shape():
    try:
        distance([0, 0, 1, 0, 0], [0, 0, 1, 0, 0])
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "dx must be of length len(phi.shape)"
    else:
        raise Exception("no exception thrown")


def test_corner_case1():
    """
    Test catching speeds which are too small. Speeds less than the
    machine epsilon are masked off to avoid an overflow
    """
    t = travel_time([-1, -1, 0, 1, 1], [1, 1, 1, 1, 0])
    assert isinstance(t, np.ma.MaskedArray)
    np.testing.assert_array_equal(t.data[:-1], [2, 1, 0, 1])
    np.testing.assert_array_equal(t.mask, [False, False, False, False, True])

    t2 = travel_time([-1, -1, 0, 1, 1], [1, 1, 1, 1, 1e-300])
    np.testing.assert_array_equal(t, t2)


def test_corner_case2():
    """
    Test that when the mask cuts off the solution, the cut off
    points are also masked.
    """
    ma    = np.ma.MaskedArray([1, 1, 1, 0], [False, True, False, False])
    d     = distance(ma)
    exact = np.ma.MaskedArray([0, 0, 1, 0], [True, True, False, False])
    np.testing.assert_array_equal(d.mask, exact.mask)
    np.testing.assert_array_equal(d, exact)


def test0():
    """circular level set"""
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    d     = distance(phi, dx, self_test=True)
    exact = np.sqrt(X ** 2 + Y ** 2) - r
    np.testing.assert_allclose(d, exact, atol=dx)


def test1():
    """planar level set"""
    N         = 50
    X, Y      = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    dx        = 2.0 / (N - 1)
    phi       = np.ones_like(X)
    phi[0, :] = -1
    d         = distance(phi, dx, self_test=True)
    exact     = Y + 1 - dx / 2.0
    np.testing.assert_allclose(d, exact)


def test2():
    """masked input"""
    N         = 50
    X, Y      = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    dx        = 2.0 / (N - 1)
    phi       = np.ones_like(X)
    phi[0, 0] = -1
    mask      = np.logical_and(abs(X) < 0.25, abs(Y) < 0.25)
    mphi      = np.ma.MaskedArray(phi.copy(), mask)
    d0        = distance(phi, dx, self_test=True)
    d         = distance(mphi, dx, self_test=True)
    d0[mask]  = 0
    d[mask]   = 0
    shadow    = d0 - d
    bsh       = abs(shadow) > 0.001
    diff      = (bsh).sum()

    assert diff > 635 and diff < 645


def test3():
    """test eikonal solution"""
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi) * 2
    t     = travel_time(phi, speed, dx)
    exact = 0.5 * np.abs(np.sqrt(X ** 2 + Y ** 2) - 0.5)

    np.testing.assert_allclose(t, exact, atol=dx)


def test4():
    """test 1d"""
    N          = 100
    X          = np.linspace(-1.0, 1.0, N)
    dx         = 2.0 / (N - 1)
    phi        = np.zeros_like(X)
    phi[X < 0] = -1
    phi[X > 0] = 1
    d          = distance(phi, dx, self_test=True)

    np.testing.assert_allclose(d, X)


def test5():
    """test 3d"""
    N            = 15
    X            = np.linspace(-1, 1, N)
    Y            = np.linspace(-1, 1, N)
    Z            = np.linspace(-1, 1, N)
    phi          = np.ones((N, N, N))
    phi[0, 0, 0] = -1.0
    dx           = 2.0 / (N - 1)
    d            = distance(phi, dx, self_test=True)

    exact        = np.sqrt((X + 1) ** 2 +
                           (Y + 1)[:, np.newaxis] ** 2 +
                           (Z + 1)[:, np.newaxis, np.newaxis] ** 2)

    np.testing.assert_allclose(d, exact, atol=dx)


def test6():
    """ test default dx """
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi) * 2
    travel_time(phi, speed, self_test=True)


def test_dx():
    """ test non-square grid and dx different in different directions """
    N      = 50
    NX, NY = N, 5 * N
    X, Y   = np.meshgrid(np.linspace(-1, 1, NY), np.linspace(-1, 1, NX))
    r      = 0.5
    phi    = X ** 2 + Y ** 2 - r ** 2
    dx     = [2.0 / (NX - 1), 2.0 / (NY - 1)]
    d      = distance(phi, dx, self_test=True)
    exact  = np.sqrt(X ** 2 + Y ** 2) - r

    np.testing.assert_allclose(d, exact, atol=max(dx))


def test7():
    """ no zero level set test """
    try:
        distance([1, 1], self_test=False)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == \
            "the array phi contains no zero contour (no zero level set)"
    else:
        raise Exception("no exception thrown")


def test8():
    """ shape mismatch test """
    try:
        travel_time([-1, 1], [2])
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "phi and speed must have the same shape"
    else:
        raise Exception("no exception raised")


def test8_2():
    """ speed wrong type test """
    try:
        travel_time([0, 0, 1, 1], 2)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "speed must be a 1D to 12-D array of doubles"
    else:
        raise Exception("no exception raised")


def test9():
    """ dx mismatch test """
    try:
        travel_time([-1, 1], [2, 2], [2, 2, 2, 2])
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "dx must be of length len(phi.shape)"
    else:
        raise Exception("no exception raised")


def test10():
    """ test c error handling """
    try:
        distance([-1, 1], self_test=44)
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
    else:
        raise Exception("no exception raised")


def test11():
    """ check array type test """
    try:
        distance(np.array(["a", "b"]))
    except Exception as ex:
        assert ValueError.__name__ == ex.__class__.__name__
        assert ex.args[0] == "phi must be a 1 to 12-D array of doubles"
    else:
        raise Exception("no exception raised")


if __name__ == '__main__':
    import nose
    nose.main()
