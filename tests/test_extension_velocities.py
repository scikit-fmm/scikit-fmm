import numpy as np
from skfmm import extension_velocities, distance


def test_1d_extension_constant():
    phi =   [-1,-1,-1,1,1,1]
    speed = [1,1,1,1,1,1]
    d, f_ext  = extension_velocities(phi, speed, self_test=True)
    np.testing.assert_allclose(speed, f_ext)

def test_1d_extension_block():
    phi =   np.ones(10)
    phi[0] =- 1
    speed = np.ones(10)
    speed[0:3] = 5
    d, f_ext  = extension_velocities(phi, speed, self_test=True)
    np.testing.assert_allclose(f_ext, 5)

def test_exntension0():
    """
    test that a uniform speed value is preserved
    """
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.25
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi)

    d, f_ext     = extension_velocities(phi, speed, dx, self_test=True)
    np.testing.assert_allclose(f_ext, 1.0)


def test_exntension1():
    """ constant value marchout test
    """
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.25
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi)
    speed[abs(Y)<0.3] = 10.0
    d, f_ext     = extension_velocities(phi, speed, dx, self_test=True)
    np.testing.assert_allclose(f_ext, 10.0)


def test_distance_from_extension():
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.25
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi)
    d, f_ext     = extension_velocities(phi, speed, dx, self_test=True)
    d2           = distance(phi, dx, self_test=True)
    np.testing.assert_allclose(d, d2)

def test_extension_glitch():
    N     = 150
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2
    speed = np.ones_like(phi)
    speed[X>0.25] = 3.0
    d2, f_ext = extension_velocities(phi, speed, dx)

    assert (f_ext <= 3.0000001).all()
    assert (f_ext >= 1).all()

    np.testing.assert_almost_equal(f_ext[137, 102], 1, 3)
    np.testing.assert_almost_equal(f_ext[103, 78], 1, 2)
    np.testing.assert_almost_equal(f_ext[72, 100], 3, 3)
    np.testing.assert_almost_equal(f_ext[72, 86], 3, 3)
    np.testing.assert_almost_equal(f_ext[110, 121], 3, 3)

