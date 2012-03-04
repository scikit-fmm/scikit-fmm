import numpy as np
from skfmm import extension_velocities


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

    d, f_ext     = extension_velocities(phi, dx, speed, self_test=True)
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
    d, f_ext     = extension_velocities(phi, dx, speed, self_test=True)
    np.testing.assert_allclose(f_ext, 10.0)

