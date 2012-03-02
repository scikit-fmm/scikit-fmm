from skfmm import distance, travel_time
from skfmm.utility import contour_lengths_2d, total_contour_length_2d
from skfmm.utility import contour_area_2d, normal_direction, curvature

import numpy as np

def test_contour_length():
    """
    """
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2

    length_array = contour_lengths_2d(X, Y, phi)
    length = total_contour_length_2d(X, Y, phi)
    assert length == length_array.sum()

    exact = 2.0 * np.pi * r
    np.testing.assert_allclose(length, exact, atol= dx / 4.0)

def test_contour_area():
    N     = 50
    X, Y  = np.meshgrid(np.linspace(-1, 1, N), np.linspace(-1, 1, N))
    r     = 0.5
    dx    = 2.0 / (N - 1)
    phi   = (X) ** 2 + (Y) ** 2 - r ** 2

    area = contour_area_2d(phi, dx)
    exact = np.pi * r ** 2.0
    np.testing.assert_allclose(area, exact, atol = dx / 4.0)

