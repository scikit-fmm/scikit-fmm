import numpy as np
from skfmm import distance, travel_time

def test_issue_1():
    phi = np.arange(-5, 5) + 0.499
    d = distance(phi)
    t = travel_time(phi, speed=np.ones_like(phi))
    np.testing.assert_allclose(t, np.abs(d))
