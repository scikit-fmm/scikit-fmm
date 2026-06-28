import numpy as np
import pytest
from sys import float_info

from skfmm.bi_cubic_init import ainv, bc_interp, bc_interp_eq2, _newton2d, BiCubicInit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_interp(seed):
    """bc_interp from a fixed random coefficient vector."""
    rng = np.random.default_rng(seed)
    a_flat = rng.standard_normal(16) * 0.5
    return bc_interp(np.matrix(a_flat).T)


def _linear_x_interp():
    """bc_interp for f(x,y) = x - 0.5.

    Zero contour is the vertical line x = 0.5.  All second and mixed
    derivatives are exactly zero; fx = 1, fy = 0 everywhere.
    """
    X_data = np.array([
        -0.5,  0.5, -0.5,  0.5,   # f at (0,0),(1,0),(0,1),(1,1)
         1.0,  1.0,  1.0,  1.0,   # gx = 1 everywhere
         0.0,  0.0,  0.0,  0.0,   # gy = 0 everywhere
         0.0,  0.0,  0.0,  0.0,   # fxy = 0 everywhere
    ])
    a = ainv * np.matrix(X_data).T
    return bc_interp(a)


# ---------------------------------------------------------------------------
# Derivative methods
# ---------------------------------------------------------------------------

class TestDerivativeMethods:
    """Analytic derivatives verified against centred finite differences."""

    def setup_method(self):
        self.interp = _make_interp(42)
        self.x, self.y = 0.3, 0.7

    def test_fx(self):
        x, y, h = self.x, self.y, 1e-6
        fd = (self.interp(x + h, y) - self.interp(x - h, y)) / (2 * h)
        assert abs(self.interp.fx(x, y) - fd) < 1e-9

    def test_fy(self):
        x, y, h = self.x, self.y, 1e-6
        fd = (self.interp(x, y + h) - self.interp(x, y - h)) / (2 * h)
        assert abs(self.interp.fy(x, y) - fd) < 1e-9

    def test_fxx(self):
        # h=1e-4 balances truncation and floating-point cancellation
        x, y, h = self.x, self.y, 1e-4
        fd = (self.interp(x + h, y) - 2 * self.interp(x, y) + self.interp(x - h, y)) / h**2
        assert abs(self.interp.fxx(x, y) - fd) < 1e-7

    def test_fyy(self):
        x, y, h = self.x, self.y, 1e-4
        fd = (self.interp(x, y + h) - 2 * self.interp(x, y) + self.interp(x, y - h)) / h**2
        assert abs(self.interp.fyy(x, y) - fd) < 1e-7

    def test_fxy(self):
        x, y, h = self.x, self.y, 1e-4
        fd = (self.interp(x+h, y+h) - self.interp(x+h, y-h)
              - self.interp(x-h, y+h) + self.interp(x-h, y-h)) / (4 * h**2)
        assert abs(self.interp.fxy(x, y) - fd) < 1e-7

    def test_linear_derivatives_exact(self):
        """For f(x,y)=x-0.5, all derivatives except fx should be zero."""
        interp = _linear_x_interp()
        for x, y in [(0.1, 0.2), (0.5, 0.5), (0.9, 0.8)]:
            assert interp.fx(x, y)  == pytest.approx(1.0, abs=1e-12)
            assert interp.fy(x, y)  == pytest.approx(0.0, abs=1e-12)
            assert interp.fxx(x, y) == pytest.approx(0.0, abs=1e-12)
            assert interp.fyy(x, y) == pytest.approx(0.0, abs=1e-12)
            assert interp.fxy(x, y) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# fj consistency
# ---------------------------------------------------------------------------

class TestFjConsistency:
    """bc_interp_eq2.fj must agree with individual __call__ / derivative methods."""

    def setup_method(self):
        interp = _make_interp(42)
        self.interp = interp
        self.eq2 = bc_interp_eq2(interp, b0=0.0, b1=0.0)
        self.x, self.y = 0.4, 0.6

    def test_f0_matches_interp_call(self):
        x, y = self.x, self.y
        f0, *_ = self.eq2.fj(x, y)
        assert f0 == pytest.approx(bc_interp.__call__(self.eq2, x, y))

    def test_f1_matches_eq2_call(self):
        x, y = self.x, self.y
        _, f1, *_ = self.eq2.fj(x, y)
        assert f1 == pytest.approx(self.eq2(x, y))

    def test_j00_is_fx(self):
        x, y = self.x, self.y
        _, _, j00, *_ = self.eq2.fj(x, y)
        assert j00 == pytest.approx(self.interp.fx(x, y))

    def test_j01_is_fy(self):
        x, y = self.x, self.y
        _, _, _, j01, *_ = self.eq2.fj(x, y)
        assert j01 == pytest.approx(self.interp.fy(x, y))


# ---------------------------------------------------------------------------
# Newton solver
# ---------------------------------------------------------------------------

class TestNewtonSolver:

    def test_linear_x_known_solution(self):
        """f(x,y) = x - 0.5: closest point to (0,0) is exactly (0.5, 0)."""
        interp = _linear_x_interp()
        eq2 = bc_interp_eq2(interp, b0=0.0, b1=0.0)
        sx, sy, converged = _newton2d(eq2)
        assert converged
        assert sx == pytest.approx(0.5, abs=1e-10)
        assert sy == pytest.approx(0.0, abs=1e-10)

    def test_residuals_small_at_convergence(self):
        """Converged solutions must satisfy both equations to tight tolerance."""
        # seed 7 produces a bicubic with a sign change in [0,1]^2
        rng = np.random.default_rng(7)
        interp = None
        for _ in range(200):
            a_flat = rng.standard_normal(16) * 0.5
            candidate = bc_interp(np.matrix(a_flat).T)
            vals = [candidate(x, y) for x, y in [(0, 0), (1, 0), (0, 1), (1, 1)]]
            if min(vals) < 0 < max(vals):
                interp = candidate
                break
        assert interp is not None, "setup: no sign change found"

        for b0, b1 in [(0, 0), (1, 0), (0, 1), (1, 1)]:
            eq2 = bc_interp_eq2(interp, b0, b1)
            sx, sy, converged = _newton2d(eq2)
            if converged:
                assert abs(interp(sx, sy)) < 1e-10
                assert abs(eq2(sx, sy))    < 1e-10

    def test_singular_jacobian_reports_not_converged(self):
        """f = 0 everywhere → Jacobian is zero → solver must not claim convergence."""
        interp = bc_interp(np.matrix(np.zeros(16)).T)
        eq2 = bc_interp_eq2(interp, b0=0.0, b1=0.0)
        _, _, converged = _newton2d(eq2)
        assert not converged


# ---------------------------------------------------------------------------
# BiCubicInit integration
# ---------------------------------------------------------------------------

class TestBiCubicInit:

    def test_circle_all_borders_initialized(self):
        N = 30
        x = np.linspace(-1.5, 1.5, N)
        X, Y = np.meshgrid(x, x)
        phi = X**2 + Y**2 - 1.0
        b = BiCubicInit(phi, 1)
        mask = b.d == float_info.max
        assert not np.any(b.aborders & mask), "some border points were not initialized"

    def test_exact_zero_phi_gets_zero_distance(self):
        # linspace(-1, 1, 21) includes x=0 exactly, so phi=X has grid points on the contour
        N = 21
        x = np.linspace(-1.0, 1.0, N)
        X, Y = np.meshgrid(x, x)
        phi = X.copy()
        b = BiCubicInit(phi, 1)
        zero_mask = phi == 0.0
        assert zero_mask.any(), "setup: no grid point exactly on zero contour"
        assert np.all(b.d[zero_mask] == 0.0)
