import numpy as np
import pytest
from sys import float_info

from skfmm.bi_cubic_init import (ainv, bc_interp, bc_interp_eq2, _newton2d,
                                 _newton2d_vec, _cubic_roots01, BiCubicInit)
import skfmm


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_interp(seed):
    """bc_interp from a fixed random coefficient vector."""
    rng = np.random.default_rng(seed)
    return bc_interp(rng.standard_normal(16) * 0.5)


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
    return bc_interp(ainv @ X_data)


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
            candidate = bc_interp(a_flat)
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
        interp = bc_interp(np.zeros(16))
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


# ---------------------------------------------------------------------------
# Vectorized helpers
# ---------------------------------------------------------------------------

class TestVectorizedHelpers:

    def test_newton_vec_matches_scalar(self):
        rng = np.random.default_rng(3)
        n = 200
        a = rng.standard_normal((n, 16)) * 0.5
        b0 = rng.integers(0, 2, n).astype(float)
        b1 = rng.integers(0, 2, n).astype(float)
        x0 = rng.random(n)
        y0 = rng.random(n)
        # a[k].reshape((4,4), order="f") is the scalar layout; the batch
        # layout is the same matrix
        A = a.reshape(n, 4, 4).transpose(0, 2, 1)
        xv, yv, cv = _newton2d_vec(A, b0, b1, x0, y0)
        for k in range(n):
            eq2 = bc_interp_eq2(bc_interp(a[k]), b0[k], b1[k])
            xs, ys, cs = _newton2d(eq2, x0=x0[k], y0=y0[k])
            assert cs == cv[k]
            if cs:
                assert xs == pytest.approx(xv[k], abs=1e-9)
                assert ys == pytest.approx(yv[k], abs=1e-9)

    def test_cubic_roots_match_numpy(self):
        rng = np.random.default_rng(5)
        c = rng.standard_normal((500, 4))
        # include some cubics with three roots in [0, 1]
        r = rng.random((50, 3))
        c[:50] = np.array([-r[:, 0]*r[:, 1]*r[:, 2],
                           r[:, 0]*r[:, 1] + r[:, 1]*r[:, 2] + r[:, 0]*r[:, 2],
                           -r.sum(axis=1), np.ones(50)]).T
        got = _cubic_roots01(c)
        for k in range(c.shape[0]):
            z = np.roots(c[k, ::-1])
            want = np.sort([w.real for w in z
                            if abs(w.imag) < 1e-9 and 1e-9 < w.real < 1 - 1e-9])
            have = np.sort(got[k][~np.isnan(got[k])])
            assert have.size == want.size
            assert np.allclose(have, want, atol=1e-7)

    def test_cubic_roots_degenerate(self):
        c = np.array([[-0.5, 1.0, 0.0, 0.0],   # linear
                      [0.0, 0.0, 0.0, 0.0],    # identically zero
                      [1.0, 0.0, 0.0, 0.0],    # constant, no root
                      [-0.25, 0.0, 1.0, 0.0]]) # quadratic, root at 0.5
        got = _cubic_roots01(c)
        assert np.nanmin(got[0]) == pytest.approx(0.5)
        assert np.all(np.isnan(got[2]))
        assert np.nanmin(got[3]) == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Robustness: every frozen point gets a distance
# ---------------------------------------------------------------------------

def _frozen_ok(phi, speed=None):
    b = BiCubicInit(phi, 1, speed)  # raises if a frozen point is missed
    d = b.d[b.aborders]
    # the zero on an edge with a sign change is at most one cell away
    assert np.all((d >= 0) & (d <= 1.0 + 1e-12))
    return b


class TestRobustness:

    def test_circle_signed_distance(self):
        # a true SDF input; the original center-start Newton failed here
        N = 61
        x, h = np.linspace(-1.5, 1.5, N, retstep=True)
        X, Y = np.meshgrid(x, x)
        R = np.sqrt(X**2 + Y**2)
        b = _frozen_ok((R - 1.0) / h)
        m = b.aborders
        assert np.abs(b.d[m] - np.abs(R - 1.0)[m] / h).max() < 0.01

    @pytest.mark.parametrize("seed", range(10))
    def test_random_rough_fields(self, seed):
        # under-resolved contours with thin slivers and kinks
        rng = np.random.default_rng(seed)
        phi = rng.standard_normal((40, 40))
        _frozen_ok(phi, speed=rng.random((40, 40)) + 1)

    def test_square_corners(self):
        x = np.linspace(-1, 1, 37)
        X, Y = np.meshgrid(x, x)
        _frozen_ok(np.maximum(np.abs(X), np.abs(Y)) - 0.5)

    def test_thin_sliver(self):
        x = np.linspace(-1, 1, 50)
        X, Y = np.meshgrid(x, x)
        _frozen_ok(np.abs(Y - 0.3*X) - 0.01)


def test_distance_second_order_convergence():
    errs = []
    for N in (41, 81, 161):
        x, h = np.linspace(-1.5, 1.5, N, retstep=True)
        X, Y = np.meshgrid(x, x)
        exact = np.sqrt(X**2 + Y**2) - 1
        d = skfmm.distance(X**2 + Y**2 - 1, dx=float(h), init_order=2)
        errs.append(np.abs(d - exact)[np.abs(exact) < 0.3].max())
    assert errs[1] < errs[0] / 2.5 and errs[2] < errs[1] / 2.5

