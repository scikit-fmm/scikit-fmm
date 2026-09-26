import numpy as np
from sys import float_info

ainv = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                 [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                 [-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0],
                 [2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                 [0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0],
                 [0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0],
                 [-3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0],
                 [0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0],
                 [9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1],
                 [-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1],
                 [2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0],
                 [0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0],
                 [-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1],
                 [4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1]], dtype=np.double)

class bc_interp(object):
    def __init__(self, a):
        self.a = a.reshape((4,4),order="f")

    def __call__(self, x, y):
        a = self.a
        return (a[0,0] + a[0,1]*y + a[0,2]*y**2 + a[0,3]*y**3 +
                a[1,0]*x + a[1,1]*x*y + a[1,2]*x*y**2 + a[1,3]*x*y**3 +
                a[2,0]*x**2 + a[2,1]*x**2*y + a[2,2]*x**2*y**2 +
                a[2,3]*x**2*y**3 + a[3,0]*x**3 + a[3,1]*x**3*y +
                a[3,2]*x**3*y**2 + a[3,3]*x**3*y**3)

    def fx(self, x, y):
        a = self.a
        return (a[1,0] + a[1,1]*y + a[1,2]*y**2 + a[1,3]*y**3 +
                2*a[2,0]*x + 2*a[2,1]*x*y + 2*a[2,2]*x*y**2 + 2*a[2,3]*x*y**3 +
                3*a[3,0]*x**2 + 3*a[3,1]*x**2*y + 3*a[3,2]*x**2*y**2 +
                3*a[3,3]*x**2*y**3)

    def fy(self, x, y):
        a = self.a
        return (a[0,1] + 2*a[0,2]*y + 3*a[0,3]*y**2 +
                a[1,1]*x + 2*a[1,2]*x*y + 3*a[1,3]*x*y**2 +
                a[2,1]*x**2 + 2*a[2,2]*x**2*y + 3*a[2,3]*x**2*y**2 +
                a[3,1]*x**3 + 2*a[3,2]*x**3*y + 3*a[3,3]*x**3*y**2)

    def fxx(self, x, y):
        a = self.a
        return (2*a[2,0] + 2*a[2,1]*y + 2*a[2,2]*y**2 + 2*a[2,3]*y**3 +
                6*a[3,0]*x + 6*a[3,1]*x*y + 6*a[3,2]*x*y**2 + 6*a[3,3]*x*y**3)

    def fyy(self, x, y):
        a = self.a
        return (2*a[0,2] + 6*a[0,3]*y +
                2*a[1,2]*x + 6*a[1,3]*x*y +
                2*a[2,2]*x**2 + 6*a[2,3]*x**2*y +
                2*a[3,2]*x**3 + 6*a[3,3]*x**3*y)

    def fxy(self, x, y):
        a = self.a
        return (a[1,1] + 2*a[1,2]*y + 3*a[1,3]*y**2 +
                2*a[2,1]*x + 4*a[2,2]*x*y + 6*a[2,3]*x*y**2 +
                3*a[3,1]*x**2 + 6*a[3,2]*x**2*y + 9*a[3,3]*x**2*y**2)


class bc_interp_eq2(bc_interp):
    def __init__(self, interp, b0, b1):
        self.a = interp.a
        self.b0, self.b1 = b0, b1

    def __call__(self, x, y):
        return -(self.b0 - x) * self.fy(x, y) + (self.b1 - y) * self.fx(x, y)

    def fj(self, x, y):
        """Returns (F0, F1, J00, J01, J10, J11) in one pass for Newton-Raphson.

        F0 = f(x,y)  [bicubic value]
        F1 = eq2(x,y) [orthogonality condition: closest-point constraint]
        J is the 2x2 Jacobian of (F0, F1) w.r.t. (x, y), computed analytically.
        """
        b0, b1 = self.b0, self.b1
        _fx  = self.fx(x, y)
        _fy  = self.fy(x, y)
        _fxx = self.fxx(x, y)
        _fyy = self.fyy(x, y)
        _fxy = self.fxy(x, y)
        f0 = bc_interp.__call__(self, x, y)
        f1 = -(b0 - x) * _fy + (b1 - y) * _fx
        # J[0,:] = grad f
        j00 = _fx
        j01 = _fy
        # J[1,:] = grad eq2, derived by product rule
        j10 = _fy - (b0 - x) * _fxy + (b1 - y) * _fxx
        j11 = -_fx - (b0 - x) * _fyy + (b1 - y) * _fxy
        return f0, f1, j00, j01, j10, j11


def _cross_derivative(f, h):
    """Cross-derivative (fxy) of a 2D field using the 4-point corner stencil."""
    xgr = np.zeros_like(f)
    nx, ny = f.shape
    xgr[1:-1,1:-1] = (f[2:,2:] - f[2:,:-2] - f[:-2,2:] + f[:-2,:-2])
    xgr[   0,   0] = -(f[ 2, 0] - f[ 2, 2] - f[ 0, 0] + f[ 0, 2])
    xgr[nx-1,ny-1] = -(f[-1,-3] - f[-1,-1] - f[-3,-3] + f[-3,-1])
    xgr[   0,ny-1] =  (f[ 0,-3] - f[ 0,-1] - f[ 2,-3] + f[ 2,-1])
    xgr[nx-1,   0] = -(f[-1, 0] - f[-1, 2] - f[-3, 0] + f[-3, 2])
    xgr[0,   1:-1] = -(f[ 2, :-2] - f[ 2, 2:] - f[0, :-2] + f[0, 2:])
    xgr[nx-1,1:-1] =  (f[ -3, :-2] - f[ -3, 2:] - f[-1, :-2] + f[-1, 2:])
    xgr[1:-1,   0] = -(f[ 2:,  0] - f[2:, 2] - f[:-2,   0] + f[:-2, 2])
    xgr[1:-1,ny-1] = -(f[ 2:, -3] - f[2:, -1] - f[:-2, -3] + f[:-2, -1])
    xgr /= (4 * h**2)
    return xgr


def _newton2d(eq2, tol=1.49e-8, max_iter=50, x0=0.5, y0=0.5, max_step=0.5):
    """Damped Newton-Raphson solver for the 2D closest-point system.

    Finds (x, y) such that f(x,y)=0 and the vector to the query point
    is parallel to grad f (i.e. the closest point on the zero contour).
    Uses the analytic Jacobian from eq2.fj; no external dependencies.
    Each step is limited to length ``max_step`` (in cell units) and the
    iteration is abandoned if it wanders more than one cell away from the
    reference cell, since a solution there is not useful to the caller.
    Returns (x, y, converged).
    """
    x, y = x0, y0
    for _ in range(max_iter):
        f0, f1, j00, j01, j10, j11 = eq2.fj(x, y)
        det = j00 * j11 - j01 * j10
        if abs(det) < 1e-15:
            return x, y, False
        # 2x2 solve via Cramer's rule
        dx = (j11 * f0 - j01 * f1) / det
        dy = (j00 * f1 - j10 * f0) / det
        step = np.hypot(dx, dy)
        if step > max_step:
            dx *= max_step / step
            dy *= max_step / step
        x -= dx
        y -= dy
        if abs(dx) < tol and abs(dy) < tol:
            return x, y, True
        if not (-1.0 < x < 2.0 and -1.0 < y < 2.0):
            return x, y, False
    return x, y, False


def _powers(t):
    """Rows (1, t, t^2, t^3) and their first and second t-derivatives."""
    one, zero = np.ones_like(t), np.zeros_like(t)
    p = np.stack((one, t, t*t, t*t*t), axis=-1)
    dp = np.stack((zero, one, 2*t, 3*t*t), axis=-1)
    ddp = np.stack((zero, zero, 2*one, 6*t), axis=-1)
    return p, dp, ddp


def _bicubic(a, px, py):
    """Evaluate sum a[k,p,q] * px[k,p] * py[k,q] for a batch of cells."""
    return np.einsum('kpq,kp,kq->k', a, px, py)


def _newton2d_vec(a, b0, b1, x, y, tol=1.49e-8, max_iter=50, max_step=0.5):
    """Vectorized version of ``_newton2d`` over a batch of problems.

    ``a`` has shape (n, 4, 4) (bicubic coefficients, ``a[k, p, q]`` multiplies
    x**p * y**q); ``b0``, ``b1`` are the query points and ``x``, ``y`` the
    starting points, all shape (n,). Same iteration, damping and stopping
    rules as ``_newton2d``. Returns (x, y, converged).
    """
    x = np.array(x, dtype=np.double)
    y = np.array(y, dtype=np.double)
    converged = np.zeros(x.shape, dtype=bool)
    active = np.ones(x.shape, dtype=bool)
    for _ in range(max_iter):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        A, xi, yi, q0, q1 = a[idx], x[idx], y[idx], b0[idx], b1[idx]
        px, dpx, ddpx = _powers(xi)
        py, dpy, ddpy = _powers(yi)
        f = _bicubic(A, px, py)
        fx = _bicubic(A, dpx, py)
        fy = _bicubic(A, px, dpy)
        fxx = _bicubic(A, ddpx, py)
        fyy = _bicubic(A, px, ddpy)
        fxy = _bicubic(A, dpx, dpy)
        f1 = -(q0 - xi) * fy + (q1 - yi) * fx
        j00, j01 = fx, fy
        j10 = fy - (q0 - xi) * fxy + (q1 - yi) * fxx
        j11 = -fx - (q0 - xi) * fyy + (q1 - yi) * fxy
        det = j00 * j11 - j01 * j10
        singular = np.abs(det) < 1e-15
        det = np.where(singular, 1.0, det)
        dx = (j11 * f - j01 * f1) / det
        dy = (j00 * f1 - j10 * f) / det
        step = np.hypot(dx, dy)
        scale = np.where(step > max_step,
                         max_step / np.maximum(step, 1e-300), 1.0)
        dx = np.where(singular, 0.0, dx * scale)
        dy = np.where(singular, 0.0, dy * scale)
        xi, yi = xi - dx, yi - dy
        x[idx], y[idx] = xi, yi
        done = ~singular & (np.abs(dx) < tol) & (np.abs(dy) < tol)
        converged[idx[done]] = True
        escaped = ~((-1.0 < xi) & (xi < 2.0) & (-1.0 < yi) & (yi < 2.0))
        active[idx[singular | done | escaped]] = False
    return x, y, converged


def _cubic_roots01(c, n_bisect=60):
    """All real roots in [0, 1] of the cubics c0 + c1 s + c2 s^2 + c3 s^3.

    ``c`` has shape (n, 4). Returns an (n, 3) array of roots, NaN where a
    cubic has fewer than three. The interval is split at the cubic's
    critical points into (at most three) monotone pieces, and each piece
    whose end values differ in sign (or touch zero) holds exactly one root,
    found by bisection -- no convergence failure is possible.
    """
    c = np.asarray(c, dtype=np.double)
    n = c.shape[0]
    c0, c1, c2, c3 = c.T

    def g(s):
        return ((c3[:, None] * s + c2[:, None]) * s + c1[:, None]) * s + \
            c0[:, None]

    # critical points: roots of 3 c3 s^2 + 2 c2 s + c1
    qa, qb, qc = 3 * c3, 2 * c2, c1
    crit = np.full((n, 2), np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        disc = qb * qb - 4 * qa * qc
        quad = (qa != 0) & (disc >= 0)
        q = -0.5 * (qb + np.copysign(np.sqrt(np.maximum(disc, 0.0)), qb))
        crit[:, 0] = np.where(quad, q / qa, np.nan)
        crit[:, 1] = np.where(quad, qc / q, np.nan)
        lin = (qa == 0) & (qb != 0)
        crit[lin, 0] = -qc[lin] / qb[lin]
    crit[~((crit > 0) & (crit < 1))] = 1.0  # NaN or outside -> empty piece
    knots = np.sort(np.column_stack((np.zeros(n), crit, np.ones(n))), axis=1)
    lo, hi = knots[:, :-1].copy(), knots[:, 1:].copy()
    glo = g(lo)
    has_root = glo * g(hi) <= 0
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        gmid = g(mid)
        left = glo * gmid <= 0
        hi = np.where(left, mid, hi)
        lo = np.where(left, lo, mid)
        glo = np.where(left, glo, gmid)
    return np.where(has_root, 0.5 * (lo + hi), np.nan)


def _edge_zero_points(a):
    """Points on the boundary of each unit reference cell where the bicubic
    is zero.

    ``a`` has shape (n, 4, 4). Returns (x, y), each shape (n, 12), NaN where
    there is no point. Along each cell edge the bicubic reduces to a cubic
    (the 1D cubic Hermite interpolant of phi and its derivative along that
    edge), so an edge whose end values differ in sign always has a root.
    These points lie on the interpolated zero contour, so they are always
    valid closest-point candidates; they are the fallback that guarantees
    every point next to a sign change gets a distance.
    """
    n = a.shape[0]
    zero, one = np.zeros((n, 3)), np.ones((n, 3))
    s_y0 = _cubic_roots01(a[:, :, 0])        # edge y=0: (s, 0)
    s_y1 = _cubic_roots01(a.sum(axis=2))     # edge y=1: (s, 1)
    s_x0 = _cubic_roots01(a[:, 0, :])        # edge x=0: (0, s)
    s_x1 = _cubic_roots01(a.sum(axis=1))     # edge x=1: (1, s)
    ex = np.hstack((s_y0, s_y1, zero, one))
    ey = np.hstack((zero, one, s_x0, s_x1))
    nan = np.isnan(np.hstack((s_y0, s_y1, s_x0, s_x1)))
    ex[nan] = np.nan
    ey[nan] = np.nan
    return ex, ey


def _cell_coefficients(f, gx, gy, gxy, ci, cj):
    """Bicubic coefficients, shape (n, 4, 4), for the cells whose lower-left
    corners are (ci, cj); ``a[k, p, q]`` multiplies x**p * y**q."""
    ii = np.stack((ci, ci+1, ci, ci+1), axis=1)
    jj = np.stack((cj, cj, cj+1, cj+1), axis=1)
    X = np.hstack((f[ii, jj], gx[ii, jj], gy[ii, jj], gxy[ii, jj]))
    return (X @ ainv.T).reshape(-1, 4, 4).transpose(0, 2, 1)


class BiCubicInit(object):
    """Second-order initialization of the frozen points next to the zero
    contour of a 2D level set function on a unit-spaced grid.

    For every cell the zero contour crosses, phi is interpolated with a
    bicubic (values, gradients and cross-derivative at the four corners),
    and for each frozen corner of the cell the distance to the closest point
    on the interpolated zero contour is computed. A frozen point's distance
    is the minimum over its adjacent cells. The closest point in a cell is
    the nearest of

    - interior solutions of the closest-point equations (damped Newton,
      started from the first-order projection of the query point, then
      from the cell center if that fails), and
    - the zero points on the cell's edges (exact cubic roots).

    The edge points make the method robust: every frozen point is next to
    an edge with a sign change, which always holds a root, so every frozen
    point gets a distance even where Newton fails (under-resolved or kinked
    contours).

    Attributes: ``d`` (distance, ``float_info.max`` where not frozen),
    ``aborders`` (frozen points), ``border_cells`` (processed cells),
    ``pdict`` (offset to the closest point for each frozen point),
    ``n_newton_fail`` (diagnostic count of (point, cell) pairs where Newton
    never converged; harmless since the edge points cover them), and with
    ``speed``: ``out_speed``/``_speed_set`` (speed interpolated at the
    closest point).
    """

    def __init__(self, phi, h, speed=None):
        self.in_speed = None
        self.out_speed = None
        phi = phi.astype(np.double)
        if speed is not None:
            self.in_speed = speed.astype(np.double)
            assert self.in_speed.shape == phi.shape
            self.out_speed = np.zeros_like(self.in_speed)
            self._speed_set = np.zeros(phi.shape, dtype=bool)

        assert h==1
        self.phi = phi
        self.d = np.ones_like(phi) * float_info.max
        self.h = h
        self.pdict = {}
        self.n_newton_fail = 0
        assert len(phi.shape)==2

        gx, gy = np.gradient(phi, edge_order=2)
        self.gx, self.gy = gx/h, gy/h # is this correct for h != 1?
        self.xgr = _cross_derivative(phi, h)

        if self.in_speed is not None:
            gx_s, gy_s = np.gradient(self.in_speed, edge_order=2)
            self.gx_speed = gx_s / h
            self.gy_speed = gy_s / h
            self.xgr_speed = _cross_derivative(self.in_speed, h)
        self.find_frozen()
        self.ini_frozen()

        if np.any(self.aborders & (self.d == float_info.max)):
            raise RuntimeError("scikit-fmm (skfmm) error, some initially "
                               "frozen points were not initialized correctly.")

    def find_frozen(self):
        """Mark the initially frozen points and the cells to process.

        A point is frozen (``aborders``) if phi is exactly zero there or
        phi changes sign between it and one of its four neighbors; this is
        the same set the first-order C++ initialization freezes.  A cell
        (indexed by its lower-left corner) is processed if its four corner
        values change sign, i.e. the zero contour crosses it.  Every edge
        with a sign change belongs to such a cell.
        """
        phi = self.phi
        aborders = phi == 0.0
        self.d[aborders] = 0.0
        sx = phi[:-1, :] * phi[1:, :] < 0
        sy = phi[:, :-1] * phi[:, 1:] < 0
        aborders[:-1, :] |= sx
        aborders[1:, :] |= sx
        aborders[:, :-1] |= sy
        aborders[:, 1:] |= sy

        corners = np.stack((phi[:-1, :-1], phi[1:, :-1],
                            phi[:-1, 1:], phi[1:, 1:]))
        border_cells = np.zeros_like(phi, dtype=bool)
        border_cells[:-1, :-1] = ((corners.min(axis=0) < 0) &
                                  (corners.max(axis=0) > 0))
        self.aborders = aborders
        self.border_cells = border_cells

    def ini_frozen(self):
        phi = self.phi
        nx, ny = phi.shape

        # points (numerically) on the contour: distance zero
        on_contour = self.aborders & (np.abs(phi) < float_info.epsilon)
        self.d[on_contour] = 0.0
        if self.in_speed is not None:
            self.out_speed[on_contour] = self.in_speed[on_contour]
            self._speed_set[on_contour] = True

        ci, cj = np.nonzero(self.border_cells)
        if ci.size == 0:
            return
        a = _cell_coefficients(phi, self.gx, self.gy, self.xgr, ci, cj)

        # one query per (cell, frozen corner) pair
        cell = np.repeat(np.arange(ci.size), 4)
        qi = np.tile(np.array([0, 1, 0, 1]), ci.size)  # corner in cell coords
        qj = np.tile(np.array([0, 0, 1, 1]), ci.size)
        pi, pj = ci[cell] + qi, cj[cell] + qj           # corner in the grid
        keep = self.aborders[pi, pj] & ~on_contour[pi, pj]
        cell, qi, qj, pi, pj = cell[keep], qi[keep], qj[keep], pi[keep], pj[keep]
        qi = qi.astype(np.double)
        qj = qj.astype(np.double)

        # candidates on the cell edges
        ex, ey = _edge_zero_points(a)
        edist = np.hypot(ex[cell] - qi[:, None], ey[cell] - qj[:, None])
        edist[np.isnan(edist)] = np.inf
        k = np.argmin(edist, axis=1)
        rows = np.arange(cell.size)
        best = edist[rows, k]
        bx, by = ex[cell, k], ey[cell, k]

        # interior candidates from Newton: first from the first-order
        # projection of the query point onto the contour, then (for queries
        # that produced nothing usable) from the cell center
        g0, g1 = self.gx[pi, pj], self.gy[pi, pj]
        gg = g0*g0 + g1*g1
        pr = np.where(gg > 0, phi[pi, pj] / np.where(gg > 0, gg, 1.0), 0.0)
        x0 = np.where(gg > 0, np.clip(qi - pr*g0, 0.0, 1.0), 0.5)
        y0 = np.where(gg > 0, np.clip(qj - pr*g1, 0.0, 1.0), 0.5)
        on_boundary = (pi == 0) | (pj == 0) | (pi == nx-1) | (pj == ny-1)
        todo = rows
        any_converged = np.zeros(cell.size, dtype=bool)
        for start in ((x0, y0), (0.5, 0.5)):
            sx0 = np.broadcast_to(start[0], rows.shape)[todo]
            sy0 = np.broadcast_to(start[1], rows.shape)[todo]
            sx, sy, conv = _newton2d_vec(a[cell[todo]], qi[todo], qj[todo],
                                         sx0, sy0)
            any_converged[todo[conv]] = True
            inside = (sx >= 0) & (sx <= 1) & (sy >= 0) & (sy <= 1)
            # for boundary points we need to accept points outside [0,1]
            # this could/should be done more cleverly/safely
            ok = conv & (inside | on_boundary[todo])
            dist = np.hypot(sx - qi[todo], sy - qj[todo])
            better = ok & (dist < best[todo])
            t = todo[better]
            best[t], bx[t], by[t] = dist[better], sx[better], sy[better]
            # converged outside this cell: the closest point is in a
            # neighboring cell, which handles it; retry only the rest
            todo = todo[~ok]
        self.n_newton_fail = int(np.count_nonzero(~any_converged))

        # reduce over the cells adjacent to each point: keep the minimum
        found = np.isfinite(best)
        order = np.lexsort((best, pj, pi))
        order = order[found[order]]
        flat = pi[order] * ny + pj[order]
        _, first = np.unique(flat, return_index=True)
        w = order[first]
        wi, wj = pi[w], pj[w]
        self.d[wi, wj] = np.minimum(self.d[wi, wj], best[w])
        for i, j, ox, oy in zip(wi, wj, bx[w] - qi[w], by[w] - qj[w]):
            self.pdict[(int(i), int(j))] = (ox, oy)
        if self.in_speed is not None:
            a_s = _cell_coefficients(self.in_speed, self.gx_speed,
                                     self.gy_speed, self.xgr_speed,
                                     ci[cell[w]], cj[cell[w]])
            px, _, _ = _powers(bx[w])
            py, _, _ = _powers(by[w])
            self.out_speed[wi, wj] = _bicubic(a_s, px, py)
            self._speed_set[wi, wj] = True


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    N = 100
    x, h = np.linspace(-2.5,2.5,N,retstep=True)
    extent = (x[0],x[-1],x[0],x[-1])
    X,Y = np.meshgrid(x,x)

    # solve the eikonal equation: |grad phi| = 1
    # given elliptic zero contour 0.5*x^2 + y^2 = 1

    phi = X**2 + Y**2 - 1.
    plt.matshow(phi)
    plt.colorbar()
    plt.show()

    a=BiCubicInit(phi, 1)

    plt.matshow(a.aborders)
    plt.colorbar()
    plt.show()

    arr=np.copy(a.d)
    mask = arr==float_info.max
    bb=np.ma.MaskedArray(arr,mask)
    bb[phi<0] *= -1

    # make sure that all point adjacent to the zero contour got a value
    assert np.logical_and(a.aborders==True, np.logical_not(mask)).sum() == \
        (a.aborders==True).sum()

    plt.matshow(bb)
    plt.colorbar()
    plt.show()
