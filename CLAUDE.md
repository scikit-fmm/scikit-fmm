# scikit-fmm — notes for Claude Code

This is an open source scientific computing Python package. There is a
c compiled module.

## Build and dev setup

Use meson-python for development installs — standard editable install doesn't work:

to do a build:
```
python -m pip install --no-build-isolation --force-reinstall -e .
```

Install test dependencies separately:
```
python -m pip install ".[test]"
```

linting Python:
```
ruff check skfmm/
```


## Running tests

With the editable dev install (`--no-build-isolation`), the compiled `cfmm` extension
is built in-place inside `skfmm/`, so pytest can be run from the repo root:

```
python -m pytest tests/
```

In CI the package is installed normally (`pip install ".[test]"`), which puts the
extension in site-packages. Running pytest from the repo root then causes the local
`skfmm/` directory to shadow the installed package on `sys.path` and `cfmm` won't
be found. CI therefore runs from the parent directory:

```
cd .. && python -m pytest scikit-fmm/tests/
```

`skfmm.test(verbose=True)` is a public API that invokes pytest internally — keep it
working when changing the test setup.

## Non-obvious code constraints

- `BiCubicInit` asserts `h==1`. It only works with unit grid spacing. Callers in
  `pfmm.py` scale the resulting distances by `dx[0]` afterward.

- `from .heap import heap` in `__init__.py` shadows the submodule name, so
  `skfmm.heap` resolves to the `heap` class, not the module. In tests, use
  `importlib.import_module('skfmm.heap')` to get the module object for doctests.

- `dinit_` is `protected` in `distanceMarcher` so derived marcher classes
  (`extensionVelocityMarcher`) can access it directly without a redundant copy.

## Public API conventions

- The public Python API lives in `pfmm.py`. The type stub `__init__.pyi` must be
  kept in sync whenever function signatures change.

- NumPy-style docstrings. Use `{1, 2}` notation for enumerated int parameters
  (e.g. `init_order : {1, 2}, optional`), `Literal[1, 2]` in the stub.

## `init_order=2` (experimental higher-order initialization)

The `init_order=2` path in `distance()`, `travel_time()`, and
`extension_velocities()` uses bicubic interpolation (`BiCubicInit`) to compute
more accurate initial distances for the narrow-band frozen points. Constraints:

- 2D arrays only
- Uniform grid spacing (`dx[0] == dx[1]`)
- `order=2` (first-order FMM update is incompatible)
- No periodic boundaries
- No masked arrays
- `extension_velocities` only: no `ext_mask`

Active development is on the `new-init` branch.

### How `BiCubicInit` finds the closest point (robustness design)

`bi_cubic_init.py` is vectorized: all (border cell, frozen corner) queries
are solved at once in NumPy (`_newton2d_vec`, `_cubic_roots01`,
`_edge_zero_points`); the scalar `bc_interp`/`bc_interp_eq2`/`_newton2d`
remain as the readable reference and are cross-checked against the
vectorized versions in the tests.

- Only cells whose four corners change sign are processed (`find_frozen`).
  The old loop also marked about as many cells with no sign change, where
  Newton cannot find a zero; those produced most of the old
  "Newton did not converge" messages.
- Candidates for a corner's closest point in a cell are (a) Newton
  solutions of the closest-point equations inside the cell and (b) the
  zeros of the bicubic on the cell's four edges. The minimum distance wins.
- (b) is the robustness guarantee: along an edge the bicubic is exactly the
  1D cubic Hermite interpolant, so an edge whose end values differ in sign
  always has a root (found by splitting at critical points + bisection,
  which cannot fail). Every frozen point is on such an edge, so the old
  "frozen points were not initialized" RuntimeError can no longer happen.
- Newton is damped (step <= 0.5 cell) and started from the first-order
  projection `q - phi*grad/|grad|^2` (clamped to the cell), then from the
  cell center only if that gives nothing usable. Newton converging
  *outside* the cell is normal (closest point is in a neighbor cell) and
  is not a failure. `n_newton_fail` counts true non-convergence
  (diagnostic only; no printing).
- Found with the salt.py dissolution example (C:/src/intrepid2026): the old
  code failed after ~300 reinitializations of a kinked/under-resolved
  front, and even on the plain circle SDF `sqrt(x^2+y^2)-1` at 61x61.
  Now: no failures; ~15x faster than the old Python loops; distance error
  converges at ~2nd order vs 1st order for `init_order=1`.

### Build gotcha on this machine

The editable install rebuilds the C++ extension whenever a .cpp/.h mtime
changes (e.g. after `git stash`/`git stash pop`), and the rebuild fails
with `CreateProcess failed` if `cl` isn't on PATH. Rebuild from an MSVC
environment:
`cmd /c "\"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat\" && cd /d build\cp314 && meson compile"`
