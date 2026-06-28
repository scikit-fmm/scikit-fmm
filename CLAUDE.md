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
