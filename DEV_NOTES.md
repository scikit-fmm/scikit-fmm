
## Development set up

To do a normal build from source see the README.md, if you are
developing and want to make changes and test them in this folder use
this approach:

- `pip install meson-python ninja numpy pipx`
Then build with
- `python -m pip install --no-build-isolation -e .`
This builds the module in place and sets up import hooks that rebuilds
the module on demand if you change the c++ files.

- force rebuild: `python -m pip install --no-build-isolation --force-reinstall -e .`
- linting `ruff check skfmm/`
- C++ linting: `pipx run clang-tidy skfmm/*.cpp -- -std=c++17`

## To make a release:

### Increment version number

Currently the version number is hard-coded in the following locations:

- `./doc/conf.py`
- `./pyproject.toml`
- `./meson.build`
- `./skfmm/__init__.py`
- `./README.txt`

The version number format is YYYY.MM.DD with leading zeros

### Tag the release

```
git tag YYYY.MM.DD
git push --tags origin master
```

### Build the documentation

To build the documentation you need Sphinx and the numpydoc Sphinx extension

```
make html
```

Check that everything looks ok by opening `build/html/index.html`

The documentation should build automatically on RTD.

Check https://scikit-fmm.readthedocs.io/en/latest/ — trigger a build if needed.

### Upload to PyPI

Make sure you have twine installed with:

```
pip install twine
```

Once the final commit has been made, wait for the GitHub workflow builds to finish. Make sure all the tests pass on all the platforms.

Copy the build artifacts from the GitHub actions to `dist/`

```
twine upload dist/*
```

## Cython Wrapper for heap class

The binary min heap C++ class has a Cython wrapper to expose it to Python. Cython is not required to build this module as the c++ files created by Cython are checked into the git repo. If changes are made to `skfmm/heap.h` or `skfmm/heap.cpp` the Cython c++ file `pheap.cpp` needs to be updated. This is done with the command:

```
cython3 --cplus -3 pheap.pyx
```

The updated `pheap.cpp` should also be checked into the git repo.
