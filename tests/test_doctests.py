import doctest
import importlib
import skfmm


def test_skfmm_doctests():
    results = doctest.testmod(skfmm, verbose=False, optionflags=doctest.ELLIPSIS)
    assert results.failed == 0, f"{results.failed} doctest(s) failed in skfmm"


def test_heap_doctests():
    # skfmm/__init__.py does `from .heap import heap`, so `skfmm.heap` is the
    # class, not the submodule.  importlib gives us the module object directly.
    heap_module = importlib.import_module('skfmm.heap')
    results = doctest.testmod(heap_module, verbose=False)
    assert results.failed == 0, f"{results.failed} doctest(s) failed in skfmm.heap"
