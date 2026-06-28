from collections.abc import Sequence
from typing import Literal, overload

import numpy as np
from numpy.typing import ArrayLike

__version__: str
__docformat__: str

# dx: scalar applied uniformly, or one value per dimension
_Dx = float | int | Sequence[float]

# periodic: True/False for all dims, or per-dim sequence
_Periodic = bool | Sequence[bool]

# ---------------------------------------------------------------------------
# distance
# ---------------------------------------------------------------------------

@overload
def distance(
    phi: np.ma.MaskedArray,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
) -> np.ma.MaskedArray: ...

@overload
def distance(
    phi: ArrayLike,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
) -> np.ndarray | np.ma.MaskedArray: ...

# ---------------------------------------------------------------------------
# travel_time
# ---------------------------------------------------------------------------
# Returns MaskedArray when phi is masked.
# Also returns MaskedArray when any speed value is below machine epsilon
# (those cells are masked internally) or when narrow leaves far-field points.

@overload
def travel_time(
    phi: np.ma.MaskedArray,
    speed: ArrayLike,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
) -> np.ma.MaskedArray: ...

@overload
def travel_time(
    phi: ArrayLike,
    speed: ArrayLike,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
) -> np.ndarray | np.ma.MaskedArray: ...

# ---------------------------------------------------------------------------
# extension_velocities
# ---------------------------------------------------------------------------
# Returns a (distance, f_ext) tuple; each array is independently masked.
# When phi is masked both outputs are guaranteed to be MaskedArray.

@overload
def extension_velocities(
    phi: np.ma.MaskedArray,
    speed: ArrayLike,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    ext_mask: ArrayLike | None = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
    init_order: Literal[1, 2] = ...,
) -> tuple[np.ma.MaskedArray, np.ma.MaskedArray]: ...

@overload
def extension_velocities(
    phi: ArrayLike,
    speed: ArrayLike,
    dx: _Dx = ...,
    self_test: bool = ...,
    order: Literal[1, 2] = ...,
    ext_mask: ArrayLike | None = ...,
    narrow: float = ...,
    periodic: _Periodic = ...,
    init_order: Literal[1, 2] = ...,
) -> tuple[np.ndarray | np.ma.MaskedArray, np.ndarray | np.ma.MaskedArray]: ...

# ---------------------------------------------------------------------------
# heap
# ---------------------------------------------------------------------------

class heap:
    def __init__(self, max_size: int, self_test: bool = ...) -> None: ...
    def push(self, addr: int, value: float) -> int: ...
    def pop(self) -> tuple[int, float]: ...
    def update(self, heap_id: int, value: float) -> None: ...
    def empty(self) -> bool: ...
    def peek(self) -> float: ...

# ---------------------------------------------------------------------------
# test utilities
# ---------------------------------------------------------------------------

def testing() -> None: ...
def test(verbose: bool | None = ...) -> int: ...
