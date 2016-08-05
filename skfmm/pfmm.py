from sys import float_info
import numpy as np

from .cfmm import cFastMarcher

FAR, NARROW, FROZEN, MASK = 0, 1, 2, 3
DISTANCE, TRAVEL_TIME, EXTENSION_VELOCITY = 0, 1, 2


def pre_process_args(phi, dx, narrow, periodic, ext_mask=None):
    """
    get input data into the correct form for calling the c extension module
    This wrapper allows for a little bit of flexibility in the input types
    """
    if not isinstance(phi, np.ndarray):
        phi = np.array(phi)

    if type(dx) is float or type(dx) is int:
        dx = [dx for x in range(phi.ndim)]
    dx = np.array(dx)

    if isinstance(phi, np.ma.MaskedArray):
        flag = np.zeros(phi.shape, dtype=np.int)
        flag[phi.mask] = MASK
        phi = phi.data
    else:
        flag = np.zeros(phi.shape, dtype=np.int)

    if ext_mask is None:
        ext_mask = np.zeros(phi.shape, dtype=np.int)

    periodic_data = 0
    if isinstance(periodic, bool):
        if periodic:
            periodic_data = int(2**phi.ndim-1)
    else:
        if hasattr(periodic, "__len__") and len(periodic) == phi.ndim:
            for i, value in enumerate(periodic):
                if value:
                    periodic_data |= 1<<i
        else:
            raise ValueError("parameter \"periodic\" must be of type bool or sequence of type bool of length phi.ndim.")

    if narrow < 0:
        raise ValueError("parameter \"narrow\" must be greater than or equal to zero.")

    return phi, dx, flag, ext_mask, periodic_data


def post_process_result(result):
    """
    post-process results from the c module (add mask)
    """
    if (result == float_info.max).any():
        mask = (result == float_info.max)
        result[mask] = 0
        result = np.ma.MaskedArray(result, mask)
    return result


def distance(phi, dx=1.0, self_test=False, order=2,
             narrow=0.0, periodic=False):
    """Return the signed distance from the zero contour of the array phi.

    Parameters
    ----------
    phi : array-like
          the zero contour of this array is the boundary location for
          the distance calculation. Phi can of 1,2,3 or higher
          dimension and can be a masked array.

    dx  : float or an array-like of len phi.ndim, optional
          the cell length in each dimension.

    self_test : bool, optional
                if True consistency checks are made on the binary min
                heap during the calculation. This is used in testing and
                results in a slower calculation.

    order : int, optional
            order of computational stencil to use in updating points during
            the fast marching method. Must be 1 or 2, the default is 2.

    narrow : float, optional
             narrow band half-width. If this optional argument is
             specified the marching algorithm is limited to within a
             given narrow band. If far-field points remain when this
             condition is met a masked array is returned. The default
             value is 0.0 which means no narrow band limit.

    periodic : bool or an array-like of len phi.ndim, optional
               specifies whether and in which directions periodic
               boundary conditions are used. True sets periodic
               boundary conditions in all directions. An array-like
               (interpreted as True or False values) specifies the
               absence or presence of periodic boundaries in
               individual directions. The default value is False,
               i.e., no periodic boundaries in any direction.

    Returns
    -------
    d : an array the same shape as phi
        contains the signed distance from the zero contour (zero level set)
        of phi to each point in the array. The sign is specified by the sign
        of phi at the given point.

    """
    phi, dx, flag, ext_mask, periodic = \
                        pre_process_args(phi, dx, narrow, periodic)
    d = cFastMarcher(phi, dx, flag, None, ext_mask,
                     int(self_test), DISTANCE, order, narrow, periodic)
    d = post_process_result(d)
    return d


def travel_time(phi, speed, dx=1.0, self_test=False, order=2,
                narrow=0.0, periodic=False):
    """Return the travel from the zero contour of the array phi given the
    scalar velocity field speed.

    Parameters
    ----------
    phi : array-like
          the zero contour of this array is the boundary location for
          the travel time calculation. Phi can of 1,2,3 or higher
          dimension and can be a masked array.

    speed : array-like, the same shape as phi
            contains the speed of interface propagation at each point
            in the domain.

    dx  : float or an array-like of len phi.ndim, optional
          the cell length in each dimension.

    self_test : bool, optional
                if True consistency checks are made on the binary min
                heap during the calculation. This is used in testing and
                results in a slower calculation.

    order : int, optional
            order of computational stencil to use in updating points during
            the fast marching method. Must be 1 or 2, the default is 2.

    narrow : float, optional
             narrow band half-width. If this optional argument is
             specified the marching algorithm is limited to travel
             times within a given value. If far-field points
             remain when this condition is met a masked array is
             returned. The default value is 0.0 which means no narrow
             band limit.

    periodic : bool or an array-like of len phi.ndim, optional
               specifies whether and in which directions periodic
               boundary conditions are used. True sets periodic
               boundary conditions in all directions. An array-like
               (interpreted as True or False values) specifies the
               absence or presence of periodic boundaries in
               individual directions. The default value is False,
               i.e., no periodic boundaries in any direction.

    Returns
    -------
    t : an array the same shape as phi
        contains the travel time from the zero contour (zero level
        set) of phi to each point in the array given the scalar
        velocity field speed. If the input array speed has values less
        than or equal to zero the return value will be a masked array.

    """
    phi, dx, flag, ext_mask, periodic \
        = pre_process_args(phi, dx, narrow, periodic)
    t = cFastMarcher(phi, dx, flag, speed, ext_mask,
                     int(self_test), TRAVEL_TIME, order, narrow, periodic)
    t = post_process_result(t)
    return t


def extension_velocities(phi, speed, dx=1.0, self_test=False, order=2,
                         ext_mask=None, narrow=0.0, periodic=False):
    """Extend the velocities defined at the zero contour of phi, in the
    normal direction, to the rest of the domain. Extend the velocities
    such that grad f_ext dot grad d = 0 where where f_ext is the
    extension velocity and d is the signed distance function.

    Parameters
    ----------
    phi : array-like
          the zero contour of this array is the boundary location for
          the travel time calculation. Phi can of 1,2,3 or higher
          dimension and can be a masked array.

    speed : array-like, the same shape as phi
            contains the speed of interface propagation at each point
            in the domain.

    dx  : float or an array-like of len phi.ndim, optional
          the cell length in each dimension.

    self_test : bool, optional
                if True consistency checks are made on the binary min
                heap during the calculation. This is used in testing and
                results in a slower calculation.

    order : int, optional
            order of computational stencil to use in updating points during
            the fast marching method. Must be 1 or 2, the default is 2.

    ext_mask : array-like, the same shape as phi, optional
               enables initial front values to be eliminated when
               calculating the value at the interface before the
               values are extended away from the interface.

    narrow : float, optional
             narrow band half-width. If this optional argument is
             specified the marching algorithm is limited to within a
             given narrow band. If far-field points remain when this
             condition is met a masked arrays are returned. The default
             value is 0.0 which means no narrow band limit.

    periodic : bool or an array-like of len phi.ndim, optional
               specifies whether and in which directions periodic
               boundary conditions are used. True sets periodic
               boundary conditions in all directions. An array-like
               (interpreted as True or False values) specifies the
               absence or presence of periodic boundaries in
               individual directions. The default value is False,
               i.e., no periodic boundaries in any direction.

    Returns
    -------
    (d, f_ext) : tuple
        a tuple containing the signed distance function d and the
        extension velocities f_ext.

    """
    phi, dx, flag, ext_mask, periodic = \
                pre_process_args(phi, dx, narrow, periodic, ext_mask)

    distance, f_ext = cFastMarcher(phi, dx, flag, speed, ext_mask,
                                   int(self_test), EXTENSION_VELOCITY,
                                   order, narrow, periodic)
    distance = post_process_result(distance)
    f_ext = post_process_result(f_ext)

    return distance, f_ext
