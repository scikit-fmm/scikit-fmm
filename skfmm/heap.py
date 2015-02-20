from .pheap import pheap as _pheap

# This wrapper exists only to add docstrings to the heap class.
# Cython cannot do this automatically at the time of writing.
# see: http://docs.cython.org/src/userguide/special_methods.html


class heap(object):
    """.. note::

         Using this class is not required to use
         :py:func:`distance`, :py:func:`travel_time` or
         :py:func:`extension_velocities`. It is provided for
         experimenting with fast marching algorithms.


    This class implements a binary min heap (or heap queue) to support
    the fast marching algorithm. A min heap is a list which has the
    property that the smallest element is always the first element.
    http://en.wikipedia.org/wiki/Binary_heap

    This class differs from the heap queue (:py:obj:`heapq`) in the
    Python standard library because it supports changing the value of
    elements in the heap and maintains forward and backward ids.

    The fast marching method uses this data structure to track
    elements in the solution narrow band. The fast marching method
    needs to know which element in the narrow band is nearest the zero
    level-set at each iteration. When a new point enters the solution
    narrow band the element is added to the heap.

    New elements (an address and an initial value) are added to the
    heap with the :py:meth:`push` method. The :py:meth:`push` method
    returns an integer (a :py:obj:`heap_id`) which is used to update
    the value of the element.

    As the solution evolves, the distance of points already in the heap
    are updated via the :py:meth:`update` method. The
    :py:meth:`update` method takes the :py:obj:`heap_id` returned by the
    :py:meth:`push` along with a new distance value.

    The smallest element is taken off the heap with the :py:meth:`pop`
    method. The address and value of the top element is returned.

    The constructor for heap needs to know the number of elements that
    will enter the heap.

    >>> from skfmm import heap
    >>> h = heap(10)
    >>> h.push(10,0.2)
    0
    >>> h.push(11,0.3)
    1
    >>> h.push(12,0.1)
    2
    >>> h.peek()
    0.1
    >>> h.update(1, 0.01)
    >>> h.peek()
    0.01
    >>> h.pop()
    (11, 0.01)
    >>> h.pop()
    (12, 0.1)
    >>> h.pop()
    (10, 0.2)
    >>> assert h.empty()

    """

    def __init__(self, max_size, self_test=False):
        """Create a new min heap.

        Parameters
        ----------
        max_heap : int

                   The maximum size of the heap.

        self_test : Boolean

                    If True a consistency check is made after each
                    heap operation. This is used in testing and
                    results in a slower calculation.

        """
        self._heap = _pheap(max_size, self_test)

    def push(self, addr, value):
        """Add a value to the heap, give an address and a value.

        Parameters
        ----------
        addr : int

               An id number which is returned when the element is
               popped off the heap.

        value : float

                An initial numerical value for this element.

        Returns
        -------
        heap_id : int

                          An id number into the heap used to update
                          the value of this element. The value of a
                          point in the heap can be updated by calling
                          :py:meth:`update` with this id and a new value.

        """
        return self._heap._push(addr, value)

    def pop(self):
        """Remove and return the address and value of the top element on the
        heap.

        Returns
        -------

        (addr, value) : tuple

                        A tuple of the address given when the element
                        was pushed onto the heap and the current value
                        of the element.

        """
        return self._heap._pop()

    def update(self, heap_id, value):
        """Update the value of a point already in the heap, the heap ordering
        is updated.

        Parameters
        ----------

        heap_id : int

                  The :py:obj:`heap_id` value returned by
                  :py:meth:`push` when this element was added to the
                  heap.

        value : float

                The new value for this element.

        """
        self._heap._set(heap_id, value)

    def empty(self):
        """
        Returns
        -------

        empty : Boolean

                True if the heap is empty.

        """
        return self._heap._empty()

    def peek(self):
        """
        Returns the top (smallest) value on the heap without removing it.

        Returns
        -------
        value : float

                The top (smallest) value on the heap.
        """
        return self._heap._peek()
