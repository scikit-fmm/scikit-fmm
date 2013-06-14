from pheap import pheap as _pheap

# This wrapper exists only to add docstrings to the heap class.
# Cython cannot do this automatically at the time of writing.
# see: http://docs.cython.org/src/userguide/special_methods.html

class heap(_pheap):
    """This class implements a binary min heap data structure to support the
    fast marching algorithm. A min heap is a list which has the property
    that the smallest element is always the first element.

    The fast marching method uses this data structure to track elements in
    the solution narrow band. The fast marching method needs to know which
    element in the narrow band is nearest the zero level-set at each
    iteration.

    When a new point enters the solution narrow band the element is added
    to the heap.

    New elements are added to the heap with the push() method. The
    push() method returns an integer (a heap index) which is used
    later to refer to the element.

    As the solution evolves the distance of points already in the heap can
    be updated via the set() method. The set() method takes the heap index
    returned by the push() along with a new distance value.

    The narrow band element nearest the zero-level set is taken off
    the heap with the pop() method. The address of the top element is
    returned to the caller.

    The constructor for heap needs to know the number of elements that
    will enter the narrow band.

    >>> from skfmm import heap
    >>> h = heap(10,False)
    >>> h.push(0,0.2)
    >>> h.push(1,0.3)
    >>> h.push(2,0.1)
    >>> h.set(1, 0.01)

    >>> addr, dist = h.pop()
    >>> assert dist==0.01
    >>> addr, dist = h.pop()
    >>> assert dist==0.1
    >>> addr, dist = h.pop()
    >>> assert dist==0.2
    >>> assert h.empty()

    """
    def __init__(self, _, __):
        """Create a new min heap.

        Parameters
        ----------
        n : int
            The maximum size of the heap.

        self_test : Boolean

                    If True a consistency check is made after each
                    heap operation.This is used in testing and results
                    in a slower calculation.
        """

    def push(self, addr, value):
        """ Add a value to the heap, give an address and a value  """
        return self._push(addr, value)

    def pop(self):
        """Remove and return the top element on the heap. The return value is
        a tuple of address and value

        """
        return self._pop()

    def set(self, addr, value):
        """Update the value of a point already in the heap, the heap
        ordering is updated if needed. Input is the address returned by push

        """
        self._set(addr, value)

    def empty(self):
        """Return True if the heap is empty

        """
        return self._empty()

    def peek(self):
        """Return the value on top of the heap without removing it.

        """
        return self._peek()
