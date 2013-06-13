# distutils: language = c++
# distutils: sources = heap.cpp

from libcpp cimport bool

cdef extern from "heap.h":
   cdef cppclass heap:
      heap(int, bool) except +
      int push(int, double)
      void pop(int *, double *)
      void set(int, double)
      bool empty()
      double peek()

cdef class pheap:
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
    cdef heap *thisptr
    def __cinit__(self, n, self_test):
        """Create a new min heap.

        Parameters
        ----------
        n : int
            The maximum size of the heap

        self_test : Boolean

                    If True a consistency check is made after each
                    heap operation.This is used in testing and results
                    in a slower calculation.

        """
        self.thisptr = new heap(n,False)
    def __dealloc__(self):
        del self.thisptr
    def push(self, addr, value):
        return self.thisptr.push(addr, value)
    def pop(self):
        cdef int addr
        cdef double dist
        self.thisptr.pop(&addr, &dist)
        return addr, dist
    def set(self, n, dist):
        self.thisptr.set(n, dist)
    def empty(self):
        """
        Return True if the heap is empty
        """
        return self.thisptr.empty()
    def peek(self):
        """
        Return the value on top of the heap without removing it.
        """
        cdef double dist
        dist = self.thisptr.peek()
