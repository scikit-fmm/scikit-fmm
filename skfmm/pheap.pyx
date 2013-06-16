# distutils: language = c++
# distutils: sources = heap.cpp

from libcpp cimport bool

cdef extern from "heap.h":
   cdef cppclass heap:
      heap(int, bool) except +
      int push(int, double) except +
      void pop(int *, double *) except +
      void set(int, double) except +
      bool empty()
      double peek() except +

cdef class pheap:
    cdef heap *thisptr
    def __cinit__(self, n, self_test):
        self.thisptr = new heap(n, self_test)

    def __dealloc__(self):
        del self.thisptr

    def _push(self, addr, value):
        return self.thisptr.push(addr, value)

    def _pop(self):
        cdef int addr
        cdef double value
        self.thisptr.pop(&addr, &value)
        return addr, value

    def _set(self, addr, value):
        self.thisptr.set(addr, value)

    def _empty(self):
        return self.thisptr.empty()

    def _peek(self):
        cdef double dist
        return self.thisptr.peek()
