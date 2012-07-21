//sDistanceMarcher.h
#pragma once

#include "no_malloc_heap.h"
#include "base_marcher.h"

extern "C" {

class sDistanceMarcher
{
public:
  sDistanceMarcher();

  virtual          ~sDistanceMarcher();
  void             set(double *phi,  double *distance, double *dx,
                       long   *flag, long   *hp,       double *hd,
                       long   *hi1,      long   *hi2,
                       long   *hi3,  int     ndim,     int    *shape,
                       bool   self_test,               int order);
  void             march();
  int              getError() const { return error_;}

private:
  virtual double           solveQuadratic(int i, const double &a,
                                          const double &b, double &c);

  void             initalizeNarrow();
  void             solve();
  void             _getIndex(int current, int coord[MaximumDimension])
  {
    int rem = current;
    for (int i=0; i<dim_; i++)
    {
      coord[i] = rem/shift_[i]; // integer assignment is like floor()
      rem -= coord[i]*shift_[i];
    }
  }

  int               order_;
  long            * heapptr_;        // heap back pointers
  nm_heap           heap_;
  int               shape_[MaximumDimension];    // size of each dimension
  int               shift_[MaximumDimension];
  bool              self_test_;

private:
  virtual void     initalizeFrozen();
  virtual double   updatePointOrderTwo(int i);
  virtual double   updatePointOrderOne(int i);

  virtual void     cleanUp() { }
  virtual void     finalizePoint(int i, double phi_i) { }

  int              _getN(int current, int dim, int dir, int flag)
  {
    // assume c order.
    // for a given point find his neighbor in the given dimension
    // and direction. Return -1 if not possible.
    // consult shape_ information
    int coord[MaximumDimension];
    _getIndex(current, coord);
    int newc = coord[dim]+dir;
    if (newc >= shape_[dim] || newc < 0) return -1;
    int newa = current + dir*shift_[dim];
    if (flag_[newa]==flag)  return -1;
    _getIndex(newa, coord);
    return newa;
  }

  double          * distance_; // return value modified in place
  double          * phi_;
  double          * dx_;
  long            * flag_;
  int               error_;
  int               dim_;            // number of dimensions
  int               size_;           // flat size
  double            idx2_[MaximumDimension];
};

} // extern "C"
