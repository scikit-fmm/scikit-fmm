//base_marcher.h
#pragma once
const unsigned int MaximumDimension  = 12;
const char Far    = 0;
const char Narrow = 1;
const char Frozen = 2;
const char Mask   = 3;

using namespace std;
#include <limits>
const double doubleEpsilon    = numeric_limits<double>::epsilon();
const double maxDouble        = numeric_limits<double>::max();

class heap;

extern "C" {

class baseMarcher
{
public:
  baseMarcher(double *phi,      double *dx,  long *flag,
              double *distance, int ndim,    int *shape,
              bool self_test,   int order);

  virtual          ~baseMarcher();
  void             march();
  int              getError() const { return error_;}

private:
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
  int             * heapptr_;        // heap back pointers
  heap            * heap_;
  int               shape_[MaximumDimension];    // size of each dimension
  int               shift_[MaximumDimension];
  bool              self_test_;

protected:
  // derived classes must implement these functions
  virtual void     initalizeFrozen() = 0;
  virtual double   updatePointOrderTwo(int i) = 0;
  virtual double   updatePointOrderOne(int i) = 0;

  // derived classes may implement these functions
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
