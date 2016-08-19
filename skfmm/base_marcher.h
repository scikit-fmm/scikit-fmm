//base_marcher.h
#pragma once
const unsigned int MaximumDimension  = 12;
const char Far    = 0;
const char Narrow = 1;
const char Frozen = 2;
const char Mask   = 3;

#include <limits>
using namespace std;
#define doubleEpsilon            numeric_limits<double>::epsilon()
#define maxDouble                numeric_limits<double>::max()

class heap;

extern "C" {

class baseMarcher
{
public:
  baseMarcher(double *phi,      double *dx,  long *flag,
              double *distance, int ndim,    int *shape,
              bool self_test,   int order,   double narrow,
              int periodic);

  virtual          ~baseMarcher();
  void             march();
  int              getError() const { return error_;}

private:
  void             initalizeNarrow();
  void             solve();
  void             _setIndex(int current, int coord[MaximumDimension])
  {
    int rem = current;
    for (int i=0; i<dim_; i++)
    {
      coord[i] = rem/shift_[i]; // integer assignment is like floor()
      rem -= coord[i]*shift_[i];
    }
  }

  int              _getIndex(int coord[MaximumDimension])
  {
    int ret = 0;
    for (int i=0; i<dim_; i++) ret += coord[i]*shift_[i];
    return ret;
  }

  double            narrow_;
  int               order_;
  int             * heapptr_;        // heap back pointers
  heap            * heap_;
  int               shape_[MaximumDimension];    // size of each dimension
  int               shift_[MaximumDimension];
  int               periodic_;
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
    // for a given point find the neighbor in the given dimension and
    // direction.
    // Return -1 if neighbor point is invalid (out of bounds)
    // Return -1 if neighbor point has flag_ value equal to the flag input
    int coord[MaximumDimension];
    _setIndex(current, coord);
    int newc = coord[dim]+dir;

    if (periodic_ & 1<<dim) {
      if (newc==-1)                  newc = shape_[dim]-1;
      else if (newc==-2)             newc = shape_[dim]-2;
      else if (newc==shape_[dim])    newc = 0;
      else if (newc==shape_[dim]+1)  newc = 1;
      coord[dim] = newc;
      return _getIndex(coord);
    }

    if (newc >= shape_[dim] || newc < 0) return -1;
    int newa = current + dir*shift_[dim];
    if (flag_[newa]==flag)  return -1;
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
