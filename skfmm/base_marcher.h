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
  baseMarcher(double *phi,      double *dx,  int *flag,
              double *distance, int ndim,    int *shape,
              bool self_test);

  virtual          ~baseMarcher();
  void             march();
  int              getError() const { return error_;}

private:
  void              initalizeNarrow();
  void              solve();
  inline void       _getIndex(int current, int coord[MaximumDimension]);

  int             * heapptr_;        // heap back pointers
  heap            * heap_;
  int               shape_[MaximumDimension];    // size of each dimension
  int               shift_[MaximumDimension];
  bool              self_test_;

protected:
  // derived classes must implement these functions
  virtual void     initalizeFrozen() = 0;
  virtual double   updatePoint(int i) = 0;
  int              _getN(int current, int dim, int dir, int flag);
  virtual void     cleanUp() {}

  double          * distance_; // return value modified in place
  double          * phi_;
  double          * dx_;
  int             * flag_;
  int               error_;
  int               dim_;            // number of dimensions
  int               size_;           // flat size
  double            idx2_[MaximumDimension];
};

} // extern "C"
