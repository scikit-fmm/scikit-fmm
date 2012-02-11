//fast_marching.h
#include "heap.h"

const unsigned int MaximumDimension  = 12;

const char Far    = 0;
const char Narrow = 1;
const char Frozen = 2;
const char Mask   = 3;

using namespace std;

extern "C" {

class fastMarcher
{
 public:
  fastMarcher(double *phi,   double *dx,       int *flag,
	      double *speed, double *distance,
	      int ndim,      int    *shape,    bool self_test);
  virtual ~fastMarcher();

  int              getError() const { return error_;}

 private:
  void             initalizeFrozen();
  void             initalizeNarrow();
  void             solve();
  double           updatePoint(int i);

  inline int       _getN(int current, int dim, int dir, int flag);
  inline void      _getIndex(int current, int coord[MaximumDimension]);

  double          * distance_; // return value modified in place
  double          * phi_;
  double          * dx_;
  int             * flag_;
  double          * speed_;
  int             * heapptr_;        // heap back pointers
  heap            * heap_;
  int               dim_;            // number of dimensions
  int               shape_[MaximumDimension];    // size of each dimension
  int               shift_[MaximumDimension];
  double            idx2_[MaximumDimension];
  int               size_;           // flat size
  bool              hasSpeed_;
  int               error_;
};

} // extern "C"
