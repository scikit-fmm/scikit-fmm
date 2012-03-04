//fast_marching.h

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


  double          * distance_; // return value modified in place
  double          * phi_;
  double          * dx_;
  int             * flag_;
  int               error_;
  int               dim_;            // number of dimensions
  int               size_;           // flat size
  double            idx2_[MaximumDimension];
};


class distanceMarcher : public baseMarcher
{
public:
  distanceMarcher(double *phi,      double *dx, int *flag,
                  double *distance, int ndim,   int *shape,
                  bool self_test) :
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test) { }
  virtual ~distanceMarcher() { }

protected:
  virtual double           solveQuadratic(int i, const double &a,
                                          const double &b, double &c);

  virtual void             initalizeFrozen();
  virtual double           updatePoint(int i);
};

class travelTimeMarcher : public distanceMarcher
{
public:
  travelTimeMarcher(double *phi,      double *dx, int *flag,
                    double *distance, int ndim,   int *shape,
                    bool self_test,
                    double *speed) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test),
    speed_(speed)
  {
    for (int i=0; i<size_; i++)
    {
      // we need to be carefull here: very small speed values can result
      // in an overflow
      if (speed_[i]<doubleEpsilon) flag_[i]=Mask;
    }
  }

  virtual ~travelTimeMarcher() { }

protected:
  virtual void             initalizeFrozen();
  virtual double           solveQuadratic(int i, const double &a,
                                          const double &b, double &c);
private:
  double *speed_;
};

class extensionMarcher : public baseMarcher
{
public:
  extensionMarcher(double *phi,      double *dx,   int *flag,
                   double *distance, int     ndim, int *shape,
                   bool self_test,
                   double *speed,
                   double *f_ext) :
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test),
    speed_(speed), f_ext_(f_ext)  { }
  virtual ~extensionMarcher() { }

protected:
  virtual void             initalizeFrozen() {}
  virtual double           updatePoint(int i) { return 0.0; }

private:
  double *speed_;
  double *f_ext_; // return value modified in place
};

} // extern "C"
