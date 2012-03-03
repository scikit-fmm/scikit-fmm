//fast_marching.h

const unsigned int MaximumDimension  = 12;
const char Far    = 0;
const char Narrow = 1;
const char Frozen = 2;
const char Mask   = 3;

class heap;

using namespace std;

extern "C" {

class baseMarcher
{
public:
  baseMarcher(double *phi,      double *dx,  int *flag,
              double *distance, int ndim,    int *shape,
              bool self_test);

  virtual          ~baseMarcher();
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


  inline int       _getN(int current, int dim, int dir, int flag);
  void             march();

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
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test)
  {
    march();
  }
  virtual ~distanceMarcher() { }

protected:
  virtual void             initalizeFrozen();
  virtual double           updatePoint(int i);
};

class travelTimeMarcher : public baseMarcher
{

  travelTimeMarcher(double *phi,      double *dx, int *flag,
                    double *distance, int ndim,   int *shape,
                    bool self_test,
                    double *speed) :
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test)
  {
    speed_ = speed;
    march();
  }
  virtual ~travelTimeMarcher() { }

protected:
  virtual void             initalizeFrozen() {}
  virtual double           updatePoint(int i) {}

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
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test)
  {
    speed_ = speed;
    f_ext_ = f_ext;
    march();
  }
  virtual ~extensionMarcher() { }

protected:
  virtual void             initalizeFrozen() {}
  virtual double           updatePoint(int i) {}

private:
  double *speed_;
  double *f_ext_; // return value modified in place
};

} // extern "C"
