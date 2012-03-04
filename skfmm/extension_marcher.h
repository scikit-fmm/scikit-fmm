//extension_marcher.h

#include "base_marcher.h"

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
