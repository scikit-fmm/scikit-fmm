//extension_velocity_marcher.h

#include "distance_marcher.h"

class extensionVelocityMarcher : public distanceMarcher
{
public:
  extensionVelocityMarcher(double *phi,      double *dx,   int *flag,
                   double *distance, int     ndim, int *shape,
                   bool self_test,
                   double *speed,
                   double *f_ext) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test),
    speed_(speed), f_ext_(f_ext)  { }
  virtual ~extensionVelocityMarcher() { }

protected:
  virtual void             initalizeFrozen();
  virtual double           updatePoint(int i);
  virtual void             cleanUp();

private:
  double *speed_;
  double *f_ext_;
};
