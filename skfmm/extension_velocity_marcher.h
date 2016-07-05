//extension_velocity_marcher.h

#include "distance_marcher.h"

class extensionVelocityMarcher : public distanceMarcher
{
public:
  extensionVelocityMarcher(double *phi,      double *dx,    long *flag,
                           double *distance, int     ndim,  int *shape,
                           bool self_test,   int order,     long *ext_mask,
                           double *speed,
                           double *f_ext,    double narrow, int periodic) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test,
                    order, narrow, periodic),
  speed_(speed), f_ext_(f_ext), ext_mask_(ext_mask) { }
  virtual ~extensionVelocityMarcher() { }

protected:
  virtual void             initalizeFrozen();
  virtual void             finalizePoint(int i, double phi_i);
  virtual void             cleanUp();

private:
  double *speed_;
  double *f_ext_;
  long *ext_mask_;
};
