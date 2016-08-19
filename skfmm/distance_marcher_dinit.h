//distance_marcher_dinit.h
#pragma once

#include "distance_marcher.h"

class distanceMarcherDInit : public distanceMarcher
{
public:
  distanceMarcherDInit(double *phi,      double *dx, long *flag,
                       double *distance, int ndim,   int *shape,
                       bool self_test,   int order,  double narrow,
                       int periodic, double *dinit) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test,
                    order, narrow, periodic)
  {
    dinit_ = dinit;
  }
  virtual ~distanceMarcherDInit() { }
  virtual void initalizeFrozen();
private:
  double * dinit_;
};
