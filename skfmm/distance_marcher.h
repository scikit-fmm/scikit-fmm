//distance_marcher.h
#pragma once
#include "base_marcher.h"

class distanceMarcher : public baseMarcher
{
public:
  distanceMarcher(double *phi,      double *dx, long *flag,
                  double *distance, int ndim,   int *shape,
                  bool self_test,   int order,  double narrow,
                  int periodic):
    baseMarcher(phi, dx, flag, distance, ndim, shape, self_test, order,
                narrow, periodic) { }
  virtual ~distanceMarcher() { }

protected:
  virtual double           solveQuadratic(int i, const double &a,
                                          const double &b, double &c);

  virtual void             initalizeFrozen();
  virtual double           updatePointOrderOne(int i);
  virtual double           updatePointOrderTwo(int i);
};
