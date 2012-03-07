//distance_marcher.h
#pragma once
#include "base_marcher.h"

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
