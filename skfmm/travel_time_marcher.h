//travel_time_marcher.h
#include "distance_marcher.h"
#include <set>

class heap;

class travelTimeMarcher : public distanceMarcher
{
public:
  travelTimeMarcher(double *phi,      double *dx, long *flag,
                    double *distance, int ndim,   int *shape,
                    bool self_test,   int order,
                    double *speed,    double narrow,
                    int periodic) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test,
                    order, narrow, periodic),
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
  virtual double           updatePointOrderTwo(int i);
  virtual double           updatePointOrderTwo(int i, std::set<int> avoid_dim);
  virtual double           solveQuadratic(int i, const double &a,
                                          const double &b, double &c);
private:
  double *speed_;
};
