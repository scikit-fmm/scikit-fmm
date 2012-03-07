//travel_time_marcher.cpp

#include "travel_time_marcher.h"
#include "math.h"

void travelTimeMarcher::initalizeFrozen()
{
  distanceMarcher::initalizeFrozen();
  for (int i=0; i<size_; i++)
  {
    if (flag_[i]==Frozen)
    {
      // convert distance to time
      distance_[i]=fabs(distance_[i]/speed_[i]);
    }
  }
}

double travelTimeMarcher::solveQuadratic(int i, const double &a,
                                         const double &b,
                                         double &c)
{
  c -= 1/pow(speed_[i],2);
  double r0=0;
  double r1=0;
  double det = pow(b,2)-4*a*c;
  if (det>0)
  {
    r0 = (-b+sqrt(det))/2.0/a;
    r1 = (-b-sqrt(det))/2.0/a;
  }
  else
  {
    return 0;
  }
  return r0;
}


