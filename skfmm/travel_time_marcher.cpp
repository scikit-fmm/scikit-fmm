//travel_time_marcher.cpp

#include "travel_time_marcher.h"
#include "math.h"
#include "heap.h"

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

void travelTimeMarcher::initalizeNarrow()
{
  // for each point in the far field check if neighbor is frozen
  // if so calculate distance and insert into heap
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] == Far)
    {
      for (int dim=0; dim<dim_; dim++)
      {
        for (int j=-1; j<2; j+=2) // each direction
        {
          int naddr = _getN(i,dim,j,Mask);
          if (naddr!=-1 && flag_[naddr]==Frozen)
          if (flag_[i]==Far)
          {
            flag_[i]     =  Narrow;
            double d     =  updatePointOrderOne(i);
            distance_[i] =  d;
            heapptr_[i]  =  heap_->push(i,fabs(d));
          }
        } // for each direction
      } // for each dimension
    } // each far field point
  }
}


double travelTimeMarcher::solveQuadratic(int i, const double &a,
                                         const double &b,
                                         double &c)
{
  c -= 1/pow(speed_[i],2);
  double r0=0;
  double det = pow(b,2)-4*a*c;
  if (det>0)
  {
    r0 = (-b+sqrt(det))/2.0/a;
  }
  else
  {
    return 0;
  }
  return r0;
}


