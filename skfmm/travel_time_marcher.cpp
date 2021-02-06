//travel_time_marcher.cpp

#include "travel_time_marcher.h"
#include "math.h"
#include "heap.h"
#include <stdexcept>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element

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

double travelTimeMarcher::updatePointOrderTwo(int i)
{
    double res = updatePointOrderTwo(i, std::set<int>());
    if(res == std::numeric_limits<double>::infinity())
        throw std::runtime_error("Unreachable voxel");
    else
        return res;
}


// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double travelTimeMarcher::updatePointOrderTwo(int i, std::set<int>avoid_dim)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  for (int dim=0; dim<dim_; dim++)
  {
    if(avoid_dim.find(dim) != avoid_dim.end()) continue; //we should avoid this dimension
    double value1 = maxDouble;
    double value2 = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value1))
        {
          value1 = distance_[naddr];
          int naddr2 = _getN(i,dim,j*2,Mask);
          if (naddr2!=-1 &&
              flag_[naddr2]==Frozen &&
              ((distance_[naddr2]<=value1 && value1 >=0) ||
               (distance_[naddr2]>=value1 && value1 <=0)))
          {
            value2=distance_[naddr2];
            if (phi_[naddr2] * phi_[naddr] < 0  || phi_[naddr2] * phi_[i] < 0)
              value2 *= -1;
          }
        }
      }
    }
    if (value2<maxDouble)
    {
      double tp = oneThird*(4*value1-value2);
      a+=idx2_[dim]*aa;
      b-=idx2_[dim]*2*aa*tp;
      c+=idx2_[dim]*aa*pow(tp,2);
    }
    else if (value1<maxDouble)
    {
      a+=idx2_[dim];
      b-=idx2_[dim]*2*value1;
      c+=idx2_[dim]*pow(value1,2);
    }
  }
  try {
    double res = solveQuadratic(i,a,b,c);
    return res;
  } catch(std::runtime_error& err) {
    //if the determinant is negative, we try to reach the voxel with one dimension less and take the minimum
    if(avoid_dim.size() == dim_) return std::numeric_limits<double>::infinity(); //end of the recursion, use inf so that it is discarded selecting the minimum
    std::vector<double> sols;
    for (int ind=0; ind<dim_; ind++){
        //remove one dimension more than what we are already doing
        std::set<int> tempset = avoid_dim;
        std::pair<std::set<int>::iterator, bool> ret = tempset.insert(ind);
        if(!ret.second) continue; //avoid recursive call on identical parameters (the set already had *ind* in it)
        sols.push_back(updatePointOrderTwo(i,tempset));
    }
    if(sols.size()==0) return std::numeric_limits<double>::infinity();//All the derivates with different dimensionalities are 0
    return *std::min_element(sols.begin(), sols.end());
  }
}


double travelTimeMarcher::solveQuadratic(int i, const double &a,
                                         const double &b,
                                         double &c)
{
  c -= 1/pow(speed_[i],2);
  double r0=0;
  double det = pow(b,2)-4*a*c;
  if (det>=0)
  {
    r0 = (-b+sqrt(det))/2.0/a;
  }
  else
  {
    throw std::runtime_error("Negative discriminant in time marcher quadratic.");
  }
  return r0;
}
