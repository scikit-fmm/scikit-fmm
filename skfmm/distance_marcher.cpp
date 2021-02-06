#include "distance_marcher.h"
#include "math.h"

double distanceMarcher::updatePointOrderOne(int i)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  for (int dim=0; dim<dim_; dim++)
  {
    double value = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value))
        {
          value = distance_[naddr];
        }
      }
    }
    if (value<maxDouble)
    {
      a+=idx2_[dim];
      b-=idx2_[dim]*2*value;
      c+=idx2_[dim]*pow(value,2);
    }
  }
  double tmp = solveQuadratic(i,a,b,c);
  return tmp;
}


// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double distanceMarcher::updatePointOrderTwo(int i)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  for (int dim=0; dim<dim_; dim++)
  {
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
  double tmp = solveQuadratic(i,a,b,c);
  return tmp;
}


// find and return the correct root
double distanceMarcher::solveQuadratic(int i, const double &a,
                                       const double &b,
                                       double &c)
{
  c-=1;
  double det = pow(b,2)-4*a*c;
  if (det > 0)
  {
    if (phi_[i] > doubleEpsilon) { return (-b + sqrt(det)) / 2.0 / a; }
    else                         { return (-b - sqrt(det)) / 2.0 / a; }
  }
  else
  {
    return 0;
  }
}


void distanceMarcher::initalizeFrozen()
{
  //loop over phi to find zero values
  //  and mark them as frozen
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] != Mask && phi_[i]==0.0)
    {
      flag_[i]=Frozen;
      distance_[i]=0.0;
    }
  }
  //loop over all of phi and for each point check each direction
  //  to see if we cross the zero level set
  //     if so calculate the minimum distance to the zero level set
  //     mark as frozen.
  for (int i=0; i<size_; i++)
  if (flag_[i] == Far)
  {
    double ldistance[MaximumDimension];
    bool borders=false;
    for (int dim=0; dim<dim_; dim++)
    {
      ldistance[dim]=0;
      for (int j=-1; j<2; j+=2) // each direction
      {
        int naddr = _getN(i,dim,j,Mask);
        if (naddr!=-1 && phi_[i] * phi_[naddr]<0)
        {
          // this cell and neighbor span the zero level set.
          borders=true;
          //calculate the distance to the zero level set.
          double d = dx_[dim]*phi_[i]/(phi_[i]-phi_[naddr]);
          if (ldistance[dim]==0 || ldistance[dim]>d)
          {
            ldistance[dim] = d;
          }
        }
      } // for each direction
    } // for each dimension
    if (borders)
    {
      double dsum = 0;
      for (int dim=0; dim<dim_; dim++)
        if (ldistance[dim]>0) dsum += 1/ldistance[dim]/ldistance[dim];
      if (phi_[i]<0)
        distance_[i] = -sqrt(1/dsum);
      else distance_[i] = sqrt(1/dsum);
      flag_[i]=Frozen;
    }
  }// for each point in the far field
}
