//sDistanceMarcher.cpp

// this file is a C++ implementation of the fast marching method. The
// code in this file is independent of Python.

#include "sDistanceMarcher.h"
#include "no_malloc_heap.h"
#include "math.h"

extern "C" {


sDistanceMarcher::sDistanceMarcher()
{
  order_      =   0;
  heapptr_    =   0;
  error_      =   0;
  phi_        =   0;
  dx_         =   0;
  flag_       =   0;
  distance_   =   0;
  dim_        =   0;
  size_       =   0;
  self_test_  =   0;
}

sDistanceMarcher::~sDistanceMarcher()
{
}

void sDistanceMarcher::set(double *phi,  double *distance, double *dx,
                           long   *flag, long   *hp,       double *hd,
                           long   *hp,   long   *hi1,      long   *hi2,
                           long   *hi3,  int     ndim,     int    *shape,
                           bool   self_test,               int order);
{
  order_      =   order;
  error_      =   1;
  phi_        =   phi;
  dx_         =   dx;
  flag_       =   flag;
  distance_   =   distance;
  dim_        =   ndim;
  size_       =   1;
  self_test_  =   self_test;
  heapptr_    =   hp;

  for (int i=0; i<dim_; i++)
  {
    shape_[i]  = shape[i];
    size_     *= shape[i];
    idx2_[i]   = 1/dx_[i]/dx_[i];
  }
  for (int i=0; i<dim_; i++)
  {
    int prod=1;
    for (int j=i+1; j<dim_; j++) prod*=shape_[j];
    shift_[i]=prod;
  }
  heap_.set(size_, hd, hi1, hi2, hi3, self_test);
}

void sDistanceMarcher::march()
{
  initalizeFrozen();
  initalizeNarrow();
  solve();
  cleanUp();
}

void sDistanceMarcher::initalizeNarrow()
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
            double d;
            if (order_ == 2)
              d =  updatePointOrderTwo(i);
            else
              d =  updatePointOrderOne(i);

            distance_[i] =  d;
            heapptr_[i]  =  heap_->push(i,fabs(d));
          }
        } // for each direction
      } // for each dimension
    } // each far field point
  }
}

// now we apply the fast marching algorithm main loop
// (1) take the smallest narrow band element and
//     freeze it.
// (2) for each neighbor of the frozen point calculate distance based
// on frozen elements
//     - mark each neighbor as narrow band and stick it
//       into the heap
//     - if the neighbor is already in the heap update the
//       distance value.
void sDistanceMarcher::solve()
{
  int frozenCount=0;
  for (int i=0; i<size_; i++)
    if (flag_[i] == Frozen) frozenCount++;
  if (!frozenCount)
  {
    error_ = 2;
    return;
  }
  int i=0;
  while (! heap_->empty())
  {
    i++;
    double  value   = 0;
    int     addr    = 0;
    heap_->pop(&addr, &value);
    flag_[addr]=Frozen;
    finalizePoint(addr, value);

    for (int dim=0; dim<dim_; dim++)
    {
      for (int j=-1; j<2; j+=2) // each direction
      {
        int naddr = _getN(addr,dim,j,Frozen);
        if (naddr!=-1 && flag_[naddr]!=Frozen)
        {
          if (flag_[naddr]==Narrow)
          {

            double d;
            if (order_ == 2)
              d =  updatePointOrderTwo(naddr);
            else
              d =  updatePointOrderOne(naddr);
            if (d)
            {
              heap_->set(heapptr_[naddr],fabs(d));
              distance_[naddr]=d;
            }
          }
          else if (flag_[naddr]==Far)
          {
            double d;
            if (order_ == 2)
              d =  updatePointOrderTwo(naddr);
            else
              d =  updatePointOrderOne(naddr);
            if (d)
            {
              distance_[naddr]=d;
              flag_[naddr]=Narrow;
              heapptr_[naddr] = heap_->push(naddr,fabs(d));
            }
          }
        }
      }
    }
  }
  // add back mask here. The python wrapper will look for elements
  // equal to mexDouble and add the mask back
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] == Mask) distance_[i] = maxDouble;
    if (flag_[i] == Far)  distance_[i] = maxDouble;
  }
  error_ = 0;
  return;
}


double sDistanceMarcher::updatePointOrderOne(int i)
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
        if (distance_[naddr]<value)
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
  return solveQuadratic(i,a,b,c);
}


// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double sDistanceMarcher::updatePointOrderTwo(int i)
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
        if (distance_[naddr]<value1)
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
  return solveQuadratic(i,a,b,c);
}


// find and return the correct root
double sDistanceMarcher::solveQuadratic(int i, const double &a,
                                       const double &b,
                                       double &c)
{
  c-=1;
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
  if (phi_[i] > doubleEpsilon) return r0;
  else return r1;
}

void sDistanceMarcher::initalizeFrozen()
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

} // extern "C"
