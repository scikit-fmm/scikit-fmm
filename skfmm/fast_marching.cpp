//fast_marching.cpp

// this file is a C++ implementation of the fast marching method. The
// code in this file is independent of Python.

#include "fast_marching.h"
#include "math.h"
#include <limits>
#include <iostream>

extern "C" {

const double doubleEpsilon    = numeric_limits<double>::epsilon();
const double maxDouble        = numeric_limits<double>::max();

fastMarcher::fastMarcher(
  double *phi,   double *dx,       int *flag,
  double *speed, double *distance,
  int     ndim,  int    *shape,    bool self_test)
{
  error_      =   1;
  phi_        =   phi;
  dx_         =   dx;
  flag_       =   flag;
  if (speed==0)
  {
    hasSpeed_ =   false;
    speed_    =   0;
  }
  else
  {
    hasSpeed_ =   true;
    speed_    =   speed;
  }
  distance_   =   distance;
  dim_        =   ndim;
  size_       =   1;
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
  initalizeFrozen();
  int maxHeap=0;
  for (int i=0; i<size_; i++)
    if (flag_[i] == Far) maxHeap++;
  heap_ = new heap(maxHeap, self_test);
  heapptr_ = new int[size_];

  if (hasSpeed_)
  {
    for (int i=0; i<size_; i++)
    {
      // we need to be carefull here: very small speed values can result
      // in an overflow
      if (speed_[i]<doubleEpsilon) flag_[i]=Mask;
      else
      {
        if (flag_[i]==Frozen || flag_[i]==Narrow)
        {
          // convert distance to time
          if (flag_[i]==Narrow)
            heap_->set(heapptr_[i],fabs(distance_[i]/speed_[i]));
          distance_[i]=fabs(distance_[i]/speed_[i]);
        }
      }
    }
  }
  initalizeNarrow();
  solve();
}

inline void fastMarcher::_getIndex(int current,
                                   int coord[MaximumDimension])
{
  int rem = current;
  for (int i=0; i<dim_; i++)
  {
    coord[i] = rem/shift_[i]; // integer assignment is like floor()
    rem -= coord[i]*shift_[i];
  }
}

inline int fastMarcher::_getN(int current, int dim, int dir, int flag)
{
  // assume c order.
  // for a given point find his neighbor in the given dimension
  // and direction. Return -1 if not possible.
  // consult shape_ information
  int coord[MaximumDimension];
  _getIndex(current, coord);
  int newc = coord[dim]+dir;
  if (newc >= shape_[dim] || newc < 0) return -1;
  int newa = current + dir*shift_[dim];
  if (flag_[newa]==flag)  return -1;
  _getIndex(newa, coord);
  return newa;
}


fastMarcher::~fastMarcher()
{
  delete   heap_;
  delete[] heapptr_;
}


void fastMarcher::initalizeFrozen()
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

void fastMarcher::initalizeNarrow()
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
            double d     =  updatePoint(i);
            distance_[i] =  d;
            heapptr_[i]  =  heap_->push(i,fabs(d));
          }
        } // for each direction
      } // for each dimension
    } // each far field point
  }
}

// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double fastMarcher::updatePoint(int i)
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
          if (naddr2!=-1 && flag_[naddr2]==Frozen && distance_[naddr2]<value1)
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
  if (hasSpeed_)
    c-= 1/pow(speed_[i],2);
  else
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
  if (hasSpeed_) return r0;
  if (phi_[i] > doubleEpsilon) return r0;
  else return r1;
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
void fastMarcher::solve()
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

    for (int dim=0; dim<dim_; dim++)
    {
      for (int j=-1; j<2; j+=2) // each direction
      {
        int naddr = _getN(addr,dim,j,Frozen);
        if (naddr!=-1 && flag_[naddr]!=Frozen)
        {
          if (flag_[naddr]==Narrow)
          {
            double d = updatePoint(naddr);
            if (d)
            {
              heap_->set(heapptr_[naddr],fabs(d));
              distance_[naddr]=d;
            }
          }
          else if (flag_[naddr]==Far)
          {
            double d = updatePoint(naddr);
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

} // extern "C"
