//fast_marching.cpp

// this file is a C++ implementation of the fast marching method. The
// code in this file is independent of Python.

#include "base_marcher.h"
#include "heap.h"
#include "math.h"
#include <vector>

extern "C" {


baseMarcher::baseMarcher(
  double *phi,      double *dx,   long *flag,
  double *distance, int     ndim, int *shape,
  bool self_test,   int order,    double narrow,
  int periodic)
{
  narrow_     =   narrow;
  order_      =   order;
  error_      =   1;
  phi_        =   phi;
  dx_         =   dx;
  flag_       =   flag;
  distance_   =   distance;
  dim_        =   ndim;
  size_       =   1;
  self_test_  =   self_test;
  heapptr_    =   0;
  heap_       =   0;
  periodic_   =   periodic;

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
}

void baseMarcher::march()
{
  initalizeFrozen();

  int maxHeap=0;
  for (int i=0; i<size_; i++)
    if (flag_[i] == Far) maxHeap++;
  heap_ = new heap(maxHeap, self_test_);
  heapptr_ = new int[size_];
  initalizeNarrow();
  solve();
  cleanUp();
}

baseMarcher::~baseMarcher()
{
  delete   heap_;
  delete[] heapptr_;
}

void baseMarcher::initalizeNarrow()
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
void baseMarcher::solve()
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
  std::vector<int> toFreeze;

  while (! heap_->empty())
  {
    toFreeze.clear();
    i++;
    double  value   = 0;
    int     addr    = 0;
    heap_->pop(&addr, &value);

    if ((narrow_ != 0) && (fabs(value) > narrow_))
    {
      break;
    }
    flag_[addr]=Frozen;
    finalizePoint(addr, value);
    toFreeze.push_back(addr);

    bool done=false;
    while (!done)
    {
      if (! heap_->empty() && value == heap_->peek())
      {
        double  l_value   = 0;
        int     l_addr    = 0;
        heap_->pop(&l_addr, &l_value);
        flag_[l_addr]=Frozen;
        finalizePoint(l_addr, l_value);
        toFreeze.push_back(l_addr);
      }
      else
        done = true;
    }

    for (unsigned int k=0; k<toFreeze.size(); k++)
    {
      int addr = toFreeze[k];

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
          //==========================================================
          // update the far point in the second order stencil
          // "jump" over a Frozen point if needed
          if (order_ == 2)
          {
            int local_naddr = _getN(addr,dim,j,Mask);
            if (local_naddr!=-1 && flag_[local_naddr]==Frozen)
            {
              int naddr2 = _getN(addr,dim,j*2,Frozen);
              if (naddr2!=-1 && flag_[naddr2]==Narrow)
              {
                double d = updatePointOrderTwo(naddr2);
                if (d)
                {
                  heap_->set(heapptr_[naddr2], fabs(d));
                  distance_[naddr2]=d;
                }
              }
            }
          }
          //==========================================================
        } // for each direction
      } // for each dimension
    }
  } // main loop of Fast Marching Method

  // add back mask here. The python wrapper will look for elements
  // equal to maxDouble and add the mask back
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] != Frozen) distance_[i] = maxDouble;
  }
  error_ = 0;
  return;
}

} // extern "C"
