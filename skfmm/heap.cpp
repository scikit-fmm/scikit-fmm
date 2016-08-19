//heap.cpp

//////////////////////////////////////////////////////////////////////////////
// The heap index (return by the push() method) is a map from the grid      //
// address to the index of heap's internal data representation.             //
//                                                                          //
// c_ is the current size of the heap.                                      //
// depth_ is the final size of the heap.                                    //
//                                                                          //
// distance_[] is the unsigned distance from the zero level set to each     //
// element in the heap.                                                     //
//                                                                          //
// address_[] is the (original) grid address of each element in the heap.   //
//                                                                          //
// backPointer_[] is a map from the index of distance_ (or address_) to the //
// current location of the element in the heap.                             //
//                                                                          //
// heap_ is an array of integer indices into the distance_ (or address_)    //
// array the heap invariant is maintained by moving elements in this list.  //
//                                                                          //
// if selfTest_ is true a consistency check is done after each operation.   //
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// Currently, the heap constructor needs to know the number of elements //
// that will enter the heap during the calculation.                     //
//////////////////////////////////////////////////////////////////////////

#include "heap.h"
#include <stdio.h>
#include <limits>
#include <stdexcept>

heap::heap(int maxLength, bool selfTest)
{
  maxLength_    = maxLength;
  heapLength_   = 0;
  listLength_   = 0;
  selfTest_     = selfTest;
  distance_     = new double [maxLength_];
  backPointer_  = new int    [maxLength_];
  heap_         = new int    [maxLength_];
  address_      = new int    [maxLength_];
}

heap::~heap()
{
  delete[] distance_;
  delete[] backPointer_;
  delete[] heap_;
  delete[] address_;
}

int heap::push(int address, double value)
{
  if (heapLength_==maxLength_)
    throw std::runtime_error("heap push error: heap full\n");
  heap_[heapLength_]            = listLength_;
  address_[listLength_]         = address;
  distance_[listLength_]        = value;
  backPointer_[listLength_]     = heapLength_;
  listLength_++;
  heapLength_++;
  _siftDown(0,heapLength_-1);
  if (selfTest_) test();
  return listLength_-1;
}

const double &heap::peek() const
{
  if (heapLength_==0)
    throw std::runtime_error("heap peek error: empty heap\n");
  int loc=heap_[0];
  return distance_[loc];
}

void heap::pop(int *address, double *value)
{
  // out: address to element with min value and value
  if (heapLength_==0)
    throw std::runtime_error("heap pop error: empty heap\n");
  int loc                 = heap_[0];
  *value                  = distance_[loc];
  *address                = address_[heap_[0]];
  heap_[0]                = heap_[heapLength_-1];
  backPointer_[heap_[0]]  = 0;
  heapLength_--;
  _siftUp(0);
  if (selfTest_) test();
}

void heap::_siftDown(int startPos, int pos)
{
  int newItem;
  int parent;
  int parentPos;
  newItem = heap_[pos];
  while (pos > startPos)
  {
    parentPos = (pos-1)>>1;
    parent = heap_[parentPos];
    if (distance_[newItem] < distance_[parent])
    {
      heap_[pos]=parent;
      backPointer_[parent]=pos;
      pos=parentPos;
      continue;
    }
    break;
  }
  heap_[pos]=newItem;
  backPointer_[newItem]=pos;
}

void heap::_siftUp(int pos)
{
  int endPos = heapLength_;
  int startPos = pos;
  int newItem;
  int rightPos;
  newItem = heap_[pos];
  int childPos = 2*pos + 1;
  while (childPos < endPos)
  {
    rightPos = childPos + 1;
    if ((rightPos < endPos) && !
        (distance_[heap_[childPos]] < distance_[heap_[rightPos]]))
    {
      childPos = rightPos;
    }
    heap_[pos]=heap_[childPos];
    backPointer_[heap_[childPos]]=pos;
    pos = childPos;
    childPos = 2*pos + 1;
  }
  heap_[pos] = newItem;
  _siftDown(startPos, pos);
}


void heap::set(int index, double newDistance)
{
  double oldDistance = distance_[index];
  int pos = backPointer_[index];
  distance_[index]=newDistance;
  if (newDistance > oldDistance)
  {
    _siftUp(pos);
  }
  if (distance_[heap_[pos]] != newDistance)
  {
    if (selfTest_) test();
    return;
  }
  _siftDown(0,pos);
  if (selfTest_) test();
}

bool heap::empty() const
{
  if (heapLength_==0) return true;
  return false;
}

void heap::test() const
{
  for (int i=0; i<heapLength_; i++)
  {
    int c[2];
    c[0]=2*i+1;
    c[1]=c[0]+1;
    for (int j=0; j<2; j++)
    {
      if (c[j] < heapLength_-1)
      {
        double dp = distance_[heap_[i]];
        double dc = distance_[heap_[c[j]]];
        if (! (dp<=dc))
        {
          throw std::runtime_error("heap invariant violation");
        }
      }
    }
  }
  for (int i=0; i<heapLength_; i++)
  {
    if (! (backPointer_[heap_[i]]==i))
    {
      printf("error %i\n",i);
      throw std::runtime_error("heap backpointer inconsistancy");
    }
  }
}

// void heap::print() const
// {
//   for (int i=0; i<heapLength_; i++)
//   {
//     printf("%i: %i ", i, heap_[i]);
//     printf("(%i)",     address_[heap_[i]]);
//     printf(" %lf", distance_[heap_[i]]);
//     printf(" [%i]\n", backPointer_[heap_[i]]);
//   }
//   printf("\n");
// }
