//no_malloc_heap.cpp

#include "no_malloc_heap.h"
#include <stdio.h>
#include <limits>
#include <stdexcept>

nm_heap::nm_heap()
{
  selfTest_     = false;
  maxLength_    = 0;
  heapLength_   = 0;
  listLength_   = 0;
  distance_     = 0;
  backPointer_  = 0;
  heap_         = 0;
  address_      = 0;
}

nm_heap::~nm_heap()
{
}

void nm_heap::set_data(int maxLength, double *hd, long *hi1, long *hi2,
                       long *hi3,  bool self_test)
{
  maxLength_    = maxLength;
  heapLength_   = 0;
  listLength_   = 0;
  distance_     = hd;
  backPointer_  = hi1;
  heap_         = hi2;
  address_      = hi3;
  selfTest_     = self_test;
}

int nm_heap::push(int address, double value)
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

void nm_heap::pop(int *address, double *value)
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

void nm_heap::_siftDown(int startPos, int pos)
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

void nm_heap::_siftUp(int pos)
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

void nm_heap::set(int index, double newDistance)
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

bool nm_heap::empty() const
{
  if (heapLength_==0) return true;
  return false;
}

void nm_heap::test() const
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
    if (! backPointer_[heap_[i]]==i)
    {
      printf("error %i\n",i);
      throw std::runtime_error("heap backpointer inconsistancy");
    }
  }
}

void nm_heap::print() const
{
  for (int i=0; i<heapLength_; i++)
  {
    printf("%i: %i ", i, heap_[i]);
    printf("(%i)",     address_[heap_[i]]);
    printf(" %lf", distance_[heap_[i]]);
    printf(" [%i]\n", backPointer_[heap_[i]]);
  }
  printf("\n");
}
