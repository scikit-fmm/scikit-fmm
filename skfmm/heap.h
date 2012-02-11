// heap.h
/**************************************************************************/
/* This class implements a binary min heap data structure to support the  */
/* fast marching algorithm. A min heap is a list which has the property	  */
/* that the smallest element is always the first element.		  */
/* 									  */
/* The fast marching method uses this data structure to track elements in */
/* the solution narrow band. The fast marching method needs to know which */
/* element in the narrow band is nearest the zero level-set at each	  */
/* iteration.								  */
/* 									  */
/* When a new point enters the solution narrow band the element is added  */
/* to the heap.								  */
/* 									  */
/* New elements are added to the heap with the push() method. The address */
/* passed to pop is the (flat) address of the element in the grid and the */
/* value passed is the distance. The push() method returns an integer (a  */
/* heap index) which is used later to refer to the element.		  */
/* 									  */
/* As the solution evolves the distance of points already in the heap can */
/* be updated via the set() method. The set() method takes the heap index */
/* returned by the push() along with a new distance value.		  */
/* 									  */
/* The narrow band element nearest the zero-level set is taken off the	  */
/* heap with the pop() method. The grid address of the top element is	  */
/* returned to the caller.						  */
/* 									  */
/* The constructor for heap needs to know the number of elements that	  */
/* will enter the narrow band. See heap.cpp for implementation details.	  */
/**************************************************************************/

class heap
{
 public:
  heap(int depth, bool selfTest=false);
  virtual ~heap();
  int        push(int address, double value);
  void       pop(int *address, double *value);
  void       set(int index, double value);
  bool       empty() const;

 private:
  void       print() const;
  void       test() const;
  void       _siftUp(int pos);
  void       _siftDown(int startPos, int pos);

  int        maxLength_;
  int        listLength_;
  int        heapLength_;

  double   * distance_;
  int      * heap_;
  int      * address_;
  int      * backPointer_;
  bool       selfTest_;
};
