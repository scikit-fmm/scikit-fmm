//no_malloc_heap.h

class nm_heap
{
 public:
  nm_heap();
  virtual ~nm_heap();
  void       set_data(int depth, double *hd, long *hi1, long *hi2,
                 long *hi3, bool self_Test);
  int        push(int address, double value);
  void       pop(int *address, double *value);
  void       set(int index, double value);
  bool       empty() const;

 private:
  void       print() const;
  void       test() const;
  inline void       _siftUp(int pos);
  inline void       _siftDown(int startPos, int pos);

  int        maxLength_;
  int        listLength_;
  int        heapLength_;

  double    * distance_;
  long      * heap_;
  long      * address_;
  long      * backPointer_;
  bool        selfTest_;
};
