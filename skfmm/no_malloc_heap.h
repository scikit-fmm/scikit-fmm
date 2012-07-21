//no_malloc_heap.h

class nm_heap
{
 public:
  nm_heap();
  virtual ~nm_heap();
  void       set(int depth, double *hd, int *hi1, int *hi2,
                 int *hi3, bool self_Test);
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

  double   * distance_;
  int      * heap_;
  int      * address_;
  int      * backPointer_;
  bool       selfTest_;
};
