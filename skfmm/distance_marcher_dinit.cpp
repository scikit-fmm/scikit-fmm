#include "distance_marcher_dinit.h"

void distanceMarcherDInit::initalizeFrozen()
{
  // this is a special case of the init frozen because the distance to
  // the zero contour of the initially frozen points has been pre-calculated
  // in the python wrapper.
  int c=0;
  for (int i=0; i<size_; i++)
  {
    if (! (dinit_[i] == maxDouble))
    {
      c++;
      flag_[i] = Frozen;
      distance_[i] = dinit_[i];
    }
  }
}
