//travel_time_marcher_genes.h
#include "distance_marcher.h"
#include <set>

class heap;

class travelTimeMarcherGenes : public distanceMarcher
{
public:
  travelTimeMarcherGenes(double* phi,      double* dx, long long int* flag,
                         double* distance, int ndim,   int* shape,
                         bool self_test,   int order,
                         double* drivers, double* speeds, double narrow,
                         int periodic) :
    distanceMarcher(phi, dx, flag, distance, ndim, shape, self_test,
                    order, narrow, periodic)
    // TODO implement passing the speeds array (it has one more dimension than
    // the speed array in travel_time_marcher).
  {
  }

  virtual ~travelTimeMarcherGenes() { }

protected:
  virtual void   initalizeFrozen();
  virtual double updatePointOrderTwo(int i);
  virtual double updatePointOrderTwo(int i, std::set<int> avoid_dim);
  virtual double solveQuadratic(int i, const double &a,
                                const double &b, double &c);
private:
  double* speed_;
};
