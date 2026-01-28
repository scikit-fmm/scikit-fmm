//travel_time_marcher_genes.cpp

#include "travel_time_marcher_genes.h"
#include "math.h"
#include "heap.h"
#include <stdexcept>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element

using std::printf;

void travelTimeMarcherGenes::initalizeFrozen()
{
  distanceMarcher::initalizeFrozen();
  for (int i=0; i<size_; i++)
  {
    printf("speeds_ = %x\n", speeds_);
    printf("speeds_[%d] = %g\n", i, speeds_[i]);
    printf("distance_[%d] = %g\n", i, distance_[i]);
    if (flag_[i]==Frozen)
    {
      // convert distance to time
      distance_[i] = fabs(distance_[i]/speeds_[i]);
    }
    printf("distance_[%d] = %g\n", i, distance_[i]);
  }
}

double travelTimeMarcherGenes::updatePointOrderTwo(int i)
{
    double res = updatePointOrderTwo(i, std::set<int>());
    if (res == std::numeric_limits<double>::infinity()) {
        throw std::runtime_error("Unreachable voxel");
    } else {
        return res;
    }
}


// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double travelTimeMarcherGenes::updatePointOrderTwo(int i, std::set<int>avoid_dim)
{
  double a,b,c;
  a=b=c=0;
  int naddr, naddr2; // addresses of neighbours
  unsigned naddr_smallest_nbr = i; // set a reasonable default value
  // Choose a "good" pair of neighbours on different axes:
  for (int dim=0; dim<dim_; dim++) {
    if (avoid_dim.find(dim) != avoid_dim.end()) {
      continue; //we should avoid this dimension
    }
    double value1 = maxDouble;
    double value2 = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value1))
        {
          value1 = distance_[naddr];
          naddr2 = _getN(i,dim,j*2,Mask);
          if (naddr2!=-1 &&
              flag_[naddr2]==Frozen &&
              ((distance_[naddr2]<=value1 && value1 >=0) ||
               (distance_[naddr2]>=value1 && value1 <=0)))
          {
            value2=distance_[naddr2];
            if (phi_[naddr2] * phi_[naddr] < 0  || phi_[naddr2] * phi_[i] < 0)
              value2 *= -1;
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
    // note the neighbour with the smallest phi/distance value:
    if (naddr != -1) {
      naddr_smallest_nbr = naddr;
        if ((naddr2 != -1) && (value1 > value2)) naddr_smallest_nbr = naddr2;
    }
  }
  printf("moe\n");
  // set an initial value for branch function at node i:
  printf("branch_[%d] = %d\n", i, branch_[i]);
  printf("drivers_[%d] = %d\n", i, drivers_[i]);
  printf("naddr_smallest_nbr = %d\n", naddr_smallest_nbr);
  if (naddr_smallest_nbr != -1) branch_[i] = branch_[naddr_smallest_nbr];
  printf("catch a\n");
  try {
    double res = solveQuadratic(i,a,b,c);
    // update branch function if a mutation is present at naddr
    // AND the mutation is not already accounted for
    if ((drivers_[i] > 0) && ((branch_[i] & drivers_[i]) == 0)) {
        branch_[i] += drivers_[i];
        printf("driver %d added at position %d, branch is now %d\n", drivers_[i], i, branch_[i]);
    }
    //TODO 28-01-26 - we're updating the branch function just on at the point we're on and not further along
    // need to do something with adding the branch to downstream unfrozen points
    printf("tiger\n");
    return res;
  } catch (std::runtime_error & err) {
    //if the determinant is negative, we try to reach the voxel with one dimension less and take the minimum
    if (avoid_dim.size() == (size_t)dim_) {
      //end of the recursion, use inf so that it is discarded selecting the minimum
      return std::numeric_limits<double>::infinity(); 
    }
    std::vector<double> sols;
    for (int ind=0; ind<dim_; ind++){
      //remove one dimension more than what we are already doing
      std::set<int> tempset = avoid_dim;
      std::pair<std::set<int>::iterator, bool> ret = tempset.insert(ind);
      //avoid recursive call on identical parameters (the set already had *ind* in it):
      if(!ret.second) continue;
      sols.push_back(updatePointOrderTwo(i,tempset));
    }
    if (sols.size()==0) {
      return std::numeric_limits<double>::infinity();
      //All the derivates with different dimensionalities are 0
    }
    // update branch function if a mutation is present at naddr
    if ((drivers_[i] > 0) && ((branch_[i] & drivers_[i]) == 0)) {
        branch_[i] += drivers_[i];
    }
    return *std::min_element(sols.begin(), sols.end());
  }
}


double travelTimeMarcherGenes::solveQuadratic(int i, const double &a,
                                         const double &b,
                                         double &c)
{
  unsigned bvalue = branch_[i];
  c -= 1/pow(speeds_[bvalue*size_+i],2); // previously speeds_[bvalue][i]:
  // change to something like speeds_[index(branch, i)]
  double r0 = 0;
  double det = pow(b, 2) - 4 * a * c;
  if (det >= 0)
  {
    r0 = (-b + sqrt(det)) / 2.0 / a;
  }
  else
  {
    throw std::runtime_error("Negative discriminant in (genetic) time marcher quadratic.");
  }
  return r0;
}
