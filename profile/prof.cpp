#include "fast_marching.h"
#include "math.h"

int main(void)
{
  int size = 1500;
  int dim  = 2;
  double dx[] = { 0.1, 0.1, 0.1 };
  int shape[3];
  shape[0]=size;
  shape[1]=size;
  shape[2]=size;

  int lsize =pow(size,dim);
  double *phi         = new double[lsize];
  double *distance    = new double[lsize];
  int    *flag        = new int[lsize];

  for (int i=0; i<lsize; i++)
  {
    phi[i]=1.0;
    distance[i]=0.0;
    flag[i]=0;
  }
  phi[0]=-1;

  fastMarcher *fm = new fastMarcher(
    phi,dx,flag,0,distance,2,shape,false);

  delete [] phi;
  delete [] distance;
  delete [] flag;
  delete fm;
  phi = distance = flag = fm = 0;

  return 0;
}

