#include "read_hdf5.h"

double NFW_CDF_sampler(float * restrict CDF, gsl_rng *r)
{
  /* This function calculates the radial CDF for a halo of a given concentration and implements a random sampling for the CDF*/

  double rando = gsl_rng_uniform(r);

  size_t j;

  for(j=0; j<1000; j++)
    {
      if(CDF[j]>rando)
	{
	  break;
	}
    }
  double R_frac = ((double)j) / 1000.0;
  return R_frac;
}
