/* $Id: randomF77.c 95 2002-11-22 13:24:41Z hothorn $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include <R.h>
#include<Rmath.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }


double F77_SUB(phid)(double *x){
  return pnorm(*x, 0.0, 1.0, 1, 0);
  
}

double F77_SUB(studnt)(int *nu, double *x){
  return pt(*x, *nu, 1, 0);
}
