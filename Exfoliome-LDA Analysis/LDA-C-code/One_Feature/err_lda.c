/*
 * file : err_lda.c
 *
 * summary : Exact computation of
 *   the true error of a designed LDA 
 *   classifier, based on a model
 *   where each class-conditional
 *   density is a spherical Gaussian
 *   distribution, with means at 
 *   (-dt,...,-dt) for class 1 and
 *   (dt,...,dt) for class 2,
 *   and standard deviations s1 and s2,
 *   respectively. The LDA classifier is
 *   assumed to have been designed from
 *   samples drawn from this model.
 *
 * function call :
 * 
 *   e = err_lda(p,p1,dt,s1,s2,a,m);
 *     double e -> true error
 *     int p -> dimensionality
 *     double p1 -> prior probability of first class
 *     double dt -> means at (dt,...,dt) and (-dt,...,-dt)
 *     double s1 -> std for class 1
 *     double s2 -> std for class 2
 *     double *a,m -> LDA coefficients (y = a^T x + m)
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 01/15/2004 
 */

#include <math.h>
#include "err_est.h"

double err_lda(int p,
	         double p1,
	         double dt,
	         double s1,
	         double s2,
	         double *a,
	         double m)
{
  double n=0,v=0,x1,x2;
  int i;

  for (i=0; i<p; i++)
    v += a[i];

  for (i=0; i<p; i++)
    n += a[i]*a[i];
  n = sqrt(n);

  dt = -dt;
  x1 = (v*dt+m)/(n*s1);
  x2 = (v*dt-m)/(n*s2);

  /* error expression */
  return p1*ncdf(x1) + (1-p1)*ncdf(x2);
}
