/*
 * function: min_dist
 *
 * summary : computes the minimum distance
 * from a given point to the a data set
 *
 * function call :
 * 
 *   d = min_dist(x,X,n,p);
 *     double d -> minimum distance
 *     double *x -> data point 
 *     double *X -> data set
 *     int n -> number of observations in X
 *     int p -> dimension of observations
 *
 * author: Ulisses Braga-Neto
 *   
 * last revision: 05/21/2003 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "err_est.h"

double min_dist(double *x,
		double *X,
		int n,
		int p)
{
  register int i,j;
  double d,dm,v;

  dm = HUGE_VAL;
  
  for (i=0; i<n; i++) {

    for (d=j=0; j<p; j++) {
      v = x[j]-X[i*p+j];
      d += v*v;
    }

    if (d<dm) dm = d; // new minimum
  }

  return sqrt(dm);
}
