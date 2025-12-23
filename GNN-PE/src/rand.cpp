#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rand.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

double gaussian (double mean, double sigma)
{
   double v1, v2;
   double s;
   double x;

   do
   {
      v1 = 2*uniform(0, 1) - 1;
      v2 = 2*uniform(0, 1) - 1;
      s = v1*v1 + v2*v2;
   }
   while (s >= 1.);

   x = v1 * sqrt ( -2. * log (s) / s);
   
   x = x * sigma + mean;

   return (x);
}

double uniform(double _min, double _max)
{

	int int_r = rand();
	long base = RAND_MAX-1;
	double f_r  = ((double) int_r) / base;
	double rlt =(_max - _min) * f_r + _min;

	if (rlt>_max)
		rlt=_max;
	if (rlt<_min)
		rlt=_min;

	return rlt;
}

double zipf(double x1, double x2, double p)
{

   double x;
   double i;
   double r, HsubV, sum;
   int V = 100;

   HsubV = 0.0;
   for(i=1; i<=V; i++) 
      HsubV += 1.0/pow( (double)i, p);

   r = uniform(0., 1.)*HsubV;
   sum = 1.0; i=uniform(0,1);
   while( sum<r){
	  i+=uniform(1, 2);
      sum += 1.0/pow( (double)i, p);
   }

   x = ( (double) i - 1. ) / ( (double) V - 1. );
   x = (x2 - x1) * x + x1;

   return(x);
}

