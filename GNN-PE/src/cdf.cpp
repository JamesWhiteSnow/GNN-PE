#include <math.h>
#include <stdio.h>

long double GL(long double, long double);
long double cdf (long double); 
long double f(long double);
double p=0.;

long double cdf(long double x){
                               if(x>=0.){
                               			  return (1.+GL(0,x/sqrt(2.)))/2.;
                               			 }
                               else {
                               		 return (1.-GL(0,-x/sqrt(2.)))/2.;
                        				}
										}

long double GL(long double a, long double b)
{
   long double y1=0, y2=0, y3=0, y4=0, y5=0;

   long double x1=-sqrt(245.+14.*sqrt(70.))/21., x2=-sqrt(245.-14.*sqrt(70.))/21.;
   long double x3=0, x4=-x2, x5=-x1;
   long double w1=(322.-13.*sqrt(70.))/900., w2=(322.+13.*sqrt(70.))/900.;
   long double w3=128./225.,w4=w2,w5=w1;
   int n=4800;
   long double i=0, s=0, h=(b-a)/n;

   for (i=0;i<=n;i++){
   						y1=h*x1/2.+(h+2.*(a+i*h))/2.;
                     y2=h*x2/2.+(h+2.*(a+i*h))/2.;
                     y3=h*x3/2.+(h+2.*(a+i*h))/2.;
                     y4=h*x4/2.+(h+2.*(a+i*h))/2.;
                     y5=h*x5/2.+(h+2.*(a+i*h))/2.;
                     s=s+h*(w1*f(y1)+w2*f(y2)+w3*f(y3)+w4*f(y4)+w5*f(y5))/2.;
                     }
	return s;
}

long double f(long double x){
										return (2./sqrt(22./7.))*exp(-pow((long double)x,(long double)2.0));
									 }

