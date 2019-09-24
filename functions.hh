
#include <iostream>
#include "TMath.h"




double fitBarEffErr(double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  
  double t = x[0];
  
  double center		= par[0];
  double width 		= par[1];
  double baseline	= par[2];
  double max     	= par[3];
  double sigma_x   	= par[4];
  
  double f;
  double edge1 = center-width/2;
  double edge2 = center+width/2;
  
  if ( t <= center)
  {
      f = baseline + max/2*(erf( (t-edge1) / (sqrt(2)*sigma_x)) + 1);
  }
  if (t > center)
  {
      f = baseline + max/2*(erf((-t+edge2) / (sqrt(2)*sigma_x)) + 1);
  }
  
  
  return f;
}
