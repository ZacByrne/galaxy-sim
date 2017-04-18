//#include <iostream>
//#include <fstream>
#include <cmath>
//#include <cassert>
//#include <boost/lexical_cast.hpp>

// NULL is gone, use nullptr
double* me = nullptr;
double pluma = 1;
//GET VALUE FOR HOLGER

double galmass = 1:
double Gconst = 1;


double beta(double v_r2, double v_theta2)
{
  betavalue = 1 - (v_theta2)/(v_r2);
  return betavalue;
}

double  rotvel(double r)
{
  //WHAT IS THE ROT VEL?
  return 0;
}


double fry(double r, double y_n)
{
  fvalue = dgravpotplum(r)/vsteldelplum(r) - 2*beta(y_n,rotvel(r))* y_n / r - dvelstelplum*y_n/vsteldelplum(r);
  return fvalue
}

//stellar density
double vsteldelplum(double r)
{
  //need density equation?
}

double dvstelplum(double r)
{
  //need density equation?
}

//Mass Function
double Mfuncplum(double r)
{
  massr = galmass * (1 + pow(pluma/r,2));
  return massr
}


double dgravpotplum(double r, )
{
  dphi = Gconst * Mfuncplum(r) * r * pow(pluma, -3) * (1 + pow(pluma/r,2));
  return dphi
}


double integrate(double start, double end, double step)
{
  double s = 0;
  double h = (end - start)/steps;
  for (int i =0; i<steps; ++i)
  {
    s += h* 
  }
}

double midpointarray(double start, double end, double step)
{
  float s =0;
  float h = (end - start)/steps;
  float yn = s;
  for (int i =0; i<steps; ++i)
   {
      yn = s;
      rn = i*h
      s = h*fry(rn+h/2, yn + h/2*fry(rn,yn));
      
   }


int main()
{
  return 0;

}


