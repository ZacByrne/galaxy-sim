//#include <iostream>
//#include <fstream>
#include <cmath>
#include <vector>
//#include <cassert>
//#include <boost/lexical_cast.hpp>

// NULL is gone, use nullptr
double* me = nullptr;
//double pluma = 1;
double pi = 3.1415926535;
//GET VALUE FOR HOLGER

//double galmass = 1:
double Gconst = 6.674*pow(10,-11);

//rotational velocity
double  rotvel(double r)
{
  //WHAT IS THE ROT VEL?
  return 1;
}

//beta
double beta(double v_r2, double v_theta2)
{
  double betavalue = 0;
    //1 - (v_theta2)/(v_r2);
  return betavalue;
}

//stellar density
double vsteldelplum(double r)
{
  double rho = 3*galmass / (4*pi*pow(pluma,3)) * pow((1 + pow(pluma/r,2)),-5/2);
  return rho;
}

double dvstelplum(double r)
{
  double drho = -15/4 * galmass/(pi*pow(pluma,5)) * pow((1 + pow(pluma/r,2)),-7/2)
}

//Mass Function
double Mfuncplum(double r)
{
  double massr = galmass * (1 + pow(pluma/r,2));
  return massr;
}

//derivative of gfravitational potential
double dgravpotplum(double r, )
{
  dphi = Gconst * Mfuncplum(r) * r * pow(pluma, -3) * pow((1 + pow(pluma/r,2)),-3/2);
  return dphi;
}

//function of r and yn for midpoint
double fry(double r, double y_n)
{
  fvalue = -1*dgravpotplum(r) - 2*beta(y_n,rotvel(r))* y_n / r - dvelstelplum*y_n/vsteldelplum(r);
  return fvalue;
}

double integrate(double start, double end, double step)
{
  double s = 0;
  double h = (end - start)/steps;
  for (int i =0; i<steps; ++i)
  {
    //s += h* 
  }
}

void midpointarray(std::vector<double>& v_rarray, double start, double end, double step)
{
  float s =0;
  float h = (end - start)/steps;
  float yn = s;
  for (int i =0; i<steps; ++i)
   {
      yn = s;
      v_rarry[i]=yn;
      rn = start + i*h;
      s = s + h/2*fry(rn,yn);
      s = yn + h*fry(rn+h/2, s);
      
      // changed to not recursive 
        //s = s + h/2*fry((rn+h/2), (yn + h/2*fry(rn,yn)));
      
   }


int main(int argc, char** argv)
{
  //input from user?
  //steps, start, end?
  double pluma;
  double start;
  double steps;
  double galmass;
  double blackmass;
  
  //input values
  cout >> "Jean's equation simulation. \n Please input core radius (double).";
  cin >> pluma;
  
  cout >> "Please input start radius.";
  cin >> start;
  
  cout >> "Please input no. steps";
  cin >> steps;

  cout >> "Please input Mass of Galaxy";
  cin >> galmass;
  
  cout >> "Please input additional black hole mass (0 for none)";
  cin >> blackmass;
  
  //Check with holger about end of int being radius
  end = pluma;
  
  //vector array
  std::vector<int> vrarray(steps);
  
  midpointarray(vrarray, start, end, steps, blackmass
  
  
  

}


