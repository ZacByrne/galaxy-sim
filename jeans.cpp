#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include <cassert>
//#include <boost/lexical_cast.hpp>

// NULL is gone, use nullptr
//double* me = nullptr;
//double pluma = 1;
double pi = 3.1415926535;
//GET VALUE FOR HOLGER

//double galmass = 1:
double Gconst = 6.674*pow(10.0,-11.0);

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
double vsteldelplum(double r, double galmass, double pluma)
{
  double rho = 3*galmass / (4*pi*pow(pluma,3)) * pow((1 + pow(r/pluma,2)),-5/2);
  return rho;
}

double dvstelplum(double r, double galmass, double pluma)
{
  double drho = -15.0/4.0 * galmass*r/(pi*pow(pluma,5)) * pow((1 + pow(r/pluma,2)),-7.0/2.0);
  return drho;
}

//Mass Function
double Mfuncplum(double r, double galmass, double pluma)
{
  double massr = galmass * (1 + pow(pluma/r,2));
  return massr;
}

//derivative of gfravitational potential
double dgravpotplum(double r, double galmass, double pluma)
{
  double dphi = Gconst * galmass * r * pow(pluma, -3) * pow((1 + pow(r/pluma,2)),-3.0/2.0);
  return dphi;
}

//function of r and yn for midpoint
double fry(double r, double y_n, double galmass, double blackmass, double pluma)
{
  double fvalue = -1*dgravpotplum(r,galmass, pluma) - Gconst*blackmass/r - 2*beta(y_n,rotvel(r))* y_n / r - dvstelplum(r, galmass, pluma)*y_n/vsteldelplum(r, galmass, pluma);
  return fvalue;
}

//double integrate(double start, double end, double step)
//{
 // double s = 0;
//  double h = (end - start)/steps;
//  for (int i =0; i<steps; ++i)
//  {
 //   //s += h* 
//  }
//}

void midpointarray(std::vector<double>& v_rarray, double start, double end, unsigned step, double galmass, double blackmass)
{
  double s =0;
  double h = (end - start)/step;
  double yn = s;
  double rn = 0;
  for (unsigned i =0; i<step; ++i)
   {
      yn = s;
      v_rarray[i]=yn;
      rn = start + i*h;
      s = s + h/2*fry(rn,yn, galmass, blackmass, end);
      s = yn + h*fry(rn+h/2, s, galmass, blackmass, end);
      
      // changed to not recursive 
        //s = s + h/2*fry((rn+h/2), (yn + h/2*fry(rn,yn)));
      
   }
}



int main(int argc, char** argv)
{
  //input from user?
  //steps, start, end?
  double pluma;
  double start;
  unsigned steps;
  double galmass;
  double blackmass;
  
  using std::cout;
  using std::cin;
  
  //input values
  cout << "Jean's equation simulation. \n Please input core radius (Kpc)." << std::endl;
  cin >> pluma;
  pluma = pluma*3.086*pow(10.0,19.0);
  //convert to metres
    
  cout << "Please input start radius. (Kpc)" << std::endl;
  cin >> start;
  start = start*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please input no. steps" << std::endl;
  cin >> steps;

  cout << "Please input Mass of Galaxy (solar masses)" << std::endl;
  cin >> galmass;
  galmass = galmass * 1.989 * pow(10,30.0);
    
  cout << "Please input additional black hole mass (0 for none)(solar masses)" << std::endl;
  cin >> blackmass;
  blackmass = blackmass* 1.989 * pow(10,30.0);
  
  std::string filename;  
  cout << "Please input output filename" << std::endl;
  cin >> filename;
  
  //Check with holger about end of int being radius
  //end = pluma;
  
  //vector array
  std::vector<double> vrarray(steps);
  
  midpointarray(vrarray, start, pluma, steps, galmass, blackmass);
  
  std::ofstream myfile;
  myfile.open (filename.c_str());
  double h = (pluma-start)/steps;
  double r = 9;
  double sigma2 = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    r = start + i*h;
    sigma2 = Gconst*galmass / (6*pluma) * pow((1 + pow(r/pluma,2)),-1.0/2.0);
    myfile << r << "   " << vrarray[i] << "   "<< sigma2  << "\n";
  }
  myfile.close();
  return 0;

}


