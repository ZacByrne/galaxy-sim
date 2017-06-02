#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

double pi = 3.1415926535;
double Gconst = 6.674*pow(10.0,-11.0);

void makeradius(std::vector<double>& rad_array, double steps, double innerr, double outerr)
{
  double h = (outerr - innerr)/steps;
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    radius = outerr - i*h;
    rad_rarray[i]= radius;
  }
}


void makerho(std::vector<double>& density_array,std::vector<double>& rad_array, double steps, double galmass, double pluma)
{
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    radius = rad_array[i];
    density_rarray[i]= 3*galmass / (4*pi*pow(pluma,3.0)) * pow((1 + pow(radius/pluma,2)),-5.0/2.0);
  }
}

void massr(std::vector<double>& density_array, std::vector<double>& mass_array, std::vector<double>& rad_array, double steps)
{
  double radius = 0;
  double massr= 0;
  for (unsigned i =0; i<steps; ++i)
  {
    maxradius = rad_array[i];
    massr = 0;
    for (unsigned j =0; j<(steps-i); ++j)
    {
      massr = massr + density_array[(steps-j)] * pow(rad_array[(steps-j)],2);
    }
    massr = 4 * pi * massr;
    mass_array[i] = massr;
  }   
}

void potental(std::vector<double>& poten_array, std::vector<double>& mass_array, std::vector<double>& rad_array, double steps)
{
  for (unsigned i =0; i<steps; ++i)
  {
    poten_array[i] = -1* Gconst * mass_array[i] / rad_array[i];
  }
  // Can remove -1 if relative potential?
}
  

int main(int argc, char** argv)
{
  double pluma;
  double outerr;
  unsigned steps;
  double galmass;
  double blackmass;
  double innerr;
  //double beta;
  
  using std::cout;
  using std::cin;
  
  /input values
  cout << "Jean's equation simulation. \n Please input core radius (Kpc)." << std::endl;
  cin >> pluma;
  pluma = pluma*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please inner integration radius (Kpc) (smallest radius for program) " << std::endl;
  cin >> innerr;
  innerr = innerr*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please input out integration radius. (Kpc) (similar to radius of galaxy)" << std::endl;
  cin >> outerr;
  start = start*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please input no. steps" << std::endl;
  cin >> steps;

  cout << "Please input Mass of Galaxy (solar masses)" << std::endl;
  cin >> galmass;
  galmass = galmass * 1.989 * pow(10.0,30.0);
  //convert into kg's
    
  cout << "Please input additional black hole mass (0 for none)(solar masses)" << std::endl;
  cin >> blackmass;
  blackmass = blackmass* 1.989 * pow(10.0,30.0);
  //convert into kg's
  
  //cout << "Please input Beta (anisotropy value)" << std::endl;
  //cin >> beta;
    
  
  std::string filename;  
  cout << "Please input output filename" << std::endl;
  cin >> filename;
  
  std::vector<double> density_array(steps);
  std::vector<double> rad_array(steps);
  std::vector<double> mass_array(steps);
  std::vector<double> poten_array(steps);
  
  
  // need to calc M(<r) from \rho


}
