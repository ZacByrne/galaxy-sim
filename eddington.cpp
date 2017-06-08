#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

double pi = 3.1415926535;
double Gconst = 6.674*pow(10.0,-11.0);

void makeradius(std::vector<double>& rad_array, unsigned steps, double innerr, double outerr)
{
  double h = (outerr - innerr)/steps;
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    radius = outerr - i*h;
    rad_array[i]= radius;
  }
}


void makerho(std::vector<double>& density_array,std::vector<double>& rad_array, unsigned steps, double galmass, double pluma)
{
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    radius = rad_array[i];
    density_array[i]= 3*galmass / (4*pi*pow(pluma,3.0)) * pow((1 + pow(radius/pluma,2)),-5.0/2.0);
  }
}

void massr(std::vector<double>& density_array, std::vector<double>& mass_array, std::vector<double>& rad_array, unsigned steps, double h)
{
  //double radius = 0;
  double massr= 0;
  for (unsigned i =0; i<steps; ++i)
  {
    //maxradius = rad_array[i];
    massr = 0;
    for (unsigned j =0; j<(steps-i+1); ++j)
    {
      massr = massr + density_array[(steps-j)] * pow(rad_array[(steps-j)],2.0) * (h);
    }
    massr = 4 * pi * massr;
    //massr = massr + 4*pi*pow(rad_array[(steps-1)],3.0)/3*density_array[(steps-1)];
    mass_array[i] = massr;
  }   
}

void potental(std::vector<double>& poten_array, std::vector<double>& mass_array, std::vector<double>& rad_array, unsigned steps, double h)
{
  double potenr = 0;
  for (unsigned i =0; i<steps; ++i)
  {                                                                                                                                                                             
    potenr = 0;
    for (unsigned j =1; j<(steps-i); ++j)
      {
	potenr = potenr + h*mass_array[(steps-j)] / pow(rad_array[(steps-j)],2.0);
	//std::cout << "Potenr  = " << potenr  << "  mass r = " << mass_array[(steps-j)] << "  rad_array  = "  << rad_array[(steps-j)] << std::endl;
      }
    potenr = -1*Gconst * potenr;
    //std::cout << "Potenr  = " << potenr  << "   phin = " << phin << std::endl; 
    poten_array[i] = potenr;
    //poten_array[i] = 
    //poten_array[i] = Gconst * mass_array[i] /(rad_array[i]) + phin;
  }
  
  // Can add -1 at 0
  double phin = poten_array[0];
  for (unsigned i =0; i<steps; ++i)
  {
    poten_array[i] = poten_array[i]-phin;
  }


}
 



void drhodphi(std::vector<double>& poten_array, std::vector<double>& drho_array, std::vector<double>& density_array, unsigned steps)
{
  drho_array[0] = 0;
  drho_array[steps] = 0;
  for (unsigned i =1; i<(steps-1); ++i)
  {
    drho_array[i] = (density_array[(i+1)] -  density_array[(i-1)])/(poten_array[(i+1)] - poten_array[(i-1)]);
  }
  
}

void dtworhodphi(std::vector<double>& poten_array, std::vector<double>& dtworho_array, std::vector<double>& density_array, unsigned steps)
{
  dtworho_array[0] = 0;
  dtworho_array[steps] = 0;
  for (unsigned i =1; i<(steps-1); ++i)
  {
    //dtworho_array[i] = (density_array[(i+1)] +  density_array[(i-1)] - 2*density_array[i])/((poten_array[(i-1)] - poten_array[i]) * (poten_array[i] -  poten_array[(i+1)]));
    dtworho_array[i] = ((density_array[(i+1)] - density_array[i])  - (density_array[i] - density_array[(i-1)]))/((poten_array[(i-1)] - poten_array[i]) * (poten_array[i] -  poten_array[(i+1)]));
  }      
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
  
  //input values
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
  outerr = outerr*3.086*pow(10.0,19.0);
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
  std::vector<double> drho_array(steps);
  std::vector<double> dtworho_array(steps);
  
  
  // need to calc M(<r) from \rho
  double h = (outerr - innerr)/steps;    
  makeradius(rad_array, steps, innerr, outerr);
  makerho(density_array,rad_array, steps, galmass, pluma);
  massr(density_array, mass_array, rad_array, steps, h);
  potental(poten_array, mass_array, rad_array, steps, h);
  drhodphi(poten_array, drho_array, density_array, steps);
  dtworhodphi(poten_array, dtworho_array, density_array, steps);
  
  
  std::ofstream myfile;
  myfile.open (filename.c_str());
  //double h = (outerr - innerr)/steps;
  //r = 0;
  double poten = 0;
  double massen = 0;
  double potennewt = 0;
  double firstder = 0;
  double secder = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    poten =  1 * Gconst * galmass/pluma * pow((1 + pow(rad_array[i]/pluma,2)),-1.0/2.0) + poten_array[0];
    massen =  galmass * pow((1 + pow(rad_array[i]/pluma,-2.0)),-3.0/2.0);
    potennewt = Gconst * mass_array[i] / pow(rad_array[i],2.0);
    firstder = -15*pow(pluma,2.0)*pow((poten_array[i]),4) / (4*pi*pow(galmass,4.0) * pow(Gconst,5.0));
    secder = 15*pow(pluma,2.0)*pow((poten_array[i]),3) / (pi*pow(galmass,4.0) * pow(Gconst,5.0));
    myfile << rad_array[i] << "   " << density_array[i] << "   "<< mass_array[i]  <<"  " << poten_array[i] << "   " << drho_array[i] << "   " << dtworho_array[i] << "    "<< poten << "   " << massen << "   " << potennewt << "   " << firstder << "    " << secder  << "\n";
    
  }

  myfile.close();

  return 0;

}
