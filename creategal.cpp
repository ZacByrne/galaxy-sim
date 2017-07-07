#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

double pi = 3.1415926535;
double Gconst = 1.0;// 6.674*pow(10.0,-11.0);    

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i,k;
  double p,qn,sig,un,u[n],lx[n],ly[n];

  for (i=0;i<n;i++) {
    lx[i+1]=x[i];
    ly[i+1]=y[i];
  }

  if (yp1 > 0.99e30)
    y2[1]=u[0]=0.0;
  else {
    y2[1] = -0.5;
    u[0]=(3.0/(lx[2]-lx[1]))*((ly[2]-ly[1])/(lx[2]-lx[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(lx[i]-lx[i-1])/(lx[i+1]-lx[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i-1]=(ly[i+1]-ly[i])/(lx[i+1]-lx[i]) - (ly[i]-ly[i-1])/(lx[i]-lx[i-1]);
    u[i-1]=(6.0*u[i-1]/(lx[i+1]-lx[i-1])-sig*u[i-2])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(lx[n]-lx[n-1]))*(ypn-(ly[n]-ly[n-1])/(lx[n]-lx[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k-1];
}

double splint(double xa[], double ya[], double y2a[], int n, double x)
{
  int klo,khi,k;
  double h,b,a;
  //std::cout << "test splint n  = "<< n <<std::endl;                                                                                                                                                      \
                                                                                                                                                                                                            
  klo=0;
  khi=n;
  while (khi-klo > 1) {
    //std::cout << "khi = " << khi << std::endl;                                                                                                                                                           \
                                                                                                                                                                                                            
    k=(khi+klo) >> 1;
    //std::cout << "k >> 1 = " << k << std::endl;                                                                                                                                                          \
                                                                                                                                                                                                            
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  //std::cout << khi << "   klo " << klo << std::endl;                                                                                                                                                     \
                                                                                                                                                                                                            
  h=xa[khi]-xa[klo];
  if (h == 0.0) std::cout<<"Bad xa input to routine splint\n";

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void generatestars(std::vector<double> mass_array, std::vector<double> rad_array, std::vector<double> draddmasso_array, unsigned starno, std::vector<double> draddmass_array,std::vector<double>& starrad_array)
{
  double randmass;
  double radius;
  

  for (unsigned i=0; i < (starno+1); ++i)
    {
      randsmass = rand()/(RAND_MAX); 
      radius = splint(mass_array.data(), rad_array.data(), draddmass_array.data(), rad_array.size(), (randmass));
      starrad_array[i]=radius;
    }

}

void starxyz(


int main() {
  unsigned steps = 200000;
  std::vector<double> mass_array;
  std::vector<double> rad_array;
  std::vector<double> distri_array;


  std::string filename;
  cout << "Please filename of galaxy data file" << std::endl;
  cin >> filename;

  std::ifstream infile(filename);
  
  float radius, density, mass, poten, drho, dtworho, analypot, analymass, distri, andrho, andtworho, andistri;
  

  //read in file
  while (infile >> radius >> density >> mass >> poten >> drho >> dtworho >> analypot >> analymass >> distri >> andrho >> andtworho >> andistri)
    {
      rad_array.push_back(radius);
      mass_array.push_back(mass);
      distri_array{distri);
    }

  //scale mass (should be 1 anyway)
  double massscale = mass_array[(mass_array.size())];
  for (unsigned i = 0; i<(mass_array.size()), ++i)
    {
      mass_array[i] = mass_array[i]/massscale;
    }
  

  unsigned starno;
  cout << "How many Stars" << std::endl;
  cin >> starno;

  std::vector<double> draddmass_array((rad_array.size()));

  double dyn1 = ;
  double dyn2 = ;

  spline(mass_array.data(), rad_array.data(), (rad_array.size()), dyn1, dyn2, draddmass_array.data());

  std::vector<double> starradius_array;

  generatestars

  

  

}
