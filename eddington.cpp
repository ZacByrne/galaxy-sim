#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

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
  //std::cout << "test splint n  = "<< n <<std::endl;                                                                                                                                                       
  klo=0;
  khi=n;
  while (khi-klo > 1) {
    //std::cout << "khi = " << khi << std::endl;                                                                                                                                                            
    k=(khi+klo) >> 1;
    //std::cout << "k >> 1 = " << k << std::endl;                                                                                                                                                           
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  //std::cout << khi << "   klo " << klo << std::endl;                                                                                                                                                      
  h=xa[khi]-xa[klo];
  if (h == 0.0) std::cout<<"Bad xa input to routine splint\n";

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


void makelogradius(std::vector<double>& rad_array, unsigned steps, double innerr, double outerr)
{
  double loginr = log10(innerr);
  double logoutr = log10(outerr); 
  double h = (logoutr - loginr)/steps;
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    
    radius = logoutr - i*h;

    rad_array[i]= pow(10.0, radius);
  }
}

void makeradius(std::vector<double>& rad_array, unsigned steps, double innerr, double outerr)
{
  // double innerr = (innerr);
  //double outerr = (outerr);
  double h = (outerr - innerr)/steps;
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
    {

      radius = outerr - i*h;

      rad_array[i]=  radius;
    }
}



void makerho(std::vector<double>& density_array,std::vector<double>& rad_array, unsigned steps, double galmass, double pluma)
{
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
  {
    radius = rad_array[i];
    density_array[i]= 3.0*galmass / (4*pi*pow(pluma,3.0)) * pow((1.0 + pow(radius/pluma,2)),-5.0/2.0);
  }
}

void massr(std::vector<double>& density_array, std::vector<double>& mass_array, std::vector<double>& rad_array, unsigned steps, double h)
{
  //double radius = 0;
  double massr= 0.0;
  for (unsigned i =0; i<steps; ++i)
  {
    //maxradius = rad_array[i];
    massr = 0.0;
    for (unsigned j =2; j<(steps-i+1); ++j)
    {
      massr = massr + density_array[(steps-j+1)] * pow((rad_array[(steps-j)]+rad_array[(steps-j+1)])/2,2.0) * ((rad_array[(steps-j)]-rad_array[(steps-j+1)]));
    }
    massr = 4.0 * pi * massr + 4.0*pi* density_array[(steps-1)] * pow(rad_array[(steps-1)],3.0)/3.0;
    //massr = massr + 4*pi*pow(rad_array[(steps-1)],3.0)/3*density_array[(steps-1)];
    mass_array[i] = massr;
    //std::cout << "test" <<std::endl;
  }   
}

double potental(std::vector<double>& poten_array, std::vector<double>& mass_array, std::vector<double>& rad_array, unsigned steps, double h)
{
  double potenr = 0.0;
  for (unsigned i =0; i<steps; ++i)
  {                                                                                                                                                                             
    potenr = 0.0;
    for (unsigned j =1; j<(steps-i); ++j)
      {
	potenr = potenr + (rad_array[(steps-j)]-rad_array[(steps-j+1)])*mass_array[(steps-j)] / pow(rad_array[(steps-j)],2.0);
	//std::cout << "Potenr  = " << potenr  << "  mass r = " << mass_array[(steps-j)] << "  rad_array  = "  << rad_array[(steps-j)] << std::endl;
      }
    potenr = potenr + mass_array[(steps-1)]/rad_array[(steps-1)];
    potenr = -1*Gconst * potenr;
    //std::cout << "Potenr  = " << potenr  << "   phin = " << phin << std::endl; 
    poten_array[i] = potenr;
    //poten_array[i] = 
    //poten_array[i] = Gconst * mass_array[i] /(rad_array[i]) + phin;
  }
  
  // Can add -1 at 0
  double phin = poten_array[0];
  std::cout << "phin = " << phin << std::endl;
  for (unsigned i =0; i<steps; ++i)
  {
    poten_array[i] = poten_array[i]-phin;
  }
  return phin;

}
 


//first derivative 
void drhodphi(std::vector<double>& poten_array, std::vector<double>& drho_array, std::vector<double>& density_array, unsigned steps)
{
  drho_array[0] = 0.0;
  //drho_array[steps] = 0;
  for (unsigned i =2; i<(steps-1); ++i)
  {
    drho_array[i] = (density_array[i-1] -  density_array[(i+1)])/(poten_array[(i-1)] - poten_array[(i+1)]);

    //std::cout << "density_array[(i+1)] = "<<  density_array[i+1] << "   den i-1  = " << density_array[(i-1)] << " pot i+1  = " <<  poten_array[(i+1)] << "  pot i-1  = " <<  poten_array[(i-1)] << "  dif top + " << (density_array[(i+1)] -  density_array[(i-1)]) << "  dif bot  =" << (poten_array[(i+1)] - poten_array[(i-1)]) << std::endl;
  }
  drho_array[steps-1] = drho_array[steps-3];
  drho_array[steps-2] = drho_array[steps-3];
}

//void dtest(std::vector<double>& drho_array, std::vector<double> density_array, std::vector<double> rad_array,  unsigned steps, double pluma,double galmass)
//{
  //drho_array[0] = 0.0;
  //for ( unsigned i =2; i<(steps-1); ++i)
    //{
      //  drho_array[i] = (density_array[(i+1)] -  density_array[(i-1)])/(Gconst * galmass/pluma * pow((1.0 + pow(rad_array[i+1]/pluma,2)),-1.0/2.0) - Gconst * galmass/pluma * pow((1.0 + pow(rad_array[i-1]/pluma,2.0)),-1.0/2.0));  
      // }
  //drho_array[steps-1] = drho_array[steps-2];

  //} 


//second derviative
void dtworhodphi(std::vector<double>& poten_array, std::vector<double>& dtworho_array, std::vector<double>& density_array, unsigned steps)
{
  // double dvalue = 0.0
  dtworho_array[0] = 0;
  //dtworho_array[steps] = 0;
  for (unsigned i =2; i<(steps-1); ++i)
  {
    
    dtworho_array[i] = (density_array[(i+1)] +  density_array[(i-1)] - 2*density_array[i])/((poten_array[(i-1)] - poten_array[i]) * (poten_array[i] -  poten_array[(i+1)]));
    //dtworho_array[i] = ((density_array[(i+1)] - density_array[i])  - (density_array[i] - density_array[(i-1)]))/((poten_array[(i-1)] - poten_array[i]) * (poten_array[i] -  poten_array[(i+1)]));
  } 
  dtworho_array[steps-1]=dtworho_array[steps-2];
  //std::cout << dtworho_array[steps] << "    "  << dtworho_array[steps-1] << std::endl;
}


 
double gausscheb(std::vector<double>& poten_array, std::vector<double>& dtworho_array, std::vector<double>& contdtworho_array, unsigned steps, double maxen)
{
  double gaussum = 0.0;
  // xi = gaussian cehb x
  // ti, rescaled variable so goes to -1 to 1
  // fxi is function(x) in gaussian
  double xi = 0.0;
  double fxi = 0.0;
  double ti = 0.0;
  steps = steps/10;
  double dtworhoxi = 0.0;

  // weighting for gauss
  double wxi = pi / steps;

  double gaussi = 0.0;
  for (unsigned i = 1; i<(steps+1);++i)
    {
      //calc xi at i
      xi = cos((2.0*i-1.0)*pi/(2*steps));
      //xi = cos((2.0i-1.0)*pi/(4*steps));
      ti = (maxen*pow(xi,1.0)+maxen)/2.0;
      //ti = maxen*xi;
      //spline(poten_array.data(), dtworho_array.data(), (poten_array.size()), dyn1, dyn2,contdtworho_array.data());
      dtworhoxi = splint(poten_array.data(), dtworho_array.data(), contdtworho_array.data(), poten_array.size(), fabs(ti));
      fxi = fabs(dtworhoxi)*(pow((maxen/2),0.5)) *pow((1+xi),0.5);
      // gauss = /sum wi * fxi
      gaussi =  fxi;
      gaussum = gaussum + gaussi;
      //std::cout << "Plint fxi = "<< fxi << "  xi = " << xi << "  ti  = "  << ti << std::endl;
    }
  //wxi is same for all values of i
  gaussum = gaussum*wxi;
  //std::cout << "sum = "<< gaussum << "  wi = " << wxi << std::endl;
  return gaussum;
}

void fedding(std::vector<double> poten_array, std::vector<double> dtworho_array, std::vector<double> contdtworho_array,std::vector<double>& fedd_array, unsigned steps, double minen)
{
  
  double fe = 0.0;
  double maxen;
  
  for (unsigned i = 0; i<(steps);++i)
    {
      //std::cout << "test"<< i << std::endl;
      //Pass through e value
      maxen = poten_array[i];
      //calc eddington formula
      fe = 1/(pow(8.0,0.5)*pow(pi,2.0))* (gausscheb(poten_array, dtworho_array, contdtworho_array, steps, maxen) + minen);
      fedd_array[i] = fe;
      //std::cout << "Maxen = "<< maxen << "  dtworho at maxen = " << (splint(poten_array.data(), dtworho_array.data(), contdtworho_array.data(), poten_array.size(),maxen)) << std::endl;
    }
  fedd_array[steps-1] = fedd_array[steps-2];
  //std::cout << "steps and steps -1 = " <<  fedd_array[steps] << "   "<< fedd_array[steps-1] << std::endl; 
}


//void rectant()
//{
  //  double rect = 0.0;
  //}


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
  cout << "Jean's equation simulation. \n Please input plummer radius (unitless)." << std::endl;
  cin >> pluma;
  //pluma = pluma*3.086*pow(10.0,16.0);
  //convert to metres
  
  cout << "Please inner integration radius (pc) (smallest radius for program) " << std::endl;
  cin >> innerr;
  //innerr = innerr*3.086*pow(10.0,16.0);
  //convert to metres
  
  cout << "Please input out integration radius. (pc) (similar to radius of galaxy)" << std::endl;
  cin >> outerr;
  //outerr = outerr*3.086*pow(10.0,16.0);
  //convert to metres
  
  cout << "Please input no. steps" << std::endl;
  cin >> steps;

  cout << "Please input Mass of Galaxy (solar masses)" << std::endl;
  cin >> galmass;
  //galmass = galmass * 1.989 * pow(10.0,30.0);
  //convert into kg's
    
  cout << "Please input additional black hole mass (0 for none)(solar masses)" << std::endl;
  cin >> blackmass;
  //blackmass = blackmass* 1.989 * pow(10.0,30.0);
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
  double phin = 0;
  
  // need to calc M(<r) from \rho
  double h = (outerr - innerr)/steps;    
  //makeradius(rad_array, steps, innerr, outerr);
  makelogradius(rad_array, steps, innerr, outerr);
  makerho(density_array,rad_array, steps, galmass, pluma);
  massr(density_array, mass_array, rad_array, steps, h);
  phin = potental(poten_array, mass_array, rad_array, steps, h);
  




  drhodphi(poten_array, drho_array, density_array, steps);
  //dtest(drho_array, density_array, rad_array,  steps, pluma, galmass);
  //unit testing dtwo function
  std::vector<double> sin_array(steps);

  for (unsigned i = 0; i<steps; ++i)
  {
    sin_array[i] = sin (rad_array[i]);
  }

  //trying first der of first der
  drhodphi(poten_array, dtworho_array, drho_array, steps);

  //dtworhodphi(rad_array, dtworho_array, sin_array, steps);
  
  //dtworhodphi(poten_array, dtworho_array, density_array, steps);
  
  double dyn1 = 0;
  double dyn2 = 1.6;
  std::vector<double> contdtworho_array(steps);

  spline(poten_array.data(), dtworho_array.data(), (poten_array.size()), dyn1, dyn2,contdtworho_array.data());

  //eddington function
  std::vector<double> fedd_array(steps);
  double minen = poten_array[3];
  std::cout << minen << std::endl;
  fedding(poten_array, dtworho_array, contdtworho_array, fedd_array, steps,minen);



  std::ofstream myfile;
  myfile.open (filename.c_str());
  //double h = (outerr - innerr)/steps;
  //r = 0;
  double poten = 0.0;
  double massen = 0.0;
  double potennewt = 0.0;
  double firstder = 0.0;
  double secder = 0.0;
  double eddan = 0.0;
  double lmass = 1.0;
  for (unsigned i =0; i<steps; ++i)
  {
    //analy potential
    poten =  Gconst * galmass/pluma * pow((1.0 + pow(rad_array[i]/pluma,2.0)),-1.0/2.0) -  Gconst * galmass/pluma * pow((1.0 + pow(rad_array[0]/pluma,2.0)),-1.0/2.0);
    //analy mass
    massen =  galmass * pow((1.0 + pow(rad_array[i]/pluma,-2.0)),-3.0/2.0);
    //potential newtonian assuming all mass in centre... debug test for michael
    potennewt = Gconst * mass_array[i] / pow(rad_array[i],1.0);
    // analytic first derivative
    firstder = 15.0*pow(pluma,2.0)*pow((poten_array[i]),4)/ (4.0*pi*pow(galmass,4.0)*pow(Gconst,5.0));
    // analytic second derivative 
    secder = 15.0*pow(pluma,2.0)*pow((poten_array[i]),3)/ (pi*pow(galmass,4.0)*pow(Gconst,5.0));
    
    //eddington
    
    eddan = 3.0*pow(2.0, 7.0/2.0)* pow(pluma, 2.0) / (7.0*pow(pi,3.0) *pow(Gconst, 5.0) * pow(galmass,4.0) * lmass) * pow(( poten),7.0/2.0); 

    // output file

    myfile << rad_array[i] << "   " << density_array[i] << "   "<< mass_array[i]  <<"  " << poten_array[i] << "   " << drho_array[i] << "   " << dtworho_array[i] << "    "<< poten << "   " << massen << "   " << fedd_array[i] << "    " << firstder << "    " << secder << "    " << eddan  << "\n";
  }

  myfile.close();

  return 0;

}
