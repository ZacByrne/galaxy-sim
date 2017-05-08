#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include <cassert>
//#include <boost/lexical_cast.hpp>

//#include <stdio.h>
//#include <math.h>
//#define NRANSI
//#include "nr.h"
//#include "nrutil.h"

// NULL is gone, use nullptr
//double* me = nullptr;
//double pluma = 1;
double pi = 3.1415926535;
//GET VALUE FOR HOLGER

//double galmass = 1:
double Gconst = 6.674*pow(10.0,-11.0);

//function from Holger for cubic spline
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

//function from holger for splint (get value given x
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


//beta
//double beta(double v_r2, double v_theta2)
//{
 // double betavalue = 0;
  //  //1 - (v_theta2)/(v_r2);
//  return betavalue;
//}

//stellar density at r in plummer model
double vsteldelplum(double r, double galmass, double pluma)
{
  double rho = 3*galmass / (4*pi*pow(pluma,3.0)) * pow((1 + pow(r/pluma,2)),-5.0/2.0);
  return rho;
}

//first deriviative of stellar density at r
double dvstelplum(double r, double galmass, double pluma)
{
  double drho = -15.0/4.0 * galmass*r/(pi*pow(pluma,5.0)) * pow((1 + pow(r/pluma,2)),-7.0/2.0);
  return drho;
}

//Mass Function at r of plummer model
double Mfuncplum(double r, double galmass, double pluma)
{
  double massr = galmass * (1 + pow(pluma/r,2.0));
  return massr;
}

//derivative of gfravitational potential
double dgravpotplum(double r, double galmass, double pluma)
{
  double dphi = Gconst * galmass * r * pow(pluma, -3.0) * pow((1 + pow(r/pluma,2)),-3.0/2.0);
  return dphi;
}

//function of r and yn for midpoint formula, calculates velocity dispersion first deriviatiuve
double fry(double r, double y_n, double galmass, double blackmass, double pluma, double beta)
{
  double fvalue = -1*dgravpotplum(r,galmass, pluma) - Gconst*blackmass/pow(r,2.0) - 2*beta* y_n / r - dvstelplum(r, galmass, pluma)*y_n/vsteldelplum(r, galmass, pluma);
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


//midpoint ODE integration function with 
void midpointarray(std::vector<double>& v_rarray, double start, double end, unsigned step, double galmass, double blackmass, double innerr, double beta)
{
  double s =0;
  double h = (innerr - start)/step;
  double yn = s;
  double rn = 0;
  for (unsigned i =0; i<step; ++i)
   {
      yn = s;
      v_rarray[i]=yn;
      rn = start + i*h;
      s = s + h/2*fry(rn,yn, galmass, blackmass, end, beta);
      s = yn + h*fry(rn+h/2, s, galmass, blackmass, end, beta);
      
      // changed to not recursive 
        //s = s + h/2*fry((rn+h/2), (yn + h/2*fry(rn,yn)));
      
   }
}

// make array of radius values for splint and spline
void radiusarray(std::vector<double>& radarray,double start, double end, unsigned step)
{
  double h = (start - end) / step;
  for (unsigned i =0; i<step; ++i)
  {
    radarray[i] = end + i*h;
    
  }
}

// get projection of 3d along a given line of sight. GIven r, integrate over z
double projection(std::vector<double> v_rarray,std::vector<double> radarray,std::vector<double> vrdarray, double start, double end, unsigned step, double galrad,double galmass, double pluma)
{
  double intertop = 0;
  double interbot = 0;
  //start = start - 20000.0;
  double hypt = 0;
  double vstel = 0;
  double popden = 0;
  double z = 0;
  double h = (start - end)/step;
  //std::cout << end << "  start " << start << std::endl;
  //std::cout << h << " step =  " << step << std::endl;
  for (unsigned i = 0; i<step; ++i)
  {
    z = end + i*h;
    //std::cout<<" h = " << h <<" i =  " <<i <<" start " << start << std::endl;
    //std::cout <<i << std::endl;
    hypt = pow((pow(galrad,2.0)+pow(z,2.0)), 0.5);
    //std::cout << "test" << std::endl;
    //    std::cout << hypt << "   " << start << std::endl;
    //std::cout << "test";
    //if (hypt > (start+200000))
    //{
    //break;
    // }
    //std::cout << " z = " << z << "  start = " << (start+2000.0) << std::endl;
    vstel = splint(radarray.data(), v_rarray.data(), vrdarray.data(), radarray.size(), hypt);
    popden = vsteldelplum(hypt, galmass, pluma);
    intertop += (vstel * popden);
    interbot += popden;
    //std::cout<< "vstel  = " << vstel <<"  intertop  = "<< intertop   <<std::endl;
  }
  double project = intertop/interbot;
  return project;

    
    
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
  double innerr;
  double beta;
  
  using std::cout;
  using std::cin;
  
  //input values
  cout << "Jean's equation simulation. \n Please input core radius (Kpc)." << std::endl;
  cin >> pluma;
  pluma = pluma*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please input end radius. (Kpc)" << std::endl;
  cin >> innerr;
  innerr = innerr*3.086*pow(10.0,19.0);
  
  cout << "Please input start radius. (Kpc)" << std::endl;
  cin >> start;
  start = start*3.086*pow(10.0,19.0);
  //convert to metres
  
  cout << "Please input no. steps" << std::endl;
  cin >> steps;

  cout << "Please input Mass of Galaxy (solar masses)" << std::endl;
  cin >> galmass;
  galmass = galmass * 1.989 * pow(10.0,30.0);
    
  cout << "Please input additional black hole mass (0 for none)(solar masses)" << std::endl;
  cin >> blackmass;
  blackmass = blackmass* 1.989 * pow(10.0,30.0);
  
  cout << "Please input Beta (anisotropy value)" << std::endl;
  cin >> beta;
    
  
  std::string filename;  
  cout << "Please input output filename" << std::endl;
  cin >> filename;
  
  //Check with holger about end of int being radius
  //end = pluma;
  
  //vector array
  std::vector<double> vrarray(steps);
  std::vector<double> radarray(steps);
  std::vector<double> vsortarray(steps);

  // Run midpoint ODE formula
  midpointarray(vrarray, start, pluma, steps, galmass, blackmass,innerr,beta);
  
  // Create sorted radius array for spline and splint
  radiusarray(radarray, start, innerr, steps);
  //cout << radarray[4];

  //resort v array for splint and spline (SHOULD REWRITE MIDPOINT CODE>>>>
  for (unsigned i =0; i<steps; ++i)
  {
    vsortarray[i] = vrarray[(vrarray.size()-i)];
  }

  //Spline function and deriviatives
  double dyn1 = 0;
  double dyn2 = 0;
  std::vector<double> vrdarray(steps);

  spline(radarray.data(), vsortarray.data(), (vsortarray.size()), dyn1, dyn2,vrdarray.data());
  // cout << vrdarray[2];

  cout << "rad array done" <<std::endl;


  //Projection array
  std::vector<double> projectarray(steps);
  double h = (start-innerr)/steps;
  double r = 0;
  for (unsigned i = 0; i < steps; ++i)
  {
    r = innerr + i*h;
    //cout << r << std::endl; 
    projectarray[i] = projection(vsortarray, radarray, vrdarray,start,innerr, steps, r,galmass, pluma);
    //cout << projectarray[i] << std::endl;
  }

  cout << "Projetion done" << std::endl;

  std::ofstream myfile;
  myfile.open (filename.c_str());
  h = (innerr-start)/steps;
  r = 0;
  double sigma2 = 0;
  double sigma3 = 0;
  //cout << "inerrer = " << innerr <<"  start  = " << start << "  h  =  "<< h <<std::endl;
  double testvd=0;
  for (unsigned i =0; i<steps; ++i)
  {
    r = start + i*h;
    sigma2 = Gconst*galmass / (6*pluma) * pow((1 + pow(r/pluma,2)),-1.0/2.0);
    sigma3 = pi*sigma2/4.0;
    testvd = splint(radarray.data(), vsortarray.data(), vrdarray.data(), radarray.size(), r);
    //    cout <<"vstel  = "<< vsteldelplum(r,galmass, pluma) << std::endl;
    myfile << r << "   " << vrarray[i] << "   "<< testvd  <<"  " << projectarray[(projectarray.size()-i)] << "   " << sigma3 << "\n";
  }
  //have changed the sigma2 to testvd for tasting purposes in output
  myfile.close();
  //  cout << "r = " << r << std::endl;
  return 0;

}



