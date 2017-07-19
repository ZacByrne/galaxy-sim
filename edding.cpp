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
  //std::cout << khi << "   klo " << klo << std::endl


  h=xa[khi]-xa[klo];
 if (h == 0.0) std::cout<<"Bad xa input to routine splint\n";

 a=(xa[khi]-x)/h;
 b=(x-xa[klo])/h;
 return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void makelogradius(std::vector<double>& rad_array, unsigned steps, double innerr, double outerr)
{
  // make an array of log spaced radius values                                                                                                                                                              
  double loginr = log10(innerr);
  double logoutr = log10(outerr);
  double h = (logoutr - loginr)/steps;
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
    {

      radius = loginr + i*h;

      rad_array[i]= pow(10.0, radius);
    }
}

double rhor(double radius,double galmass, double pluma)
{
  //give plummer density given radius, galaxy mass and plummer radius                                                                                                                                       
  double rhor;
  rhor = 3.0*galmass / (4.0*pi*pow(pluma,3.0)) * pow((1.0 + pow(radius/pluma,2)),-5.0/2.0);
  return rhor;
}

void makerho(std::vector<double>& density_array,std::vector<double>& rad_array, unsigned steps, double galmass, double pluma)
{
  // takes a radius array and creates density array. Currently using analytic plummer model                                                                                                                 
  double radius = 0;
  for (unsigned i =0; i<steps; ++i)
    {
      radius = rad_array[i];
      density_array[i] = rhor(radius, galmass, pluma);
      //density_array[i]=  3.0*galmass / (4*pi*pow(pluma,3.0)) * pow((1.0 + pow(radius/pluma,2)),-5.0/2.0);
    }
}


void massr(std::vector<double>& mass_array, std::vector<double>& rad_array, unsigned steps, double galmass, double pluma)
{
  //double radius = 0;                                                                                                                                                                                      
  double massr= 0.0;
  double radius, rgap;
  mass_array[0] = 4.0*pi*rhor(rad_array[0],galmass, pluma) * pow(rad_array[0],3.0)/3.0;
  
  for (unsigned i =1; i<steps; ++i)
    {
      //maxradius = rad_array[i]; 
      massr = 0.0;
      rgap = (rad_array[i]-rad_array[i-1])/1000.0;
      //std::cout << "rad i = "<< rad_array[i] << "  radi-1 = " << rad_array[i-1] << "   rgap  " << rgap << std::endl;
      for (unsigned j = 0; j < 1000; ++j)
	{
	  // r = r(i-1) + dr 
	  radius = rad_array[i-1] + (j+0.5) * rgap;
	  //std::cout <<"radius = " <<  radius << std::endl;
	  // mass is sum of rho * r^2 * dr
	  massr = massr + 4.0*pi*rhor(radius,galmass,pluma)*pow(radius, 2.0)*rgap;
	}
      //std::cout << "Massr  =  " << massr << " mass_array[i-1]  =  " << mass_array[i-1] << std::endl;
      mass_array[i] = mass_array[i-1]+massr;
      // poten_array[i] = Gconst*mass_array[i]/rad_array[i]; 
    }
}


void oldpotential(std::vector<double>& poten_array, std::vector<double> mass_array, std::vector<double> rad_array, unsigned steps, double galmass, double pluma)
{
  //makes a potential array from mass array and rad array                                                                                                                                                   
  double potenr = -1.0*Gconst* mass_array[0]/rad_array[0];
  for (unsigned i =1; i<steps; ++i)
    {
      potenr = 0.0;
      for (unsigned j =1; j<i+1; ++j)
	{
	  // potential is sum over m(r)/r^2 * dr                                                                                                                                                              
	  potenr = potenr + (rad_array[j]-rad_array[(j-1)])*mass_array[(j)] / pow(rad_array[(j)],2.0);
	  //std::cout << "Potenr  = " << potenr  << "  mass r = " << mass_array[(steps-j)] << "  rad_array  = "  << rad_array[(steps-j)] << std::endl;                                                        
	}
      //                                                                                                                                                                                                      
      potenr = potenr + mass_array[(steps-1)]/rad_array[(steps-1)];
      potenr = -1.0*Gconst * potenr;
      //std::cout << "Potenr  = " << potenr  << "   phin = " << phin << std::endl;                                                                                                                            
      poten_array[i] = potenr;
      //poten_array[i] =                                                                                                                                                                                      
      //poten_array[i] = Gconst * mass_array[i] /(rad_array[i]) + phin;                                                                                                                                       
    }

  // Can add -1 at 0                                                                                                                                                                                        
  double phin = poten_array[steps-1];
  std::cout << "phin = " << phin << std::endl;
  for (unsigned i =0; i<steps; ++i)
    {
      poten_array[i] = poten_array[i]-phin;
    }

}

void potential(std::vector<double>& poten_array, std::vector<double> mass_array, std::vector<double> rad_array, unsigned steps, double galmass, double pluma)
{
  double potenr, radius, rgap;
  
  potenr = 0.0;
  poten_array[0] = -1.0*Gconst* mass_array[0]/rad_array[0];
  for (int i = (steps-1); i >= 0; --i)
    {
      //potenr = 0.0;
      rgap = (rad_array[i]-rad_array[i-1])/1000.0;
      for (int j = 1000; j >= 0; --j)
        {
	  
          // r = r(i-1) + dr 
          radius = rad_array[i-1] + (j+0.5) * rgap;
	  
	  //potenial = 4piG/r * sum(rho(r')*r'^2*dr)
	  potenr = potenr + rhor(radius, galmass, pluma)*pow(radius,1.0)*rgap;
	}
      // poten_array[i-1] = 
      poten_array[i] = -1.0*Gconst * mass_array[i]/rad_array[i] + -4.0*pi*Gconst * potenr;
    }

  double phin = poten_array[steps-1];
  // std::cout << "phin = " << phin << std::endl;
  for (unsigned i =0; i<steps; ++i)
    {
      poten_array[i] = poten_array[i]-phin;
    }
  //return phin;
}

void olddrhodphi(std::vector<double>& poten_array, std::vector<double>& drho_array, std::vector<double>& density_array, unsigned steps, double galmass)
{
  double denom;
  unsigned derbreak = 0;
  //analytic for first. Should change to something else                                                                                                                                                     
  drho_array[0] = 15.0*pow(1.0,2.0)*pow((poten_array[0]),4)/ (4.0*pi*pow(galmass,4.0)*pow(Gconst,5.0));;
  //drho_array[steps] = 0;                                                                                                                                                                                  

  for (int i =1; i<(steps); ++i)
    {
      derbreak = 0;
      int j = 1;

      //loop until solution is found                                                                                                                                                                          
      while (derbreak == 0){
	// if in bounds, check denominator                                                                                                                                                                    
	if ( ((i+j)<=(steps-1)) && ((i-j)>=0)  ){
	  denom = poten_array[(i-j)] - poten_array[(i+j)];
	  //std::cout << denom << std::endl;                                                                                                                                                                  
	  if ( fabs(denom) < (pow(10.0,-10.0))){
	    // if demon is too small, expand delta x range                                                                                                                                                    
            j = j+1;
            //std::cout << "I  = " << i << "  j = "<< j << std::endl;                                                                                                                                       
	  }else {
	    // if demon is ok size, input into dhro array and break loop                                                                                                                                      
	    drho_array[i] = (density_array[i-j] - density_array[(i+j)])/(poten_array[(i-j)] - poten_array[(i+j)]);
	    //std::cout << drho_array[i] << std::endl;                                                                                                                                                        
	    derbreak = 1;
	  }
	}else{
	  drho_array[i] = drho_array[i-1];
	  derbreak = 1;
	  //if all else fails, use previous answer                                                                                                                                                            
	}
      }
  }
  //drho_array[steps-1] = drho_array[steps-3];
  //drho_array[steps-2] = drho_array[steps-3];
}

void drhodphi(std::vector<double>& poten_array, std::vector<double>& drho_array, std::vector<double>& density_array, unsigned steps, double galmass)
{
  //drho_array[0] = 0.0;
  for (unsigned i = 1; i < (steps-1); ++i)
    {
      drho_array[i] = (density_array[i-1] -  density_array[(i+1)])/(poten_array[(i-1)] - poten_array[(i+1)]);
    }
  
  drho_array[steps-1] = 0.0;

  drho_array[0] = drho_array[1];
  //drho_array[steps-2] = drho_array[steps-3];
  
}


void splined(std::vector<double> poten_array, std::vector<double>& drho_array,std::vector<double> rad_array, std::vector<double> density_array, unsigned steps, double galmass, double dyn1, double dyn2)
{
  // Sam attempt 2, spline den - pot, solve drho with equidistant potential, spline answer with equidistant potential, recall with desired potential                                                        
  std::vector<double> midspline_array(steps);
  std::vector<double> middrho_array(steps);
  spline(poten_array.data(), density_array.data(), (poten_array.size()), dyn1, dyn2, midspline_array.data());

  std::vector<double> midpot_array(steps);
  double basepoten = poten_array[0];
  double delp = (poten_array[steps-1] - poten_array[0])/steps;
  double den1, den2, pon1, pon2;
  midpot_array[0] = basepoten;

  for (unsigned i = 1; i < (steps-1); ++i)
    {
      pon1 = basepoten + delp*(i-1);
      pon2 = basepoten + delp*(i+1);
      den1 = splint(poten_array.data(), density_array.data(), midspline_array.data(), poten_array.size(), pon1);
      den2 = splint(poten_array.data(), density_array.data(), midspline_array.data(), poten_array.size(), pon2);
      middrho_array[i] = (den2 - den1)/(pon1 - pon2);
      midpot_array[i] = basepoten + delp*i;

      //std::cout << " den 1 = " << den1 << "  den 2  = " << den2 << "   pot1 = " << pon1 << "   pot 2 = " << pon2 << std::endl;

      //std::cout << " poten = " <<  midpot_array[i] << "    den = " << middrho_array[i] << std::endl;
    }
  middrho_array[0] = middrho_array[1];

  middrho_array[steps-2]= middrho_array[steps-3];
  middrho_array[steps-1] = middrho_array[steps-2];
  midpot_array[steps-1] = poten_array[steps-1];

  std::vector<double> dtwomid_array(steps);

  dyn1 = (middrho_array[1]-middrho_array[0])/delp;
  dyn2 = (middrho_array[steps-1]-middrho_array[steps-2])/delp;

  spline(midpot_array.data(), middrho_array.data(), (poten_array.size()), dyn1, dyn2, dtwomid_array.data());

  for (unsigned i = 0; i < steps; ++i)
    {
      drho_array[i] = -1.0* splint(midpot_array.data(), middrho_array.data(), dtwomid_array.data(), poten_array.size(), poten_array[i]);
    }
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
      fe = 1/(pow(8.0,0.5)*pow(pi,2.0))* (gausscheb(poten_array, dtworho_array, contdtworho_array, steps, maxen)); // adding in density  = 0 point?  + minen);                     
      fedd_array[i] = fe;
      //std::cout << "Maxen = "<< maxen << "  dtworho at maxen = " << (splint(poten_array.data(), dtworho_array.data(), contdtworho_array.data(), poten_array.size(),maxen)) << std::endl; 
    }
  fedd_array[steps-1] = fedd_array[steps-2];
  //std::cout << "steps and steps -1 = " <<  fedd_array[steps] << "   "<< fedd_array[steps-1] << std::endl;                                                                       
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
  // double phin = 0;

  // need to calc M(<r) from \rho                                                                                                                                                                        
  //double h = (outerr - innerr)/steps;
  //makeradius(rad_array, steps, innerr, outerr);                                                                                                                                                   
  makelogradius(rad_array, steps, innerr, outerr);
  makerho(density_array,rad_array, steps, galmass, pluma);
  massr(mass_array, rad_array, steps, galmass, pluma);
  potential(poten_array, mass_array, rad_array, steps, galmass, pluma);
  
  //oldpotential(poten_array, mass_array, rad_array, steps, galmass, pluma);


  double dyn1 = 1.189;
  double dyn2 = 0;
  for (unsigned i =0; i<steps; ++i)
    {
      poten_array[i] = poten_array[i];
      // std::cout << poten_array[i] << std::endl;                                                                                                                                     
    }



  //drhodphi(poten_array, drho_array, density_array, steps, galmass);

  splined(poten_array, drho_array,rad_array, density_array,steps, galmass, dyn1, dyn2);

  dyn1 = 4.776;
  dyn2 = 0;

  drhodphi(poten_array, dtworho_array, drho_array, steps, galmass);

  for (unsigned i =0; i<steps; ++i)
    {
      drho_array[i] = -1.0*drho_array[i];
      // std::cout << poten_array[i] << std::endl;                                                                                                                                         
    }


  //splined(poten_array, dtworho_array,rad_array, drho_array,steps, galmass, dyn1, dyn2);

  
  for (unsigned i =0; i<steps; ++i)
    {
      dtworho_array[i] = -1.0*dtworho_array[i];
      // std::cout << poten_array[i] << std::endl;                                                                                                                                                          
    }


  std::vector<double> contdtworho_array(steps);


  spline(poten_array.data(), dtworho_array.data(), (poten_array.size()), dyn1, dyn2,contdtworho_array.data());

  //eddington function                                                                                                                                                      
  std::vector<double> fedd_array(steps);
  double minen = poten_array[3];
  //std::cout << minen << std::endl;
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
      poten =  Gconst * galmass/pluma * pow((1.0 + pow(rad_array[i]/pluma,2.0)),-1.0/2.0) -  Gconst * galmass/pluma * pow((1.0 + pow(rad_array[steps-1]/pluma,2.0)),-1.0/2.0);
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

      myfile << rad_array[i] << "   " << density_array[i] << "   "<< mass_array[i]  <<"  " << -1.0*poten_array[i] << "   " << drho_array[i] << "   " << dtworho_array[i] << "    "<< poten << "   " << massen << "   " << fedd_array[i] << "    " << firstder << "    " << secder << "    " << eddan  << "\n";
    }

  myfile.close();

  return 0;

}






