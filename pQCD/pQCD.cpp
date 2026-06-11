#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <string>
#include "definitions.h"
using namespace std;

// function prototypes
double pressureQCD(double m, double X);
double baryondensityQCD(double m, double X);
double cs2QCD(double m,double X);

// pQCD running of alpha_s (N3LO)
double alpha_s(double lbar){
  double L = 2.*log(lbar/LMSbar);
  return 4.*pi/beta0*(1./L - beta1/pow(beta0,2.)/pow(L,2.)*log(L) + pow(beta1,2.)/pow(beta0,4.)/pow(L,3.)*(pow(log(L),2.) - log(L) - 1) + beta2/pow(beta0,3.)/pow(L,3.));
}
// Free pressure in GeV^4
double pfree(double mub){
  return pow(mub/3., 4.)*9./12./pi/pi;
}
// pQCD EOS as function of the chemical potential in MeV/fm^3
double pressureQCD(double mub, double X){
  double lbar = 2.*(mub/3.)*X;
  double as = alpha_s(lbar);

  return 1000.0*pow(gevtoinvfm, 3.)*pfree(mub) * (1. + (as/pi)*a11 + 3.*pow(as/pi, 2.)*(a21*log(3.*(as/pi)) + a22*log(X) + a23) + pow(3., 2.)*pow(as/pi, 3.)*a3);
}

// dpfree/dmub in GeV^3
double dpfree_dmub(double mub){
  return 4./3.*pow(mub/3., 3.)*9./12./pi/pi;
}  
double d2pfree_dmub2(double mub){
  return 4.*3./3./3.*pow(mub/3., 2.)*9./12./pi/pi;
}  
// dalpha_s/dmub and helper beta function and derivative of beta function. Here, das/dlbar = beta/lbar
double beta2_func_X_lbar(double lbar, double as){
  return -2.*as*((as/4./pi)*beta0 + pow(as/4./pi,2.)*beta1 + pow(as/4/pi,3.)*beta2)/lbar;
}
double Dbeta2_das_func_X_lbar(double lbar, double as){
  return -2.*(2.*(as/4./pi)*beta0 + 3.*pow(as/4./pi,2.)*beta1 + 4.*pow(as/4/pi,3.)*beta2)/lbar;
}
// pQCD baryon number density in 1/fm^3
double baryondensityQCD(double mub, double X){
  double lbar = 2.*(mub/3.)*X;
  double as = alpha_s(lbar);
  double das_dmub = (lbar/mub)*beta2_func_X_lbar(lbar, as);
  double p_as = (1. + (as/pi)*a11 + 3.*pow(as/pi, 2.)*(a21*log(3.*(as/pi)) + a22*log(X) + a23) + pow(3., 2.)*pow(as/pi, 3.)*a3);
  double dp_das = a11/pi + 3.*as/pi*(a21/2. + a21*log(3.*(as/pi)) + a22*log(X) + a23)*(2./pi) + pow(3.,2.)*pow(as/pi,2.)*a3*(3./pi);
  return pow(gevtoinvfm, 3.)*(dp_das*das_dmub*pfree(mub) + p_as*dpfree_dmub(mub));
}

// cs2 at pQCD matching point
double cs2QCD(double mub,double X){ 
  double lbar = 2.*(mub/3.)*X;
  double as = alpha_s(lbar);
  double das_dmub = (lbar/mub)*beta2_func_X_lbar(lbar, as);
  double d2as_dmub2 = -(lbar/mub/mub)*beta2_func_X_lbar(lbar, as) + (lbar/mub)*Dbeta2_das_func_X_lbar(lbar, as)*das_dmub;
  double p_as = (1. + (as/pi)*a11 + 3.*pow(as/pi, 2.)*(a21*log(3.*(as/pi)) + a22*log(X) + a23) + pow(3., 2.)*pow(as/pi, 3.)*a3);
  double dp_das = a11/pi + 3.*as/pi*(a21/2. + a23 + a21*log(3.*(as/pi)) + a22*log(X))*(2./pi) + pow(3.,2.)*pow(as/pi,2.)*a3*(3./pi);
  double d2p_das2 = 3.*(3.*a21/2. + a23 + a21*log(3.*(as/pi)) + a22*log(X))*(2./pi/pi) + pow(3.,2.)*(as/pi)*a3*(3.*2./pi/pi);
  double n = (dp_das*das_dmub*pfree(mub) + p_as*dpfree_dmub(mub));
  double dn_dmub = pfree(mub)*(d2p_das2*pow(das_dmub, 2.) + dp_das*d2as_dmub2) + 2.*dp_das*das_dmub*dpfree_dmub(mub) + p_as*d2pfree_dmub2(mub);
  
  return n / mub / dn_dmub;
  // units
}

// Speed of sound sampling program.
// compilation:
// g++ EOS.cpp -o EOS
// run:
// EOS 1 $RANDOM 
// creates outputfile: 1.eos
//
int main(int argc, char** argv){

    // opening the cs2 file for writing
    ofstream cs2_pQCD;
    cs2_pQCD.open("cs2_pQCD.dat", ios::out);
    cs2_pQCD<<"muB[GeV] X p[MeV.fm^-3] n[fm^-3] cs2"<<endl;

    // Output some pQCD data!
    int n_mu = 100;
    int n_X = 100;
    for(int cnt_mu = 0; cnt_mu < n_mu; cnt_mu++){
        double mumu = 2.25 + cnt_mu/(double (n_mu - 1))*(2.6 - 2.25);
        for(int cnt_X=0; cnt_X<n_X; cnt_X++){
            double XX=exp(log(Xmin)+(cnt_X/(double (n_X - 1)))*(log(Xmax)-log(Xmin)));
            cs2_pQCD<<std::scientific<<std::setprecision(6)<<mumu<<" "<<XX<<" "<<pressureQCD(mumu, XX)<<" "<<baryondensityQCD(mumu, XX)<<" "<<cs2QCD(mumu, XX)<<endl;
        }
    }
    cs2_pQCD.close();

    // opening the cs2 file for writing
    ofstream alphas;
    alphas.open("alphas_N3LO.dat", ios::out);
    alphas<<"lbar alphas dasdlbar"<<endl;

    //Output the running of alpha
    int n_alpha = 1000;
    for(int cnt_alpha = 0; cnt_alpha < n_alpha; cnt_alpha++){
        double lbar = 0.6 + cnt_alpha/(double (n_alpha - 1))*(5.-0.6);
        alphas<<std::scientific<<std::setprecision(6)<<lbar<<" "<<alpha_s(lbar)<<" "<<beta2_func_X_lbar(lbar, alpha_s(lbar))<<endl;
    }
    alphas.close();
 

    return 0;
}