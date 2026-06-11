// pi
const double pi=3.1415926535;
// GeV to fm^-3
const double gevtoinvfm=5.067730756672362;

// pQCD parameters (terms from 1807.04120 and 2103.05658, using notation of 2307.08734])
const double a11 = -2.;
const double a21 = -1.;
const double a22 = -3.;
const double a23 = -5.0021;
const double a3 = 1.67026; // this is the alpha_s^3 coefficient after using the PMS scale in 2103.05658
// alpha_s parameters from Schwartz QFT and https://arxiv.org/pdf/hep-ph/9701390.pdf
const double nf = 3.;
const double nc = 3.;
const double da = 8.;
const double ca = nc;
const double cf = da/2./nc;
const double tf = 1./2.;
const double beta0 = 11./3.*ca-4./3.*tf*nf;
const double beta1 = 34./3*pow(ca,2.) - 4.*cf*tf*nf - 20./3.*ca*tf*nf;
const double beta2 = 2857./54.*pow(ca,3.) + 2.*pow(cf,2.)*tf*nf - 205./9.*cf*ca*tf*nf - 1415./27.*pow(ca,2.)*tf*nf + 44./9.*cf*pow(tf,2.)*pow(nf,2.) + 158./27.*ca*pow(tf,2.)*pow(nf,2.);
const double LMSbar = 0.34408323550919756; // in GeV. This fixes alpha_s(2GeV) = 0.2994; see 0912.1856 but note that we are running at a higher order so the value of the scale is different
// range of renormalization scale parameter X, in new normalization
const double Xmin=0.5;
const double Xmax=2.0;