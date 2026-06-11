from numpy import log, pi, exp, interp, array, zeros_like, concatenate, linspace, vectorize
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
GeV3_to_fm3 = 1.0e3/1.9732705**3

################
# T = 0, high mu
################
# X = Lbar/(2*mu) !!!
nf = 3.
nc = 3.
da = 8.
delta = -0.8563832093269428068483102329159
zeta3 = 1.2020569031595942853997381615114499907649

ca = nc
cf = da/2./nc
tf = 1./2.

# from Schwartz QFT and https://arxiv.org/pdf/hep-ph/9701390.pdf
beta0 = 11./3.*ca-4./3.*tf*nf
beta1 = 34./3*ca**2. - 4.*cf*tf*nf - 20./3.*ca*tf*nf 
beta2 = 2857./54.*ca**3. + 2.*cf**2.*tf*nf - 205./9.*cf*ca*tf*nf - 1415./27.*ca**2*tf*nf + 44./9.*cf*tf**2.*nf**2. + 158./27.*ca*tf**2.*nf**2
beta3 = (149753./6.+3564.*zeta3) - (1078361./162. + 6508./27.*zeta3)*nf + (50065./162. + 6472./81.*zeta3)*nf**2. + (1093./729.)*nf**3.

c11 = -2.
#
c21 = -3.
def c22(X):
    return -3.*(3.*log(X)+5.0021)
#
c32 = 9.*11./12.
def c31(X):
    return 9.*(-6.5968 - 3*log(X))
def c30(X,c0):
    return 9.*(5.1342 + 2./3.*c0 - 18.284*log(X) - 9./2.*log(X)**2)

def p0h(X,c0):
    return -9./4.*log(X)**2 - 26.367*log(X) + c0
def c30sm(X,c0):
    return c30(X,c0) - 9.*2./3.*p0h(X,c0) #this cancels the c0 dependence, so it must be correct

def dc30dc0(X,c0):
    return 9.*2./3.

# beta function
# these are da_n/dlbar = beta_n*lbar
def beta0_func_X_lbar(lbar, a_s):
    return -2*a_s*((a_s/4/pi)*beta0)/lbar
def beta1_func_X_lbar(lbar, a_s):
    return -2*a_s*((a_s/4/pi)*beta0 + (a_s/4/pi)**2*beta1)/lbar
def beta2_func_X_lbar(lbar, a_s):
    return -2*a_s*((a_s/4/pi)*beta0 + (a_s/4/pi)**2*beta1 + (a_s/4/pi)**3*beta2)/lbar
def beta3_func_X_lbar(lbar, a_s):
    return -2*a_s*((a_s/4/pi)*beta0 + (a_s/4/pi)**2*beta1 + (a_s/4/pi)**3*beta2 + (a_s/4/pi)**4*beta3)/lbar


# pressure
def PLO(a_s):
    return 1.

def PNLO(a_s):
    return c11*(a_s/pi)

def PN2LO(a_s,X):
    return (a_s/pi)**2*(c21*log(nf*a_s/pi) + c22(X))

def PN3LO(a_s, X, c0):
    return (a_s/pi)**3*(c32*log(nf*a_s/pi)**2
                        + c31(X)*log(nf*a_s/pi) 
                        + c30(X, c0))

def dPN3LOdc0(a_s, X, c0):
    return (a_s/pi)**3*(dc30dc0(X, c0))


def PN3LO_sm(a_s, X, c0): 
    return (a_s/pi)**3*(c32*log(nf*a_s/pi)**2
                        + c31(X)*log(nf*a_s/pi) 
                        + c30sm(X, c0))

def dPNLO(a_s):
    return c11*(1/pi)

def dPN2LO(a_s,X):
    return 2*a_s/pi**2*(c21*(log(nf*a_s/pi)+1/2) + c22(X))

def dPN3LO(a_s,X,c0):
    return 3*a_s**2/pi**3*(
        c32*(log(nf*a_s/pi)**2+2/3*log(nf*a_s/pi)) +
        c31(X)*(log(nf*a_s/pi)+1/3) + 
        c30(X, c0))

# The running coupling, see Schwartz 2013
def alpha_s_NLO_pert(lbar):
    L = log((lbar)**2/(0.19424752576669174)**2)
    return 4*pi/beta0*(1/L)

def alpha_s_N2LO_pert(lbar):
    L = log((lbar)**2/(0.3778938798585379)**2)
    return 4*pi/beta0*(1/L - beta1/beta0**2/L**2*log(L))

def alpha_s_N3LO_pert(lbar):
    L = log((lbar)**2/(0.34408323550919756)**2)
    return 4*pi/beta0*(1/L - beta1/beta0**2/L**2*log(L) + beta1**2/beta0**4/L**3*(log(L)**2 - log(L) - 1) + beta2/beta0**3/L**3)

## ODE running
def alpha_s_NLO(lbar):
    return solve_ivp(beta0_func_X_lbar, [2, lbar], [0.2994]).y[0][-1]

def alpha_s_N2LO(lbar):
    return solve_ivp(beta1_func_X_lbar, [2, lbar], [0.2994]).y[0][-1]

def alpha_s_N3LO(lbar):
    return solve_ivp(beta2_func_X_lbar, [2, lbar], [0.2994]).y[0][-1]

def alpha_s_N4LO(lbar):
    return solve_ivp(beta3_func_X_lbar, [2, lbar], [0.2994]).y[0][-1]

#lbar_1_tab =  exp(linspace(log(0.02),log(10),1000))
#as_1_tab = vectorize(alpha_s_NLO_ODE)(lbar_1_tab)
#lbar_2_tab =  exp(linspace(log(0.02),log(10),1000))
#as_2_tab = vectorize(alpha_s_N2LO_ODE)(lbar_2_tab)
#lbar_3_tab =  exp(linspace(log(0.02),log(10),5000))
#as_3_tab = vectorize(alpha_s_N3LO_ODE)(lbar_3_tab)

#alpha_s_NLO  = interp1d(lbar_1_tab, as_1_tab)
#alpha_s_N2LO = interp1d(lbar_2_tab, as_2_tab)
#alpha_s_N3LO = interp1d(lbar_3_tab, as_3_tab)

# this is the full dalpha/dmu
def das_N2LO_dmu(mu,X):
	numerator = -2.20644 - 2.79253*log(0.777632*mu**2*(2*X)**2) + 4.41288*log(log(0.777632*mu**2*(2*X)**2))
	denominator = mu*(log(0.777632*mu**2*(2*X)**2))**3
	return numerator/denominator

def pFD(mu):
	return (mu)**4/(108*pi**2)

def dpFD(mu):
	return mu**3/(27*pi**2)

def d2pFD(mu):
	return mu**2/(9*pi**2)


# pressure    
def p(mu,X,order,c0): #Output in MeV/fm3
    lbar = 2*(mu/3)*X
    if order == 0:
        a_s = 1
        return (PLO(a_s))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 1: 
        a_s = alpha_s_N3LO(lbar) 
        return (PLO(a_s) + PNLO(a_s))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 2: 
        a_s = alpha_s_N3LO(lbar)
        return (PLO(a_s) + PNLO(a_s) + PN2LO(a_s, X))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 3:
        a_s = alpha_s_N3LO(lbar)
        return (PLO(a_s) + PNLO(a_s) + PN2LO(a_s, X) + PN3LO(a_s, X, c0))*pFD(mu)*GeV3_to_fm3*1.e3
    
# pressure    
def p_sm(mu,X,order,c0): #Output in MeV/fm3
    lbar = 2*(mu/3)*X
    if order == 0:
        a_s = 1
        return (PLO(a_s))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 1: 
        a_s = alpha_s_N3LO(lbar) 
        return (PLO(a_s) + PNLO(a_s))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 2: 
        a_s = alpha_s_N3LO(lbar)
        return (PLO(a_s) + PNLO(a_s) + PN2LO(a_s, X))*pFD(mu)*GeV3_to_fm3*1.e3
    if order == 3:
        a_s = alpha_s_N3LO(lbar)
        return (PLO(a_s) + PNLO(a_s) + PN2LO(a_s, X) + PN3LO_sm(a_s, X, c0))*pFD(mu)*GeV3_to_fm3*1.e3
    
def dpdc0(mu,X,order,c0): #Output in MeV/fm3
    lbar = 2*(mu/3)*X
    if order == 0:
        return 0.
    if order == 1: 
        return 0.
    if order == 2: 
        return 0.
    if order == 3:
        a_s = alpha_s_N3LO(lbar)
        return dPN3LOdc0(a_s, X, c0)*pFD(mu)*GeV3_to_fm3*1.e3


def n(mu,X,order,c0): #Output in 1/fm3
    lbar = 2*(mu/3)*X
    if order == 0:
        return dpFD(mu)*GeV3_to_fm3
    if order == 1: 
        a_s = alpha_s_N3LO(lbar) 
        das_dmu = (lbar/mu)*beta2_func_X_lbar(lbar, a_s)

        p_as = PLO(a_s) + PNLO(a_s)
        dp_das = dPNLO(a_s)
        return (dp_das*das_dmu*pFD(mu) + p_as*dpFD(mu))*GeV3_to_fm3
    if order == 2:
        a_s = alpha_s_N3LO(lbar)
        das_dmu = (lbar/mu)*beta2_func_X_lbar(lbar, a_s)

        p_as = PLO(a_s) + PNLO(a_s) + PN2LO(a_s,X)
        dp_das = dPNLO(a_s) + dPN2LO(a_s, X)
        return (dp_das*das_dmu*pFD(mu) + p_as*dpFD(mu))*GeV3_to_fm3
    if order == 3: 
        a_s = alpha_s_N3LO(lbar)
        das_dmu = (lbar/mu)*beta2_func_X_lbar(lbar, a_s)

        p_as = PLO(a_s) + PNLO(a_s) + PN2LO(a_s,X) + PN3LO(a_s, X, c0)
        dp_das = dPNLO(a_s) + dPN2LO(a_s, X) + dPN3LO(a_s, X, c0)
        return (dp_das*das_dmu*pFD(mu) + p_as*dpFD(mu))*GeV3_to_fm3
    
def n_pert(mu,X,order,c0):
    lbar = 2*(mu/3)*X
    if order == 0:
        return dpFD(mu)*asGeV3_to_fm3
    if order == 1: 
        a_s = alpha_s_N3LO(lbar) 
        p_as = PLO(a_s) + PNLO(a_s)
        dp_das = dPNLO(a_s)
        return (p_as*dpFD(mu))*GeV3_to_fm3
    if order == 2:
        a_s = alpha_s_N3LO(lbar)
        p_as = PLO(a_s) + PNLO(a_s) + PN2LO(a_s,X)
        dp_das_das_term = dPNLO(a_s)*beta0_func_X_lbar(lbar, a_s)
        return (p_as*dpFD(mu) + pFD(mu)*dp_das_das_term*(lbar/mu))*GeV3_to_fm3
    if order == 3: 
        a_s = alpha_s_N3LO(lbar)
        p_as = PLO(a_s) + PNLO(a_s) + PN2LO(a_s,X) + PN3LO(a_s, X, c0)
        dp_das_das_term = dPNLO(a_s)*beta1_func_X_lbar(lbar, a_s) + dPN2LO(a_s,X)*beta0_func_X_lbar(lbar, a_s)
        return (p_as*dpFD(mu) + pFD(mu)*dp_das_das_term*(lbar/mu))*GeV3_to_fm3