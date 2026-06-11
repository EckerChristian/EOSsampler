import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pQCD_N3LO import *

# Formatting 
from matplotlib.ticker import NullFormatter
import matplotlib as mpl
import matplotlib.pyplot as plt

###################################
# setup default plotting parameters
###################################

mpl.rcParams.update({
    'mathtext.fontset'            : 'cm',
    'axes.unicode_minus'          : False,
    'axes.formatter.use_mathtext' :  True,
})
plt.rcParams.update({
    'font.size': 14,
    'figure.figsize'    : [8.0, 5.0],
    #'text.usetex'       : True, 
    'font.family'       : "serif",
    'font.serif'        : "cmr10",
    'xtick.major.size'  : 6,
    'xtick.minor.size'  : 3,
    'ytick.major.size'  : 6,
    'ytick.minor.size'  : 3,
})
LabelSize=14
MarkerSize=75
MarkerSize=95
TickSize=20
LegendSize=14
nullfmt=NullFormatter() # no labels

#plt.rc('font', **font)

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
BIGGEST_SIZE = 24

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

###################################
# Data and plots
###################################

############
## pQCD 
############

cs2 = pd.read_csv("../build/cs2_pQCD.dat", sep=" ")
print(cs2)

x_vals = np.sort(list(set(cs2['X'])))
print(x_vals)
mu_vals = np.sort(list(set(cs2['muB[GeV]'])))
print(mu_vals)

for X in x_vals:
    datum = cs2[cs2["X"] == X]
    plt.plot(datum['muB[GeV]'], datum['cs2'])
plt.ylim([0.29, 0.53])
plt.xlabel(r'$\mu_B$ [GeV]')
plt.ylabel(r'$c_s^2$')

plt.tight_layout()
plt.savefig("../plot/pQCD_cs2_muB.pdf")
plt.close()

cs2_aves = []
cs2_stds = []
for mu in mu_vals:
    datum = cs2[cs2["muB[GeV]"] == mu]
    cs2_aves.append(np.average(datum['cs2']))
    cs2_stds.append(np.sqrt(np.var(datum['cs2'])))

np.savetxt("../pQCD/cs2_mu_mean_std_scale_ave.dat",np.column_stack((mu_vals, cs2_aves, cs2_stds)), header="muB[GeV] ave(cs2) std(cs2)", delimiter=" ", comments="")

for X in x_vals:
    datum = cs2[cs2["X"] == X]
    plt.plot(datum['muB[GeV]'], datum['cs2'], lw='0.1')

plt.plot(datum['muB[GeV]'], cs2_aves, c='black', label=r'$c_s^2$ scale ave.')
plt.plot(datum['muB[GeV]'], np.array(cs2_aves)+2.*np.array(cs2_stds), c='black', ls='--', label = r'$\pm 2\sigma$')
plt.plot(datum['muB[GeV]'], np.array(cs2_aves)-2.*np.array(cs2_stds), c='black', ls='--')
plt.ylim([0.19, 0.43])
plt.xlabel(r'$\mu_B$ [GeV]')
plt.ylabel(r'$c_s^2$')
plt.legend()

plt.tight_layout()
plt.savefig("../plot/pQCD_cs2_muB.pdf")
plt.close()

#for X in x_vals:
#    datum = cs2[cs2["X"] == X]
#    plt.plot(datum['n[fm^-3]']/0.16, datum['cs2'])
#plt.ylim([0.29, 0.53])
#plt.xlabel(r'$n_B [n_\mathrm{sat}]$')
#plt.ylabel(r'$c_s^2$')
#
#plt.tight_layout()
#plt.savefig("../plot/pQCD_cs2_nB.pdf")
#plt.close()

##for X in x_vals:
##    datum = cs2[cs2["X"] == X]
##    plt.plot(datum['muB[GeV]'], datum['p[MeV.fm^-3]']/datum['muB[GeV]']**4.)
###plt.ylim([0, 10])
##plt.xlabel(r'$\mu_B$ [GeV]')
##plt.ylabel(r'$p/\mu_B^4$')
##plt.tight_layout()
##plt.savefig("../plot/pQCD_pomu4_muB.pdf")
##plt.close()
##
##for X in x_vals:
##    datum = cs2[cs2["X"] == X]
##    plt.plot(datum['muB[GeV]'], datum['n[fm^-3]']/datum['muB[GeV]']**3./GeV3_to_fm3)
###plt.ylim([0, 10])
##plt.xlabel(r'$\mu_B$ [GeV]')
##plt.ylabel(r'$n_B/\mu_B^3$')
##plt.tight_layout()
##plt.savefig("../plot/pQCD_nBomu3_muB.pdf")
##plt.close()

############
## alpha_s
############

#alphas = pd.read_csv("../build/alphas_N3LO.dat", sep=" ")
#print(alphas)
#
#plt.plot(alphas['lbar'], alphas['alphas'])
#plt.xlabel(r'$\bar\Lambda$ [GeV]')
#plt.ylabel(r'$\alpha_s(\bar\Lambda)$')
#plt.title('Running of the coupling')
#plt.ylim([0.2,0.7])
#plt.tight_layout()
#plt.savefig("../plot/alphas_N3LO.pdf")
#plt.close()
#
#plt.plot(alphas['lbar'], np.gradient(alphas['alphas'], alphas['lbar']), label = 'numerical')
#plt.plot(alphas['lbar'], alphas['dasdlbar'], label = r'$\beta$ function')
#plt.xlabel(r'$\bar\Lambda$ [GeV]')
#plt.ylabel(r'$\alpha_s(\bar\Lambda)$')
#plt.title('Derivative of running coupling')
#plt.xlim([0.5,1.5])
#plt.legend()
#plt.tight_layout()
#plt.savefig("../plot/Dalphas_N3LO.pdf")
#plt.close()