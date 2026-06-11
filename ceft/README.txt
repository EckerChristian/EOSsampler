Correspond to Figure 1 in arXiv:2009.06441

The file names specify the chiral order (N2LO or N3LO) and the central ("cen" or mean) values and the plus ("pls) one-sigma boundary as well as the minus ("mns") one-sigma boundary of the EOS at this confidence level.

Since the underlying GP's correlation length is long, the results for "cen," "mns," and "pls" can be treated as three different random draws from the GP. We discuss this in more detail in the paper mentioned above.

The files are simple text files. The first line specifies the quantities and their units: pressure, energy density, and baryon density.

The file EoS_NSM.dat contains the combined information in the following way:
nB[1/fm^3], EcenN2LO[MeV/fm^3], EsigmaN2LO[MeV/fm^3], PcenN2LO[MeV/fm^3], PsigmaN2LO[MeV/fm^3], EcenN3LO[MeV/fm^3], EsigmaN3LO[MeV/fm^3], PcenN3LO[MeV/fm^3], PsigmaN3LO[MeV/fm^3]

In addition, the files ceft2NLO.dat and ceft3NLO.dat contain the minimal information necessary to evaluate the probabilities in the poltropic parametrization at low densities.
 nB[1/fm^3], Pcen[MeV/fm^3], Psigma[MeV/fm^3]
