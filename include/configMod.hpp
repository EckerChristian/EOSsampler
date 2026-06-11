#pragma once

#include <string>

namespace eos {

struct Config {
    std::string output_prefix = "./";
    int file_start_index = 0;
    int n_eos = 100;
    bool verbose = false;

    int n_mu = 10;
    int n_mu_pqcd = 2;
    int n_cs2 = 300;
    int n_low_cs2 = 4;
    double low_cs2_end_min = 0.001;
    double low_cs2_end_max = 0.1;
    int n_pqcd = 3;

    bool include_phase_transition = true;
    bool impose_pqcd = true;
    bool pqcd_cs2_condition = true;
    bool sample_cmax = false;

    double cmin = 1.0e-10;
    double cmax = 1.0;

    double n_cet = 1.1;
    double p_cet_soft = 1.5;
    double p_cet_stiff = 7.0;

    double n_qcd = 40.0;
    double mu_qcd = 2.6;
    double mu_cs2_cond = 2.25;
    double dp_qcd = 5.0;
    double dn_qcd = 0.25;

    double n_pt_min_nsat = 2.0;
    double n_pt_max_nsat = 7.0;
    double dn_min_nsat = 0.1;
    double dn_max_nsat = 2.0;


    // TOV sequence settings.
    int n_tov = 200;
    int n_radial = 2000;
    double pc_initial = 5.0;        // MeV/fm^3
    double pc_step = 0.5;           // MeV/fm^3
    double p_surface_cut = 1.0e-16; // GeV^4
    double r_initial = 1.0e-5;      // km
    double r_final = 18.0;          // km

    // TOV-side constraints.
    double mtov_min = 2.0;
    double mtov_max = 10.0;

    std::string crust_eos_file = "../crust/BPS2.txt";
    int n_crust = 77;

    std::string ceft_file = "../ceft/ceft2NLO_1.1ns.dat";
    int n_ceft = 25;

    std::string pqcd_file = "../pQCD/cs2_muB_mean_std_scale_ave.dat";
    int n_pqcd_file = 100;

    int len_eos() const {
        return n_crust + n_low_cs2 + n_cs2 + n_pqcd;
    }

    // Likelihood switches.
    bool enable_mass_likelihood = true;
    bool enable_gw_likelihood = true;
    bool enable_xray_likelihoods = true;
    
    // Likelihood parameters.
    double pulsar_mass_mu = 1.908;
    double pulsar_mass_sigma = 0.016;
    double chirp_mass = 1.186;
    int mass_likelihood_points = 200;
    
    // NICER paths.
    std::string nicer0030_file = "../constraints/mrdata/nicer_J0030/0030_maryland_2spts_regrid2.h5";
    std::string nicer0437_file = "../constraints/mrdata/nicer_J0437/J0437_regrid_new.h5";
    std::string nicer0740_file = "../constraints/mrdata/nicer_J0740_update/J0740_update_regrid_new.h5";
    std::string nicer0614_file = "../constraints/mrdata/nicer_J0614/J0614_regrid_new.h5";
    
    // X-ray paths.
    std::string nat17_file = "../constraints/mrdata/nat17/1702_D_X_int.o2";
    std::string nks1724b_file = "../constraints/mrdata/nks15/1724b.o2";
    std::string nks1810b_file = "../constraints/mrdata/nks15/1810b.o2";
    
    std::string shsM13_file = "../constraints/mrdata/shs18/M13_H_rs.o2";
    
    std::string shb6304_file = "../constraints/mrdata/shb18/6304_He_nopl_syst_wilm.o2";
    std::string shb6397_file = "../constraints/mrdata/shb18/6397_He_syst_wilm3.o2";
    std::string shbM28_file = "../constraints/mrdata/shb18/M28_He_syst_wilm.o2";
    std::string shbM30_file = "../constraints/mrdata/shb18/M30_H_syst_wilm.o2";
    std::string shbX7_file = "../constraints/mrdata/shb18/X7_H_syst_wilm.o2";
    std::string shbX5_file = "../constraints/mrdata/shb18/X5_H_syst_wilm.o2";
    std::string shbwCen_file = "../constraints/mrdata/shb18/wCen_H_syst_wilm.o2";
    
    std::string shb_dataset = "rescaled/data/like";
    
    // GW path.
    std::string gw170817_file = "../constraints/LV/LV_prior.h5";


};

Config load_config(const std::string& path);

} // namespace eos
