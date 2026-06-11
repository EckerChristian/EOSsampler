#include "config.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace eos {

namespace {

std::string trim(std::string s) {
    auto not_space = [](unsigned char ch) {
        return !std::isspace(ch);
    };

    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

bool parse_bool(const std::string& value) {
    if (value == "true" || value == "True" || value == "1") return true;
    if (value == "false" || value == "False" || value == "0") return false;
    throw std::runtime_error("Invalid boolean value: " + value);
}

std::string strip_quotes(std::string value) {
    value = trim(value);
    if (value.size() >= 2) {
        if ((value.front() == '"' && value.back() == '"') ||
            (value.front() == '\'' && value.back() == '\'')) {
            return value.substr(1, value.size() - 2);
        }
    }
    return value;
}

std::unordered_map<std::string, std::string> read_key_values(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open config file: " + path);
    }

    std::unordered_map<std::string, std::string> values;
    std::string line;

    while (std::getline(in, line)) {
        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        line = trim(line);
        if (line.empty()) continue;

        const auto colon = line.find(':');
        if (colon == std::string::npos) {
            throw std::runtime_error("Invalid config line: " + line);
        }

        std::string key = trim(line.substr(0, colon));
        std::string value = strip_quotes(line.substr(colon + 1));

        values[key] = value;
    }

    return values;
}

template <class T>
void set_number(const std::unordered_map<std::string, std::string>& values,
                const std::string& key,
                T& target) {
    auto it = values.find(key);
    if (it == values.end()) return;

    std::istringstream ss(it->second);
    ss >> target;
    if (!ss) {
        throw std::runtime_error("Invalid numeric value for key: " + key);
    }
}

void set_string(const std::unordered_map<std::string, std::string>& values,
                const std::string& key,
                std::string& target) {
    auto it = values.find(key);
    if (it != values.end()) {
        target = it->second;
    }
}

void set_bool(const std::unordered_map<std::string, std::string>& values,
              const std::string& key,
              bool& target) {
    auto it = values.find(key);
    if (it != values.end()) {
        target = parse_bool(it->second);
    }
}

} // namespace

Config load_config(const std::string& path) {
    Config cfg;
    const auto values = read_key_values(path);

    set_string(values, "output_prefix", cfg.output_prefix);
    set_number(values, "file_start_index", cfg.file_start_index);
    set_number(values, "n_eos", cfg.n_eos);
    set_bool(values, "verbose", cfg.verbose);

    set_number(values, "n_mu", cfg.n_mu);
    set_number(values, "n_mu_pqcd", cfg.n_mu_pqcd);
    set_number(values, "n_cs2", cfg.n_cs2);
    set_number(values, "n_poly", cfg.n_poly);
    set_number(values, "n_pqcd", cfg.n_pqcd);

    set_bool(values, "include_phase_transition", cfg.include_phase_transition);
    set_bool(values, "write_eos_extension", cfg.write_eos_extension);
    set_number(values, "eos_extension_max_delta_mu_mev", cfg.eos_extension_max_delta_mu_mev);

    set_bool(values, "impose_pqcd", cfg.impose_pqcd);
    set_bool(values, "pqcd_cs2_condition", cfg.pqcd_cs2_condition);
    set_bool(values, "sample_cmax", cfg.sample_cmax);

    set_number(values, "cmin", cfg.cmin);
    set_number(values, "cmax", cfg.cmax);

    set_number(values, "n_cet", cfg.n_cet);
    set_number(values, "p_cet_soft", cfg.p_cet_soft);
    set_number(values, "p_cet_stiff", cfg.p_cet_stiff);

    set_number(values, "n_qcd", cfg.n_qcd);
    set_number(values, "mu_qcd", cfg.mu_qcd);
    set_number(values, "mu_cs2_cond", cfg.mu_cs2_cond);
    set_number(values, "dp_qcd", cfg.dp_qcd);
    set_number(values, "dn_qcd", cfg.dn_qcd);

    set_number(values, "n_pt_min_nsat", cfg.n_pt_min_nsat);
    set_number(values, "n_pt_max_nsat", cfg.n_pt_max_nsat);
    set_number(values, "dn_min_nsat", cfg.dn_min_nsat);
    set_number(values, "dn_max_nsat", cfg.dn_max_nsat);

    set_string(values, "crust_eos_file", cfg.crust_eos_file);
    set_number(values, "n_crust", cfg.n_crust);

    set_string(values, "ceft_file", cfg.ceft_file);
    set_number(values, "n_ceft", cfg.n_ceft);

    set_string(values, "pqcd_file", cfg.pqcd_file);
    set_number(values, "n_pqcd_file", cfg.n_pqcd_file);

    if (cfg.n_mu < 2) {
        throw std::runtime_error("n_mu must be >= 2");
    }

    if (cfg.file_start_index < 0) {
        throw std::runtime_error("file_start_index must be non-negative");
    }

    if (cfg.n_eos <= 0) {
        throw std::runtime_error("n_eos must be positive");
    }

    if (cfg.len_eos() <= 0) {
        throw std::runtime_error("Invalid EOS length");
    }

    set_number(values, "n_tov", cfg.n_tov);
    set_number(values, "n_radial", cfg.n_radial);
    set_number(values, "pc_initial", cfg.pc_initial);
    set_number(values, "pc_step", cfg.pc_step);
    set_number(values, "p_surface_cut", cfg.p_surface_cut);
    set_number(values, "r_initial", cfg.r_initial);
    set_number(values, "r_final", cfg.r_final);
    
    set_number(values, "mtov_min", cfg.mtov_min);
    set_number(values, "mtov_max", cfg.mtov_max);
    set_number(values, "pulsar_mass_mu", cfg.pulsar_mass_mu);
    set_number(values, "pulsar_mass_sigma", cfg.pulsar_mass_sigma);
    set_number(values, "mass_likelihood_points", cfg.mass_likelihood_points);

    if (cfg.n_tov <= 1) {
        throw std::runtime_error("n_tov must be greater than 1");
    }
    
    if (cfg.n_radial <= 10) {
        throw std::runtime_error("n_radial must be greater than 10");
    }
    
    if (cfg.r_final <= cfg.r_initial) {
        throw std::runtime_error("r_final must be larger than r_initial");
    }
    
    if (cfg.pc_initial <= 0.0 || cfg.pc_step <= 0.0) {
        throw std::runtime_error("pc_initial and pc_step must be positive");
    }

    set_bool(values, "enable_mass_likelihood", cfg.enable_mass_likelihood);
    set_bool(values, "enable_gw_likelihood", cfg.enable_gw_likelihood);
    set_bool(values, "enable_xray_likelihoods", cfg.enable_xray_likelihoods);
    
    set_number(values, "pulsar_mass_mu", cfg.pulsar_mass_mu);
    set_number(values, "pulsar_mass_sigma", cfg.pulsar_mass_sigma);
    set_number(values, "chirp_mass", cfg.chirp_mass);
    set_number(values, "mass_likelihood_points", cfg.mass_likelihood_points);
    
    set_string(values, "nicer0030_file", cfg.nicer0030_file);
    set_string(values, "nicer0437_file", cfg.nicer0437_file);
    set_string(values, "nicer0740_file", cfg.nicer0740_file);
    set_string(values, "nicer0614_file", cfg.nicer0614_file);
    
    set_string(values, "nat17_file", cfg.nat17_file);
    set_string(values, "nks1724b_file", cfg.nks1724b_file);
    set_string(values, "nks1810b_file", cfg.nks1810b_file);
    
    set_string(values, "shsM13_file", cfg.shsM13_file);
    
    set_string(values, "shb6304_file", cfg.shb6304_file);
    set_string(values, "shb6397_file", cfg.shb6397_file);
    set_string(values, "shbM28_file", cfg.shbM28_file);
    set_string(values, "shbM30_file", cfg.shbM30_file);
    set_string(values, "shbX7_file", cfg.shbX7_file);
    set_string(values, "shbX5_file", cfg.shbX5_file);
    set_string(values, "shbwCen_file", cfg.shbwCen_file);
    
    set_string(values, "shb_dataset", cfg.shb_dataset);
    
    set_string(values, "gw170817_file", cfg.gw170817_file);

    return cfg;
}

} // namespace eos
