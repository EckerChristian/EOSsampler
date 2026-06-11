#include "tov.hpp"

#include "constants.hpp"
#include "hdf5_io.hpp"
#include "likelihoods.hpp"
#include "schema.hpp"

#include <mpi.h>
#include <hdf5.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace eos {
namespace {

struct EOSTable {
    std::vector<double> e;     // GeV^4 for TOV internals.
    std::vector<double> p;     // GeV^4 for TOV internals.
    std::vector<double> n;     // fm^-3.
    std::vector<double> mu;    // MeV.
    std::vector<double> cs2;
};

struct TransitionInfo {
    double p_pt_mevfm3 = 0.0;
    double e_left_mevfm3 = 0.0;
    double n_left_fm3 = 0.0;
    double mu_pt_gev = 0.0;
    double dn_nsat = 0.0;

    bool has_phase_transition() const {
        return p_pt_mevfm3 > 0.0 && e_left_mevfm3 > 0.0 &&
               n_left_fm3 > 0.0 && mu_pt_gev > 0.0 && dn_nsat > 0.0;
    }

    double e_right_mevfm3() const {
        const double dn_fm3 = dn_nsat * ns;
        const double mu_pt_mev = mu_pt_gev * 1000.0;
        return -p_pt_mevfm3 + mu_pt_mev * (n_left_fm3 + dn_fm3);
    }

    double delta_e_mevfm3() const {
        return e_right_mevfm3() - e_left_mevfm3;
    }

    bool exceeds_seidov_limit() const {
        if (!has_phase_transition()) return false;
        const double delta = delta_e_mevfm3();
        if (!finite(delta) || delta <= 0.0) return false;
        const double critical = 0.5 + 1.5 * p_pt_mevfm3 / e_left_mevfm3;
        return delta / e_left_mevfm3 > critical;
    }
};

bool finite(double x) {
    return std::isfinite(x);
}

std::size_t find_interval(const std::vector<double>& x, double xq) {
    if (x.size() < 2) {
        throw std::runtime_error("Interpolation table is too short");
    }

    if (xq <= x.front()) return 0;
    if (xq >= x.back()) return x.size() - 2;

    auto it = std::upper_bound(x.begin(), x.end(), xq);
    std::size_t i = static_cast<std::size_t>(std::distance(x.begin(), it));
    if (i == 0) return 0;
    if (i >= x.size()) return x.size() - 2;
    return i - 1;
}

bool is_strictly_increasing(const std::vector<double>& x) {
    if (x.size() < 2) return false;
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) return false;
    }
    return true;
}

bool is_non_decreasing(const std::vector<double>& x) {
    if (x.size() < 2) return false;
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (x[i] < x[i - 1]) return false;
    }
    return true;
}

void validate_eos_table(const EOSTable& tab) {
    const std::size_t len = tab.e.size();
    if (len < 2 || tab.p.size() != len || tab.n.size() != len ||
        tab.mu.size() != len || tab.cs2.size() != len) {
        throw std::runtime_error("EOS datasets have inconsistent lengths");
    }

    for (std::size_t i = 0; i < len; ++i) {
        if (!finite(tab.e[i]) || !finite(tab.p[i]) || !finite(tab.n[i]) ||
            !finite(tab.mu[i]) || !finite(tab.cs2[i])) {
            throw std::runtime_error("EOS table contains non-finite values");
        }
    }

    // A first-order phase transition has a constant-pressure segment, so p must
    // be non-decreasing, not strictly increasing.  Energy density must still be
    // strictly increasing because n(e), mu(e), and cs2(e) are interpolated.
    if (!is_non_decreasing(tab.p)) {
        throw std::runtime_error("EOS pressure table must be non-decreasing for TOV interpolation");
    }

    if (!is_strictly_increasing(tab.e)) {
        throw std::runtime_error("EOS energy-density table must be strictly increasing for TOV interpolation");
    }
}

double linear_interp(const std::vector<double>& x, const std::vector<double>& y, double xq) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("Invalid interpolation table");
    }

    const std::size_t i = find_interval(x, xq);
    const double dx = x[i + 1] - x[i];
    if (dx == 0.0) return y[i];

    const double t = (xq - x[i]) / dx;
    return (1.0 - t) * y[i] + t * y[i + 1];
}

double log_interp(const std::vector<double>& x, const std::vector<double>& y, double xq) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::runtime_error("Invalid interpolation table");
    }

    if (xq <= 0.0) return y.front();

    const std::size_t i = find_interval(x, xq);
    if (x[i] <= 0.0 || x[i + 1] <= 0.0 || y[i] <= 0.0 || y[i + 1] <= 0.0) {
        return linear_interp(x, y, xq);
    }

    const double lx = std::log(xq);
    const double lx1 = std::log(x[i]);
    const double lx2 = std::log(x[i + 1]);
    const double ly1 = std::log(y[i]);
    const double ly2 = std::log(y[i + 1]);

    const double dx = lx2 - lx1;
    if (dx == 0.0) return y[i];

    const double ly = ((lx2 - lx) * ly1 + (lx - lx1) * ly2) / dx;
    return std::exp(ly);
}

double energy_from_pressure(double p, const EOSTable& eos_table) {
    return log_interp(eos_table.p, eos_table.e, p);
}

double dedp_from_pressure(double p, const EOSTable& eos_table) {
    std::size_t i = find_interval(eos_table.p, p);

    // Skip the zero-width pressure interval at a first-order phase transition.
    if (eos_table.p[i + 1] == eos_table.p[i]) {
        if (i + 2 < eos_table.p.size() && eos_table.p[i + 2] > eos_table.p[i + 1]) {
            ++i;
        } else if (i > 0 && eos_table.p[i] > eos_table.p[i - 1]) {
            --i;
        } else {
            return 0.0;
        }
    }

    const double dp = eos_table.p[i + 1] - eos_table.p[i];
    if (dp == 0.0) return 0.0;
    return (eos_table.e[i + 1] - eos_table.e[i]) / dp;
}

double number_from_energy(double e, const EOSTable& eos_table) {
    return linear_interp(eos_table.e, eos_table.n, e);
}

double mu_from_energy(double e, const EOSTable& eos_table) {
    return linear_interp(eos_table.e, eos_table.mu, e);
}

double cs2_from_energy(double e, const EOSTable& eos_table) {
    return linear_interp(eos_table.e, eos_table.cs2, e);
}

void rhs(
    const std::array<double, 4>& y,
    std::array<double, 4>& f,
    double r,
    double p_surface_cut,
    const EOSTable& eos_table,
    bool& surface_reached
) {
    const double p = y[0];
    const double m = y[1];
    const double H = y[2];
    const double dH = y[3];

    if (p < p_surface_cut || r <= 0.0 || !finite(p) || !finite(m)) {
        f = {0.0, 0.0, 0.0, 0.0};
        surface_reached = true;
        return;
    }

    surface_reached = false;

    const double e = energy_from_pressure(p, eos_table);
    const double dedp = dedp_from_pressure(p, eos_table);
    const double denom = r * (r - 2.0 * Gh * m);
    const double metric = 1.0 - 2.0 * Gh * m / r;

    if (denom == 0.0 || metric <= 0.0 || !finite(e) || !finite(dedp)) {
        f = {0.0, 0.0, 0.0, 0.0};
        surface_reached = true;
        return;
    }

    f[0] = -Gh * (p + e) * (4.0 * pi * p * std::pow(r, 3.0) + m) / denom;
    f[1] = 4.0 * pi * e * r * r;
    f[2] = dH;

    const double term1 =
        -(2.0 / r + Gh * (2.0 * m / (r * r) + 4.0 * pi * r * (p - e)) / metric) * dH;
    const double term2 = 6.0 / (r * r * metric);
    const double term3 =
        std::pow(2.0 * Gh * (m + 4.0 * pi * std::pow(r, 3.0) * p) / (r * r * metric), 2.0);
    const double term4 =
        4.0 * Gh * pi * (9.0 * p + 5.0 * e + (p + e) * dedp) / metric;

    f[3] = term1 + (term2 + term3 - term4) * H;
}

std::array<double, 4> rk4_step(
    const std::array<double, 4>& y,
    double r,
    double dr,
    double p_surface_cut,
    const EOSTable& eos_table,
    bool& surface_reached
) {
    std::array<double, 4> k1{}, k2{}, k3{}, k4{}, yt{};

    rhs(y, k1, r, p_surface_cut, eos_table, surface_reached);
    if (surface_reached) return y;

    for (int i = 0; i < 4; ++i) yt[i] = y[i] + 0.5 * dr * k1[i];
    rhs(yt, k2, r + 0.5 * dr, p_surface_cut, eos_table, surface_reached);
    if (surface_reached) return y;

    for (int i = 0; i < 4; ++i) yt[i] = y[i] + 0.5 * dr * k2[i];
    rhs(yt, k3, r + 0.5 * dr, p_surface_cut, eos_table, surface_reached);
    if (surface_reached) return y;

    for (int i = 0; i < 4; ++i) yt[i] = y[i] + dr * k3[i];
    rhs(yt, k4, r + dr, p_surface_cut, eos_table, surface_reached);
    if (surface_reached) return y;

    std::array<double, 4> out{};
    for (int i = 0; i < 4; ++i) {
        out[i] = y[i] + (dr / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return out;
}

double tidal_lambda(const std::array<double, 4>& y, double r) {
    if (r <= 0.0 || y[2] == 0.0) return 0.0;

    const double C = Gh * y[1] / r;
    if (C <= 0.0 || C >= 0.5) return 0.0;

    const double yy = r * y[3] / y[2];
    const double one_minus_2c = 1.0 - 2.0 * C;

    const double numerator =
        16.0 / 15.0 * std::pow(one_minus_2c, 2.0) *
        (2.0 + 2.0 * C * (yy - 1.0) - yy);

    const double denominator =
        3.0 * std::pow(one_minus_2c, 2.0) *
            (2.0 - yy + 2.0 * C * (yy - 1.0)) * std::log(one_minus_2c)
        + 4.0 * std::pow(C, 3.0) *
            (13.0 - 11.0 * yy + C * (3.0 * yy - 2.0) + 2.0 * C * C * (1.0 + yy))
        + 2.0 * C * (6.0 - 3.0 * yy + 3.0 * C * (5.0 * yy - 8.0));

    if (denominator == 0.0 || !finite(numerator) || !finite(denominator)) return 0.0;

    const double lambda = numerator / denominator;
    return finite(lambda) && lambda > 0.0 ? lambda : 0.0;
}

void compute_stability(TOVSequence& seq, bool force_unstable_at_pt) {
    const std::size_t n = seq.mass.size();
    seq.stability.assign(n, 0.0);
    if (n == 0) return;

    constexpr double abs_tol = 1.0e-4;  // solar masses; suppresses one-step noise.
    constexpr double rel_tol = 1.0e-4;
    constexpr double min_branch_mass_span = 1.0e-2;
    constexpr std::size_t min_branch_points = 3;

    auto dm_sign = [&](std::size_t i) {
        const double dm = seq.mass[i] - seq.mass[i - 1];
        const double scale = std::max({1.0, std::abs(seq.mass[i]), std::abs(seq.mass[i - 1])});
        const double tol = std::max(abs_tol, rel_tol * scale);
        if (dm > tol) return 1;
        if (dm < -tol) return -1;
        return 0;
    };

    seq.stability[0] = 1.0;
    int branch = 1;
    bool stable = true;
    bool seidov_forced_unstable = false;
    bool saw_decreasing_segment_after_pt = false;

    for (std::size_t i = 1; i < n; ++i) {
        const int sign = dm_sign(i);
        const bool has_qm_core = i < seq.r_qm.size() && seq.r_qm[i] > 0.0;

        // If the density jump exceeds the Seidov criterion, the connected branch
        // becomes unstable at the onset of an infinitesimal QM core.  A stable
        // hybrid branch is allowed only after the sequence has first gone through
        // a decreasing-mass segment and then turns upward again.
        if (force_unstable_at_pt && has_qm_core && branch == 1 && !seidov_forced_unstable) {
            stable = false;
            seidov_forced_unstable = true;
            seq.stability[i] = -1.0;
            continue;
        }

        if (seidov_forced_unstable && !stable && branch == 1) {
            if (sign < 0) saw_decreasing_segment_after_pt = true;
            if (sign > 0 && saw_decreasing_segment_after_pt) {
                stable = true;
                branch = 2;
                seq.stability[i] = 2.0;
            } else {
                seq.stability[i] = -1.0;
            }
            continue;
        }

        if (sign > 0 && stable) {
            seq.stability[i] = static_cast<double>(branch);
        } else if (sign < 0 && stable) {
            stable = false;
            seq.stability[i] = -static_cast<double>(branch);
        } else if (sign > 0 && !stable) {
            stable = true;
            ++branch;
            seq.stability[i] = static_cast<double>(branch);
        } else if (sign < 0 && !stable) {
            seq.stability[i] = -static_cast<double>(branch);
        } else {
            seq.stability[i] = stable ? static_cast<double>(branch) : -static_cast<double>(branch);
        }
    }

    // Remove tiny one- or two-point "stable branches" caused by sampling noise.
    const int max_branch = branch;
    for (int b = 1; b <= max_branch; ++b) {
        std::vector<std::size_t> idx;
        for (std::size_t i = 0; i < n; ++i) {
            if (static_cast<int>(std::round(seq.stability[i])) == b) idx.push_back(i);
        }
        if (idx.empty()) continue;

        const double span = seq.mass[idx.back()] - seq.mass[idx.front()];
        if (idx.size() < min_branch_points || span < min_branch_mass_span) {
            for (std::size_t i : idx) seq.stability[i] = -static_cast<double>(b);
        }
    }
}


double max_stable_mass(const TOVSequence& seq) {
    double mtov = 0.0;
    for (std::size_t i = 0; i < seq.mass.size(); ++i) {
        if (seq.stability.empty() || seq.stability[i] > 0.0) {
            mtov = std::max(mtov, seq.mass[i]);
        }
    }
    return mtov;
}

EOSTable read_eos_table(hid_t group_id) {
    hid_t eos_group = hdf5::open_subgroup(group_id, std::string(schema::eos_group));

    EOSTable tab;
    tab.e = hdf5::read_double_array(eos_group, "e");
    tab.p = hdf5::read_double_array(eos_group, "p");
    tab.n = hdf5::read_double_array(eos_group, "n");
    tab.mu = hdf5::read_double_array(eos_group, "mu");
    tab.cs2 = hdf5::read_double_array(eos_group, "cs2");

    H5Gclose(eos_group);

    const std::size_t len = tab.e.size();
    if (len < 2 || tab.p.size() != len || tab.n.size() != len ||
        tab.mu.size() != len || tab.cs2.size() != len) {
        throw std::runtime_error("EOS datasets have inconsistent lengths");
    }

    // EOS files store e and p in MeV/fm^3. The TOV equations here use GeV^4.
    for (std::size_t i = 0; i < len; ++i) {
        tab.e[i] *= MeVbyfm3ToGev4;
        tab.p[i] *= MeVbyfm3ToGev4;
    }

    validate_eos_table(tab);
    return tab;
}

TOVSequence solve_sequence(const EOSTable& eos_table, const Config& cfg, const TransitionInfo& transition) {
    TOVSequence seq;
    seq.mass.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.radius.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.lambda.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.e_cent.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.mu_cent.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.n_cent.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.cs2_cent.reserve(static_cast<std::size_t>(cfg.n_tov));
    seq.r_qm.reserve(static_cast<std::size_t>(cfg.n_tov));

    const double dr = (cfg.r_final - cfg.r_initial) / static_cast<double>(cfg.n_radial - 1);
    const double p_pt = transition.p_pt_mevfm3 * MeVbyfm3ToGev4;

    double pc = cfg.pc_initial * MeVbyfm3ToGev4;
    double dpc = cfg.pc_step * MeVbyfm3ToGev4;

    for (int k = 0; k < cfg.n_tov; ++k) {
        if (k > 0) pc += dpc;
        if (pc <= 0.0) continue;

        const double ec = energy_from_pressure(pc, eos_table);
        if (!finite(ec) || ec <= 0.0) break;

        std::array<double, 4> y{
            pc,
            (4.0 / 3.0) * pi * std::pow(cfg.r_initial, 3.0) * ec,
            cfg.r_initial * cfg.r_initial,
            2.0 * cfg.r_initial
        };

        double r = cfg.r_initial;
        bool surface_reached = false;
        const bool crosses_pt = cfg.include_phase_transition && pc > p_pt && p_pt > 0.0;
        bool crossed_pt = false;
        double r_qm = 0.0;

        std::array<double, 4> previous = y;
        double previous_r = r;

        for (int i = 0; i < cfg.n_radial - 1; ++i) {
            previous = y;
            previous_r = r;

            std::array<double, 4> next = rk4_step(y, r, dr, cfg.p_surface_cut, eos_table, surface_reached);
            if (surface_reached) break;

            if (crosses_pt && !crossed_pt && y[0] >= p_pt && next[0] < p_pt) {
                crossed_pt = true;
                const double denom = y[0] - next[0];
                r_qm = denom != 0.0 ? r + dr * (y[0] - p_pt) / denom : r;
            }

            y = next;
            r += dr;

            if (!finite(y[0]) || !finite(y[1]) || !finite(y[2]) || !finite(y[3])) {
                surface_reached = true;
                break;
            }
        }

        const auto& surface_state = surface_reached ? previous : y;
        const double surface_r = surface_reached ? previous_r : r;
        const double mass = surface_state[1] / M_solarh;
        const double lambda = tidal_lambda(surface_state, surface_r);

        if (!finite(mass) || !finite(surface_r) || mass <= 0.0 || surface_r <= 0.0) break;

        seq.mass.push_back(mass);
        seq.radius.push_back(surface_r);
        seq.lambda.push_back(lambda);
        seq.e_cent.push_back(ec / MeVbyfm3ToGev4);
        seq.mu_cent.push_back(mu_from_energy(ec, eos_table));
        seq.n_cent.push_back(number_from_energy(ec, eos_table));
        seq.cs2_cent.push_back(cs2_from_energy(ec, eos_table));
        seq.r_qm.push_back(r_qm);

        if (k > 0 && seq.mass.size() >= 2) {
            const double dm = std::abs(seq.mass.back() - seq.mass[seq.mass.size() - 2]);
            if (dm < 0.05) dpc *= 1.1;
            if (dm > 0.06) dpc /= 1.1;
        }

        if (!seq.n_cent.empty() && seq.n_cent.back() >= cfg.n_qcd * ns) break;
    }

    compute_stability(seq, transition.exceeds_seidov_limit());
    return seq;
}

TransitionInfo read_transition_info(hid_t group_id) {
    TransitionInfo transition;
    if (!hdf5::link_exists(group_id, std::string(schema::params_group))) {
        return transition;
    }

    hid_t params_group = hdf5::open_subgroup(group_id, std::string(schema::params_group));
    if (hdf5::link_exists(params_group, "pPT")) {
        transition.p_pt_mevfm3 = hdf5::read_double_scalar(params_group, "pPT");
    }
    if (hdf5::link_exists(params_group, "ePTl")) {
        transition.e_left_mevfm3 = hdf5::read_double_scalar(params_group, "ePTl");
    }
    if (hdf5::link_exists(params_group, "nPTl")) {
        transition.n_left_fm3 = hdf5::read_double_scalar(params_group, "nPTl");
    }
    if (hdf5::link_exists(params_group, "muPT")) {
        transition.mu_pt_gev = hdf5::read_double_scalar(params_group, "muPT");
    }
    if (hdf5::link_exists(params_group, "dn")) {
        transition.dn_nsat = hdf5::read_double_scalar(params_group, "dn");
    }
    H5Gclose(params_group);
    return transition;
}

void write_tov_sequence(hid_t group_id, const TOVSequence& seq) {
    hid_t tov_group;
    if (hdf5::link_exists(group_id, std::string(schema::tov_group))) {
        tov_group = hdf5::open_subgroup(group_id, std::string(schema::tov_group));
    } else {
        tov_group = hdf5::create_subgroup(group_id, std::string(schema::tov_group));
    }

    const hsize_t len = static_cast<hsize_t>(seq.mass.size());

    hdf5::recreate_dataset_double(tov_group, "M", len);
    hdf5::recreate_dataset_double(tov_group, "R", len);
    hdf5::recreate_dataset_double(tov_group, "Lambda", len);
    hdf5::recreate_dataset_double(tov_group, "stab", len);
    hdf5::recreate_dataset_double(tov_group, "e_cent", len);
    hdf5::recreate_dataset_double(tov_group, "mu_cent", len);
    hdf5::recreate_dataset_double(tov_group, "n_cent", len);
    hdf5::recreate_dataset_double(tov_group, "cs2_cent", len);
    hdf5::recreate_dataset_double(tov_group, "Rqm", len);

    hdf5::write_double_array(tov_group, "M", seq.mass);
    hdf5::write_double_array(tov_group, "R", seq.radius);
    hdf5::write_double_array(tov_group, "Lambda", seq.lambda);
    hdf5::write_double_array(tov_group, "stab", seq.stability);
    hdf5::write_double_array(tov_group, "e_cent", seq.e_cent);
    hdf5::write_double_array(tov_group, "mu_cent", seq.mu_cent);
    hdf5::write_double_array(tov_group, "n_cent", seq.n_cent);
    hdf5::write_double_array(tov_group, "cs2_cent", seq.cs2_cent);
    hdf5::write_double_array(tov_group, "Rqm", seq.r_qm);

    H5Gclose(tov_group);
}

void ensure_param_datasets(hid_t params_group) {
    if (!hdf5::link_exists(params_group, "pM")) hdf5::create_dataset_double(params_group, "pM", 1);
    if (!hdf5::link_exists(params_group, "pGW")) hdf5::create_dataset_double(params_group, "pGW", 1);
    if (!hdf5::link_exists(params_group, "ptot")) hdf5::create_dataset_double(params_group, "ptot", 1);
    if (!hdf5::link_exists(params_group, "nbranches")) hdf5::create_dataset_int(params_group, "nbranches", 1);
    if (!hdf5::link_exists(params_group, "pXray")) hdf5::create_dataset_double(params_group, "pXray", 15);
}

void write_likelihoods(
    hid_t group_id,
    const TOVSequence& seq,
    const LikelihoodEvaluator& likelihood_evaluator
) {
    hid_t params_group = hdf5::open_subgroup(group_id, std::string(schema::params_group));

    const double p_cet = hdf5::link_exists(params_group, "pCET")
        ? hdf5::read_double_scalar(params_group, "pCET")
        : 1.0;
    const double p_pqcd = hdf5::link_exists(params_group, "pPQCD")
        ? hdf5::read_double_scalar(params_group, "pPQCD")
        : 1.0;

    BranchInfo branches = find_branches(seq);
    LikelihoodResult like = likelihood_evaluator.evaluate(seq, branches);

    const double p_total = p_cet * p_pqcd * like.p_mass * like.p_gw * like.p_xray;
    const int nbranches = static_cast<int>(branches.idx_start.size());

    ensure_param_datasets(params_group);

    hdf5::write_double_scalar(params_group, "pM", like.p_mass);
    hdf5::write_double_scalar(params_group, "pGW", like.p_gw);
    hdf5::write_double_scalar(params_group, "ptot", p_total);
    hdf5::write_int_scalar(params_group, "nbranches", nbranches);

    std::vector<double> p_xray_each = like.p_xray_each;
    if (p_xray_each.size() < 15) p_xray_each.resize(15, 0.0);
    if (p_xray_each.size() > 15) p_xray_each.resize(15);
    hdf5::write_double_array(params_group, "pXray", p_xray_each);

    H5Gclose(params_group);
}

} // namespace

int run_tov_solver(const Config& cfg) {
    int mpi_was_initialized = 0;
    MPI_Initialized(&mpi_was_initialized);
    if (!mpi_was_initialized) MPI_Init(nullptr, nullptr);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    (void)size;

    try {
	hid_t file_id = hdf5::open_file_rw(cfg.output_prefix, rank, cfg.file_start_index);

        // Construct once so external GW/NICER/X-ray tables are loaded once and reused.
        LikelihoodEvaluator likelihood_evaluator(cfg);

        for (int eos_index = 0; eos_index < cfg.n_eos; ++eos_index) {
            if (!hdf5::link_exists(file_id, std::to_string(eos_index))) {
                if (cfg.verbose && rank == 0) {
                    std::cerr << "Skipping missing EOS group " << eos_index << '\n';
                }
                continue;
            }

            hid_t group_id = hdf5::open_group(file_id, eos_index);
            EOSTable eos_table = read_eos_table(group_id);

            const TransitionInfo transition = read_transition_info(group_id);
            TOVSequence seq = solve_sequence(eos_table, cfg, transition);
            write_tov_sequence(group_id, seq);
            write_likelihoods(group_id, seq, likelihood_evaluator);

            if (cfg.verbose) {
                BranchInfo branches = find_branches(seq);
                std::cout << "rank=" << rank
                          << ", EOS=" << eos_index
                          << ", points=" << seq.mass.size()
                          << ", Mtov=" << max_stable_mass(seq)
                          << ", nbranches=" << branches.idx_start.size()
                          << '\n';
            }

            H5Gclose(group_id);
        }

        H5Fclose(file_id);
        if (!mpi_was_initialized) MPI_Finalize();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[rank " << rank << "] TOV failed: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    } catch (...) {
        std::cerr << "[rank " << rank << "] TOV failed with an unknown error" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}

} // namespace eos
