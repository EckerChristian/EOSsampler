#include "eos.hpp"

#include "constants.hpp"
#include "hdf5_io.hpp"
#include "schema.hpp"

#include <mpi.h>
#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace eos {
namespace {

struct CrustTable {
    std::vector<double> e;
    std::vector<double> p;
    std::vector<double> n_over_ns;
    std::vector<double> mu;
    std::vector<double> cs2;
};

struct CEFTTable {
    std::vector<double> n;
    std::vector<double> median_p;
    std::vector<double> sigma_p;
};

struct PQCDTable {
    std::vector<double> mu_b;
    std::vector<double> cs2_mean;
    std::vector<double> cs2_std;
};

struct EOSSample {
    // EOS that includes the PT
    std::vector<double> e;
    std::vector<double> p;
    std::vector<double> n;
    std::vector<double> mu;
    std::vector<double> cs2;
    // EOS that extends beyond the PT, but does not include the PT.
    std::vector<double> e_ext;
    std::vector<double> p_ext;
    std::vector<double> n_ext;
    std::vector<double> mu_ext;
    std::vector<double> cs2_ext;

    double dn = 0.0;
    double ePTl = 0.0;
    double pPT = 0.0;
    double nPTl = 0.0;
    double muPT = 0.0;
    double cs2PTr = 0.0;
    double cs2PTl = 0.0;

    double pCET = 1.0;
    double pPQCD = 1.0;
};

CrustTable load_crust(const Config& cfg) {
    if (cfg.n_crust < 2) {
        throw std::runtime_error("n_crust must be at least 2");
    }

    std::ifstream in(cfg.crust_eos_file);
    if (!in) {
        throw std::runtime_error("Could not open crust EOS file: " + cfg.crust_eos_file);
    }

    CrustTable table;
    table.e.resize(cfg.n_crust);
    table.p.resize(cfg.n_crust);
    table.n_over_ns.resize(cfg.n_crust);
    table.mu.resize(cfg.n_crust);
    table.cs2.resize(cfg.n_crust);

    for (int i = 0; i < cfg.n_crust; ++i) {
        in >> table.e[i] >> table.p[i] >> table.n_over_ns[i];
        if (!in) {
            throw std::runtime_error("Failed while reading crust EOS file");
        }

        table.mu[i] = (table.e[i] + table.p[i]) / (ns * table.n_over_ns[i]);
    }

    for (int i = 1; i < cfg.n_crust - 1; ++i) {
        table.cs2[i] =
            0.5 * (table.p[i + 1] - table.p[i - 1]) /
            (table.e[i + 1] - table.e[i - 1]);
    }

    table.cs2[0] =
        (table.p[1] - table.p[0]) /
        (table.e[1] - table.e[0]);

    table.cs2[cfg.n_crust - 1] =
        (table.p[cfg.n_crust - 1] - table.p[cfg.n_crust - 2]) /
        (table.e[cfg.n_crust - 1] - table.e[cfg.n_crust - 2]);

    return table;
}

CEFTTable load_ceft(const Config& cfg) {
    if (cfg.n_ceft < 2) {
        throw std::runtime_error("n_ceft must be at least 2");
    }

    std::ifstream in(cfg.ceft_file);
    if (!in) {
        throw std::runtime_error("Could not open CEFT file: " + cfg.ceft_file);
    }

    CEFTTable table;
    table.n.resize(cfg.n_ceft);
    table.median_p.resize(cfg.n_ceft);
    table.sigma_p.resize(cfg.n_ceft);

    for (int i = 0; i < cfg.n_ceft; ++i) {
        in >> table.n[i] >> table.median_p[i] >> table.sigma_p[i];
        if (!in) {
            throw std::runtime_error("Failed while reading CEFT file");
        }
    }

    return table;
}

PQCDTable load_pqcd(const Config& cfg) {
    if (cfg.n_pqcd_file < 1) {
        throw std::runtime_error("n_pqcd_file must be positive");
    }

    std::ifstream in(cfg.pqcd_file);
    if (!in) {
        throw std::runtime_error("Could not open pQCD file: " + cfg.pqcd_file);
    }

    std::string header;
    std::getline(in, header);

    PQCDTable table;
    table.mu_b.resize(cfg.n_pqcd_file);
    table.cs2_mean.resize(cfg.n_pqcd_file);
    table.cs2_std.resize(cfg.n_pqcd_file);

    for (int i = 0; i < cfg.n_pqcd_file; ++i) {
        in >> table.mu_b[i] >> table.cs2_mean[i] >> table.cs2_std[i];
        if (!in) {
            throw std::runtime_error("Failed while reading pQCD file");
        }
    }

    return table;
}

bool is_finite(double x) {
    return std::isfinite(x);
}

double uniform(std::mt19937& rng, double lo, double hi) {
    std::uniform_real_distribution<double> dist(lo, hi);
    return dist(rng);
}

double linear_interp(
    double x,
    const std::vector<double>& xs,
    const std::vector<double>& ys
) {
    if (xs.size() != ys.size() || xs.size() < 2) {
        throw std::runtime_error("Invalid interpolation table");
    }

    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();

    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        if (x >= xs[i] && x <= xs[i + 1]) {
            const double dx = xs[i + 1] - xs[i];
            if (dx == 0.0) return ys[i];
            const double t = (x - xs[i]) / dx;
            return (1.0 - t) * ys[i] + t * ys[i + 1];
        }
    }

    return ys.back();
}

// pQCD running of alpha_s, N3LO.
double alpha_s(double lbar) {
    const double L = 2.0 * std::log(lbar / LMSbar);
    return 4.0 * pi / beta0 *
           (1.0 / L
            - beta1 / std::pow(beta0, 2.0) / std::pow(L, 2.0) * std::log(L)
            + std::pow(beta1, 2.0) / std::pow(beta0, 4.0) / std::pow(L, 3.0) *
                  (std::pow(std::log(L), 2.0) - std::log(L) - 1.0)
            + beta2 / std::pow(beta0, 3.0) / std::pow(L, 3.0));
}

double pfree(double mub) {
    return std::pow(mub / 3.0, 4.0) * 9.0 / 12.0 / pi / pi;
}

double dpfree_dmub(double mub) {
    return 4.0 / 3.0 * std::pow(mub / 3.0, 3.0) * 9.0 / 12.0 / pi / pi;
}

double d2pfree_dmub2(double mub) {
    return 4.0 * 3.0 / 3.0 / 3.0 * std::pow(mub / 3.0, 2.0) * 9.0 / 12.0 / pi / pi;
}

double beta2_func_X_lbar(double lbar, double as) {
    return -2.0 * as *
           ((as / 4.0 / pi) * beta0
            + std::pow(as / 4.0 / pi, 2.0) * beta1
            + std::pow(as / 4.0 / pi, 3.0) * beta2) /
           lbar;
}

double Dbeta2_das_func_X_lbar(double lbar, double as) {
    return -2.0 *
           (2.0 * (as / 4.0 / pi) * beta0
            + 3.0 * std::pow(as / 4.0 / pi, 2.0) * beta1
            + 4.0 * std::pow(as / 4.0 / pi, 3.0) * beta2) /
           lbar;
}

double pressure_qcd(double mub, double X) {
    const double lbar = 2.0 * (mub / 3.0) * X;
    const double as = alpha_s(lbar);
    return 1000.0 * std::pow(gevtoinvfm, 3.0) * pfree(mub) *
           (1.0
            + (as / pi) * a11
            + 3.0 * std::pow(as / pi, 2.0) *
                  (a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X) + a23)
            + std::pow(3.0, 2.0) * std::pow(as / pi, 3.0) * a3);
}

double baryon_density_qcd(double mub, double X) {
    const double lbar = 2.0 * (mub / 3.0) * X;
    const double as = alpha_s(lbar);
    const double das_dmub = (lbar / mub) * beta2_func_X_lbar(lbar, as);
    const double p_as =
        (1.0
         + (as / pi) * a11
         + 3.0 * std::pow(as / pi, 2.0) *
               (a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X) + a23)
         + std::pow(3.0, 2.0) * std::pow(as / pi, 3.0) * a3);
    const double dp_das =
        a11 / pi
        + 3.0 * as / pi *
              (a21 / 2.0 + a23 + a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X)) *
              (2.0 / pi)
        + std::pow(3.0, 2.0) * std::pow(as / pi, 2.0) * a3 * (3.0 / pi);

    return std::pow(gevtoinvfm, 3.0) *
           (dp_das * das_dmub * pfree(mub) + p_as * dpfree_dmub(mub));
}

double cs2_qcd(double mub, double X) {
    const double lbar = 2.0 * (mub / 3.0) * X;
    const double as = alpha_s(lbar);
    const double das_dmub = (lbar / mub) * beta2_func_X_lbar(lbar, as);
    const double d2as_dmub2 =
        -(lbar / mub / mub) * beta2_func_X_lbar(lbar, as)
        + (lbar / mub) * Dbeta2_das_func_X_lbar(lbar, as) * das_dmub;

    const double p_as =
        (1.0
         + (as / pi) * a11
         + 3.0 * std::pow(as / pi, 2.0) *
               (a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X) + a23)
         + std::pow(3.0, 2.0) * std::pow(as / pi, 3.0) * a3);

    const double dp_das =
        a11 / pi
        + 3.0 * as / pi *
              (a21 / 2.0 + a23 + a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X)) *
              (2.0 / pi)
        + std::pow(3.0, 2.0) * std::pow(as / pi, 2.0) * a3 * (3.0 / pi);

    const double d2p_das2 =
        3.0 * (3.0 * a21 / 2.0 + a23 + a21 * std::log(3.0 * (as / pi)) + a22 * std::log(X)) *
              (2.0 / pi / pi)
        + std::pow(3.0, 2.0) * (as / pi) * a3 * (3.0 * 2.0 / pi / pi);

    const double n = dp_das * das_dmub * pfree(mub) + p_as * dpfree_dmub(mub);
    const double dn_dmub =
        pfree(mub) * (d2p_das2 * std::pow(das_dmub, 2.0) + dp_das * d2as_dmub2)
        + 2.0 * dp_das * das_dmub * dpfree_dmub(mub)
        + p_as * d2pfree_dmub2(mub);

    return n / mub / dn_dmub;
}

double cs_at_mu(
    double mu,
    const std::vector<double>& m,
    const std::vector<double>& c
) {
    if (m.size() != c.size() || m.size() < 2) {
        throw std::runtime_error("Invalid cs2 segment table");
    }

    if (mu <= m.front()) return c.front();
    if (mu >= m.back()) return c.back();

    for (std::size_t i = 0; i + 1 < m.size(); ++i) {
        if (mu >= m[i] && mu <= m[i + 1]) {
            const double dmu = m[i + 1] - m[i];
            if (dmu == 0.0) return c[i + 1];
            return (c[i] * (m[i + 1] - mu) + (mu - m[i]) * c[i + 1]) / dmu;
        }
    }

    return c.back();
}

double density_at_mu(
    double mu,
    const std::vector<double>& m,
    const std::vector<double>& c,
    double n0,
    double dn
) {
    if (m.size() != c.size() || m.size() < 2) {
        throw std::runtime_error("Invalid density segment table");
    }

    const int n_segments = static_cast<int>(m.size()) - 1;
    double value = n0 * ns;

    if (mu <= m.front()) return value;

    for (int i = 0; i < n_segments; ++i) {
        const bool in_segment = (mu >= m[i] && mu < m[i + 1]) ||
                                (i == n_segments - 1 && mu <= m[i + 1]);
        if (!in_segment) continue;

        for (int k = 0; k < i; ++k) {
            if (m[k] == m[k + 1]) {
                value += dn * ns;
                continue;
            }

            const double denom = c[k + 1] * m[k] - c[k] * m[k + 1];
            if (denom == 0.0) return std::numeric_limits<double>::quiet_NaN();

            value *= std::pow(
                (m[k] * c[k + 1]) / (m[k + 1] * c[k]),
                (m[k + 1] - m[k]) / denom
            );
        }

        if (m[i] == m[i + 1]) {
            value += dn * ns;
            return value;
        }

        const double denom = c[i + 1] * m[i] - c[i] * m[i + 1];
        if (denom == 0.0) return std::numeric_limits<double>::quiet_NaN();

        const double cs2_interp = c[i] * (m[i + 1] - mu) + (mu - m[i]) * c[i + 1];
        value *= std::pow(
            (m[i] * cs2_interp) / (mu * c[i] * (m[i + 1] - m[i])),
            (m[i + 1] - m[i]) / denom
        );

        return value;
    }

    // Above table range: return density at the last point.
    for (int k = 0; k < n_segments; ++k) {
        if (m[k] == m[k + 1]) {
            value += dn * ns;
            continue;
        }

        const double denom = c[k + 1] * m[k] - c[k] * m[k + 1];
        if (denom == 0.0) return std::numeric_limits<double>::quiet_NaN();

        value *= std::pow(
            (m[k] * c[k + 1]) / (m[k + 1] * c[k]),
            (m[k + 1] - m[k]) / denom
        );
    }

    return value;
}

double pressure_at_mu(
    double mu,
    const std::vector<double>& m,
    const std::vector<double>& c,
    double n0,
    double p0,
    double dn
) {
    constexpr int n_steps = 180;
    const double dx = (mu - m.front()) / static_cast<double>(n_steps);

    if (dx <= 0.0) return p0;

    double value = 0.5 *
                   (density_at_mu(m.front(), m, c, n0, dn) +
                    density_at_mu(mu, m, c, n0, dn));

    for (int i = 1; i < n_steps; ++i) {
        value += density_at_mu(m.front() + i * dx, m, c, n0, dn);
    }

    return dx * value * 1000.0 + p0;
}

void append_row_to_extension(EOSSample& sample,
                             double e,
                             double p,
                             double n,
                             double mu,
                             double cs2) {
    sample.e_ext.push_back(e);
    sample.p_ext.push_back(p);
    sample.n_ext.push_back(n);
    sample.mu_ext.push_back(mu);
    sample.cs2_ext.push_back(cs2);
}

bool valid_extension_vectors(const EOSSample& sample) {
    const std::size_t n = sample.e_ext.size();

    if (sample.p_ext.size() != n || sample.n_ext.size() != n ||
        sample.mu_ext.size() != n || sample.cs2_ext.size() != n) {
        return false;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!is_finite(sample.e_ext[i]) || !is_finite(sample.p_ext[i]) ||
            !is_finite(sample.n_ext[i]) || !is_finite(sample.mu_ext[i]) ||
            !is_finite(sample.cs2_ext[i])) {
            return false;
        }

        if (sample.e_ext[i] <= 0.0 || sample.p_ext[i] < 0.0 ||
            sample.n_ext[i] <= 0.0) {
            return false;
        }
    }

    for (std::size_t i = 1; i < n; ++i) {
        if (!(sample.e_ext[i] > sample.e_ext[i - 1])) return false;
        if (sample.p_ext[i] < sample.p_ext[i - 1]) return false;
        if (!(sample.mu_ext[i] > sample.mu_ext[i - 1])) return false;
    }

    return true;
}

bool valid_eos_vectors(const EOSSample& sample) {
    const std::size_t n = sample.e.size();

    if (n < 2 || sample.p.size() != n || sample.n.size() != n ||
        sample.mu.size() != n || sample.cs2.size() != n) {
        return false;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!is_finite(sample.e[i]) || !is_finite(sample.p[i]) ||
            !is_finite(sample.n[i]) || !is_finite(sample.mu[i]) ||
            !is_finite(sample.cs2[i])) {
            return false;
        }

        if (sample.e[i] <= 0.0 || sample.p[i] < 0.0 || sample.n[i] <= 0.0) {
            return false;
        }
    }

    for (std::size_t i = 1; i < n; ++i) {
        // Energy density must be strictly increasing for TOV-side interpolation.
        if (!(sample.e[i] > sample.e[i - 1])) {
            return false;
        }

        // Pressure may be constant across a first-order phase transition,
        // but it must not decrease.
        if (sample.p[i] < sample.p[i - 1]) {
            return false;
        }
    }

    return true;
}


EOSSample make_candidate(
    const Config& cfg,
    const CrustTable& crust,
    const CEFTTable& ceft,
    const PQCDTable& pqcd,
    std::mt19937& rng
) {
    EOSSample sample;

    const int n_mu_points = cfg.n_mu + 1;
    const double cmax_sample = cfg.sample_cmax ? uniform(rng, 0.0, 1.0) : cfg.cmax;

    std::vector<double> m(n_mu_points, 0.0);
    std::vector<double> cseg(n_mu_points, 0.0);
    std::vector<double> n_boundary(n_mu_points, 0.0);
    std::vector<double> p_boundary(n_mu_points, 0.0);

    const double eCrust = crust.e.back();
    const double pCrust = crust.p.back();
    const double nCrust = crust.n_over_ns.back();

    n_boundary[0] = cfg.n_cet;
    p_boundary[0] = uniform(rng, cfg.p_cet_soft, cfg.p_cet_stiff);

    const double gamma = std::log(pCrust / p_boundary[0]) / std::log(nCrust / n_boundary[0]);
    m[0] =
        (p_boundary[0]
         + ((eCrust - (pCrust / (gamma - 1.0))) * (n_boundary[0] / nCrust)
            + (p_boundary[0] / (gamma - 1.0)))) /
        (n_boundary[0] * ns) / 1000.0;

    cseg[0] = gamma / (1.0 + ((1000.0 * m[0] * n_boundary[0] * ns - p_boundary[0]) / p_boundary[0]));

    n_boundary[cfg.n_mu] = cfg.n_qcd;
    m[cfg.n_mu] = cfg.mu_qcd;

    const double X = std::exp(std::log(Xmin) + uniform(rng, 0.0, 1.0) * (std::log(Xmax) - std::log(Xmin)));
    cseg[cfg.n_mu] = cs2_qcd(cfg.mu_qcd, X);

    double mu_max = m[cfg.n_mu];
    if (cfg.pqcd_cs2_condition) {
        mu_max = cfg.mu_cs2_cond;
    }

    for (int i = 1; i < cfg.n_mu; ++i) {
        cseg[i] = uniform(rng, cfg.cmin, cmax_sample);

        if (cfg.pqcd_cs2_condition && i >= cfg.n_mu - cfg.n_mu_pqcd) {
            if (pqcd.mu_b.empty()) return sample;

            int index = static_cast<int>(
                static_cast<double>(i - cfg.n_mu + cfg.n_mu_pqcd) /
                static_cast<double>(cfg.n_mu_pqcd) *
                static_cast<double>(cfg.n_pqcd_file)
            );
            index = std::clamp(index, 0, cfg.n_pqcd_file - 1);

            std::normal_distribution<double> normal(
                pqcd.cs2_mean[index],
                2.0 * pqcd.cs2_std[index]
            );

            m[i] = pqcd.mu_b[index];
            cseg[i] = normal(rng);
        }
    }

    const int n_pt_points = cfg.include_phase_transition ? 2 : 0;
    for (int i = 1 + n_pt_points; i < cfg.n_mu - cfg.n_mu_pqcd; ++i) {
        m[i] = m[0] + uniform(rng, 0.0, 1.0) * (mu_max - m[0]);
    }

    const double muPT = m[0] + uniform(rng, 0.0, 1.0) * (mu_max - m[0]);
    if (cfg.include_phase_transition) {
        m[1] = muPT;
        m[2] = muPT;
    }

    std::sort(m.begin(), m.end());

    int PTpoint = 0;
    bool found_pt = false;
    if (cfg.include_phase_transition) {
        for (int i = 0; i < n_mu_points; ++i) {
            if (m[i] == muPT && !found_pt) {
                PTpoint = i;
                found_pt = true;
            }
        }
        if (!found_pt || PTpoint + 1 >= n_mu_points) return sample;
    }

    double dn = uniform(rng, cfg.dn_min_nsat, cfg.dn_max_nsat);
    if (!cfg.include_phase_transition) {
        dn = 0.0;
    }

    if (cfg.include_phase_transition) {
        const double nPTl_check = density_at_mu(muPT, m, cseg, n_boundary[0], dn) - dn * ns;
        const double nPTmin = cfg.n_pt_min_nsat * ns;
        const double nPTmax = cfg.n_pt_max_nsat * ns;

        if (!is_finite(nPTl_check) || nPTl_check < nPTmin || nPTl_check > nPTmax) {
            return sample;
        }
    }

    const double cs2PTl = cfg.include_phase_transition ? cseg[PTpoint] : 0.0;
    const double cs2PTr = cfg.include_phase_transition ? cseg[PTpoint + 1] : 0.0;

    bool pqcd_ok = false;
    const double sampled_n_qcd = density_at_mu(m[cfg.n_mu], m, cseg, n_boundary[0], dn);
    const double sampled_p_qcd = pressure_at_mu(m[cfg.n_mu], m, cseg, n_boundary[0], p_boundary[0], dn);
    const double target_n_qcd = baryon_density_qcd(cfg.mu_qcd, X);
    const double target_p_qcd = pressure_qcd(cfg.mu_qcd, X);

    if (sampled_n_qcd <= target_n_qcd + cfg.dn_qcd &&
        sampled_n_qcd >= target_n_qcd - cfg.dn_qcd &&
        sampled_p_qcd <= target_p_qcd + cfg.dp_qcd &&
        sampled_p_qcd >= target_p_qcd - cfg.dp_qcd) {
        pqcd_ok = true;
    }

    if (!cfg.impose_pqcd) {
        pqcd_ok = true;
    }

    if (!pqcd_ok) {
        return sample;
    }

    // EOS-side pQCD likelihood. The original code stored 1 for accepted samples.
    sample.pPQCD = 1.0;

    struct Row {
        double e = 0.0;
        double p = 0.0;
        double n = 0.0;
        double mu = 0.0;
        double cs2 = 0.0;
    };

    std::vector<Row> rows;
    rows.reserve(static_cast<std::size_t>(cfg.len_eos()));

    int pt_left_row_index = -1;

    for (int i = 0; i < cfg.n_crust; ++i) {
        rows.push_back({
            crust.e[i],
            crust.p[i],
            crust.n_over_ns[i] * ns,
            crust.mu[i],
            crust.cs2[i]
        });
    }

    double pCET_accum = 0.0;
    for (int i = 0; i < cfg.n_poly; ++i) {
        const double n_ratio = nCrust + (i + 1) * ((n_boundary[0] - nCrust) / (cfg.n_poly + 1.0));
        const double e_poly =
            (eCrust - (pCrust / (gamma - 1.0))) * (n_ratio / nCrust)
            + (pCrust * std::pow(n_ratio / nCrust, gamma)) / (gamma - 1.0);
        const double p_poly = pCrust * std::pow(n_ratio / nCrust, gamma);
        const double n_poly = n_ratio * ns;
        const double mu_poly = (e_poly + p_poly) / n_poly;
        const double cs2_poly =
            gamma /
            (1.0 +
             (((eCrust - (pCrust / (gamma - 1.0))) * (n_ratio / nCrust)
               + (pCrust * std::pow(n_ratio / nCrust, gamma)) / (gamma - 1.0)) /
              (pCrust * std::pow(n_ratio / nCrust, gamma))));

        const double med = linear_interp(n_poly, ceft.n, ceft.median_p);
        const double sig = linear_interp(n_poly, ceft.n, ceft.sigma_p);
        if (sig > 0.0) {
            pCET_accum += std::exp(-std::pow(p_poly - med, 2.0) / (2.0 * std::pow(sig, 2.0)));
        }

        rows.push_back({e_poly, p_poly, n_poly, mu_poly, cs2_poly});
    }

    sample.pCET = cfg.n_poly > 0 ? pCET_accum / static_cast<double>(cfg.n_poly) : 1.0;

    const double dmu = (m[cfg.n_mu] - m[0]) / static_cast<double>(cfg.n_cs2 - 1);
    int icrit = cfg.include_phase_transition ? static_cast<int>(std::floor((muPT - m[0]) / dmu)) : -1000;
    icrit = std::clamp(icrit, 0, std::max(0, cfg.n_cs2 - 2));

    for (int i = 0; i < cfg.n_cs2; ++i) {
        double mNow = m[0] + i * (m[cfg.n_mu] - m[0]) / static_cast<double>(cfg.n_cs2 - 1);
        if (cfg.include_phase_transition && (i == icrit || i == icrit + 1)) {
            mNow = muPT;
        }

        const double p_now = pressure_at_mu(mNow, m, cseg, n_boundary[0], p_boundary[0], dn);
        const double n_now_base = density_at_mu(mNow, m, cseg, n_boundary[0], dn);
        double e_now = -p_now + mNow * 1000.0 * n_now_base;
        double n_now = n_now_base;
        double cs2_now = cs_at_mu(mNow, m, cseg);
        const double mu_now = (e_now + p_now) / n_now;

        if (cfg.include_phase_transition && i == icrit) {
            e_now -= mNow * 1000.0 * dn * ns;
            n_now -= dn * ns;
            cs2_now = cs2PTl;

            sample.ePTl = e_now;
            sample.nPTl = n_now;
            sample.pPT = p_now;
        }

        if (cfg.include_phase_transition && i == icrit + 1) {
            cs2_now = cs2PTr;
        }

        rows.push_back({e_now, p_now, n_now, mu_now, cs2_now});

        if (cfg.include_phase_transition && i == icrit) {
            pt_left_row_index = static_cast<int>(rows.size()) - 1;
        }
    }
    
    if (cfg.include_phase_transition && pt_left_row_index > 0) {
        const double mu_left = rows[static_cast<std::size_t>(pt_left_row_index - 1)].mu / 1000.0;
        const double cs2_left = rows[static_cast<std::size_t>(pt_left_row_index - 1)].cs2;
        const double dmu_left = muPT - mu_left;

        if (dmu_left > 0.0 && cs2PTl > 0.0 && cs2PTl < 1.0) {
            const double slope = (cs2PTl - cs2_left) / dmu_left;
            double mu_stop = cfg.mu_qcd;

            if (slope > 0.0) {
                const double mu_causal = muPT + (1.0 - cs2PTl) / slope;
                if (mu_causal > muPT) {
                    mu_stop = std::min(mu_stop, mu_causal);
                }
            } else if (slope < 0.0) {
                const double mu_zero = muPT - cs2PTl / slope;
                if (mu_zero > muPT) {
                    mu_stop = std::min(mu_stop, mu_zero);
                }
            }

            if (mu_stop > muPT) {
                for (int i = 0; i <= pt_left_row_index; ++i) {
                    append_row_to_extension(
                        sample,
                        rows[static_cast<std::size_t>(i)].e,
                        rows[static_cast<std::size_t>(i)].p,
                        rows[static_cast<std::size_t>(i)].n,
                        rows[static_cast<std::size_t>(i)].mu,
                        rows[static_cast<std::size_t>(i)].cs2
                    );
                }

                double mu_prev = muPT;
                double n_prev = sample.nPTl;
                double p_prev = sample.pPT;
                const double step = std::max(dmu, 1.0e-6);
                double mu_next = muPT + step;

                while (mu_next < mu_stop) {
                    const double cs2_prev = cs2PTl + slope * (mu_prev - muPT);
                    const double cs2_next = cs2PTl + slope * (mu_next - muPT);

                    if (cs2_prev <= 0.0 || cs2_next <= 0.0 ||
                        cs2_prev > 1.0 || cs2_next > 1.0) {
                        break;
                    }

                    const double integrand_prev = 1.0 / (mu_prev * cs2_prev);
                    const double integrand_next = 1.0 / (mu_next * cs2_next);
                    const double integral = 0.5 * (integrand_prev + integrand_next) * (mu_next - mu_prev);

                    const double n_next = n_prev * std::exp(integral);
                    const double p_next =p_prev + 0.5 * (n_prev + n_next) * (mu_next - mu_prev) * 1000.0;
                    const double e_next = -p_next + mu_next * 1000.0 * n_next;

                    append_row_to_extension(
                        sample,
                        e_next,
                        p_next,
                        n_next,
                        mu_next * 1000.0,
                        cs2_next
                    );

                    mu_prev = mu_next;
                    n_prev = n_next;
                    p_prev = p_next;
                    mu_next += step;
                }

                if (mu_prev < mu_stop) {
                    const double cs2_prev = cs2PTl + slope * (mu_prev - muPT);
                    const double cs2_next = cs2PTl + slope * (mu_stop - muPT);

                    if (cs2_prev > 0.0 && cs2_prev <= 1.0 &&
                        cs2_next > 0.0 && cs2_next <= 1.0) {
                        const double integrand_prev = 1.0 / (mu_prev * cs2_prev);
                        const double integrand_next = 1.0 / (mu_stop * cs2_next);
                        const double integral =0.5 * (integrand_prev + integrand_next) * (mu_stop - mu_prev);

                        const double n_next = n_prev * std::exp(integral);
                        const double p_next =p_prev + 0.5 * (n_prev + n_next) * (mu_stop - mu_prev) * 1000.0;
                        const double e_next = -p_next + mu_stop * 1000.0 * n_next;

                        append_row_to_extension(
                            sample,
                            e_next,
                            p_next,
                            n_next,
                            mu_stop * 1000.0,
                            cs2_next
                        );
                    }
                }

                if (!valid_extension_vectors(sample)) {
                    //std::cout << "EOSext failed validation: "
                    //<< "size=" << sample.e_ext.size()
                    //<< ", PTpoint=" << PTpoint
                    //<< ", pt_left_row_index=" << pt_left_row_index
                    //<< ", muPT=" << muPT
                    //<< '\n';
                    sample.e_ext.clear();
                    sample.p_ext.clear();
                    sample.n_ext.clear();
                    sample.mu_ext.clear();
                    sample.cs2_ext.clear();
                }
            }
        }
    }

    if (cfg.impose_pqcd && cfg.n_pqcd > 0) {
        for (int i = 0; i < cfg.n_pqcd; ++i) {
            const double mub = 2.63 + i * (5.0 - 2.63) / 3.0;
            const double p_q = pressure_qcd(mub, X);
            const double n_q = baryon_density_qcd(mub, X);
            const double e_q = -p_q + n_q * mub * 1000.0;
            const double mu_q = (e_q + p_q) / n_q;
            const double cs2_q = cs2_qcd(mub, X);
            rows.push_back({e_q, p_q, n_q, mu_q, cs2_q});
        }
    }

    sample.e.reserve(rows.size());
    sample.p.reserve(rows.size());
    sample.n.reserve(rows.size());
    sample.mu.reserve(rows.size());
    sample.cs2.reserve(rows.size());

    for (const Row& row : rows) {
        sample.e.push_back(row.e);
        sample.p.push_back(row.p);
        sample.n.push_back(row.n);
        sample.mu.push_back(row.mu);
        sample.cs2.push_back(row.cs2);
    }

    sample.dn = dn;
    sample.muPT = cfg.include_phase_transition ? muPT : 0.0;
    sample.cs2PTr = cs2PTr;
    sample.cs2PTl = cs2PTl;

    if (!valid_eos_vectors(sample)) {
        sample = EOSSample{};
    }

    return sample;
}

void create_and_write_sample(hid_t file_id, int accepted, const EOSSample& sample) {
    hid_t group_id = hdf5::create_group(file_id, accepted);
    hid_t eos_group = hdf5::create_subgroup(group_id, std::string(schema::eos_group));
    hid_t tov_group = hdf5::create_subgroup(group_id, std::string(schema::tov_group));
    hid_t params_group = hdf5::create_subgroup(group_id, std::string(schema::params_group));

    const hsize_t len = static_cast<hsize_t>(sample.e.size());

    hdf5::create_dataset_double(eos_group, "e", len);
    hdf5::create_dataset_double(eos_group, "p", len);
    hdf5::create_dataset_double(eos_group, "n", len);
    hdf5::create_dataset_double(eos_group, "mu", len);
    hdf5::create_dataset_double(eos_group, "cs2", len);

    hdf5::write_double_array(eos_group, "e", sample.e);
    hdf5::write_double_array(eos_group, "p", sample.p);
    hdf5::write_double_array(eos_group, "n", sample.n);
    hdf5::write_double_array(eos_group, "mu", sample.mu);
    hdf5::write_double_array(eos_group, "cs2", sample.cs2);
    
    hid_t eosext_group = -1;
    if (!sample.e_ext.empty()) {
        eosext_group = hdf5::create_subgroup(group_id, "EOSext");
        const hsize_t len_ext = static_cast<hsize_t>(sample.e_ext.size());

        hdf5::create_dataset_double(eosext_group, "e", len_ext);
        hdf5::create_dataset_double(eosext_group, "p", len_ext);
        hdf5::create_dataset_double(eosext_group, "n", len_ext);
        hdf5::create_dataset_double(eosext_group, "mu", len_ext);
        hdf5::create_dataset_double(eosext_group, "cs2", len_ext);

        hdf5::write_double_array(eosext_group, "e", sample.e_ext);
        hdf5::write_double_array(eosext_group, "p", sample.p_ext);
        hdf5::write_double_array(eosext_group, "n", sample.n_ext);
        hdf5::write_double_array(eosext_group, "mu", sample.mu_ext);
        hdf5::write_double_array(eosext_group, "cs2", sample.cs2_ext);
    }

    for (int i = 0; i < 14; ++i) {
        const std::string dataset_name(schema::params_datasets[static_cast<std::size_t>(i)]);
        if (i == 12) {
            hdf5::create_dataset_int(params_group, dataset_name, 1);
        } else if (i == 13) {
            hdf5::create_dataset_double(params_group, dataset_name, 15);
        } else {
            hdf5::create_dataset_double(params_group, dataset_name, 1);
        }
    }

    hdf5::write_double_scalar(params_group, "dn", sample.dn);
    hdf5::write_double_scalar(params_group, "ePTl", sample.ePTl);
    hdf5::write_double_scalar(params_group, "pPT", sample.pPT);
    hdf5::write_double_scalar(params_group, "nPTl", sample.nPTl);
    hdf5::write_double_scalar(params_group, "muPT", sample.muPT);
    hdf5::write_double_scalar(params_group, "cs2PTr", sample.cs2PTr);
    hdf5::write_double_scalar(params_group, "cs2PTl", sample.cs2PTl);
    hdf5::write_double_scalar(params_group, "pCET", sample.pCET);
    hdf5::write_double_scalar(params_group, "pPQCD", sample.pPQCD);

    // TOV-side placeholders. These datasets exist so the schema is complete, but
    // the TOV executable should overwrite pM, pGW, ptot, nbranches, and pXray.
    hdf5::write_double_scalar(params_group, "pM", 1.0);
    hdf5::write_double_scalar(params_group, "pGW", 1.0);
    hdf5::write_double_scalar(params_group, "ptot", sample.pCET * sample.pPQCD);
    hdf5::write_int_scalar(params_group, "nbranches", 0);
    hdf5::write_double_array(params_group, "pXray", std::vector<double>(15, 0.0));

    H5Gclose(params_group);
    H5Gclose(tov_group);
    if (eosext_group >= 0) H5Gclose(eosext_group);
    H5Gclose(eos_group);
    H5Gclose(group_id);
}

} // namespace

int run_eos_sampler(const Config& cfg, unsigned int seed) {
    int mpi_was_initialized = 0;
    MPI_Initialized(&mpi_was_initialized);

    if (!mpi_was_initialized) {
        MPI_Init(nullptr, nullptr);
    }

    int rank = 0;
    int size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try {
        if (cfg.n_mu < 4) {
            throw std::runtime_error("n_mu must be at least 4 for the current EOS sampler");
        }
        if (cfg.n_mu_pqcd < 0 || cfg.n_mu_pqcd >= cfg.n_mu) {
            throw std::runtime_error("n_mu_pqcd must satisfy 0 <= n_mu_pqcd < n_mu");
        }
        if (cfg.n_cs2 < 2) {
            throw std::runtime_error("n_cs2 must be at least 2");
        }
        if (cfg.n_pqcd < 0) {
            throw std::runtime_error("n_pqcd must be non-negative");
        }

        const unsigned int rank_seed = seed + static_cast<unsigned int>(rank);

        if (cfg.verbose) {
            std::cout << "rank " << rank << ", seed " << rank_seed << '\n';
        }

        std::mt19937 rng(rank_seed);

        const CrustTable crust = load_crust(cfg);
        const CEFTTable ceft = load_ceft(cfg);

        PQCDTable pqcd;
        if (cfg.pqcd_cs2_condition) {
            pqcd = load_pqcd(cfg);
        }

	hid_t file_id = hdf5::create_file(cfg.output_prefix,rank,cfg.file_start_index);

        int accepted = 0;
        int attempts = 0;
        const int max_attempts = std::max(1000, cfg.n_eos * 100000);

        while (accepted < cfg.n_eos) {
            if (++attempts > max_attempts) {
                throw std::runtime_error(
                    "EOS sampler exceeded maximum attempts. Relax constraints or check input tables."
                );
            }

            EOSSample sample = make_candidate(cfg, crust, ceft, pqcd, rng);

            if (sample.e.empty()) {
                continue;
            }

            create_and_write_sample(file_id, accepted, sample);

            ++accepted;

	    if (cfg.verbose) {
    		std::cout << "rank " << rank
              		  << ", seed " << rank_seed
               		  << ", output file index " << (cfg.file_start_index + rank)
			  << '\n';
	    }
        }

        H5Fclose(file_id);

        int mpi_finalized = 0;
        MPI_Finalized(&mpi_finalized);
        if (!mpi_was_initialized && !mpi_finalized) {
            MPI_Finalize();
        }

        return 0;
    } catch (...) {
        int mpi_finalized = 0;
        MPI_Finalized(&mpi_finalized);
        if (!mpi_finalized) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        throw;
    }
}

} // namespace eos
