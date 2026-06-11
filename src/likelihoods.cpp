#include "likelihoods.hpp"

#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eos {
namespace {

constexpr int kXrayOutputLength = 15;
constexpr double kPi = 3.141592653589793238462643383279502884;

struct Grid2D {
    int nx = 0;
    int ny = 0;
    std::vector<double> x;
    std::vector<double> y;

    // Stored in x-major order, matching GSL's common convention:
    //   w[ix * ny + iy]
    std::vector<double> w;
};

struct LikelihoodData {
    Grid2D gw;
    std::vector<Grid2D> xray;
};

bool is_finite(double x) {
    return std::isfinite(x);
}

void check_hdf5_id(hid_t id, const std::string& what) {
    if (id < 0) {
        throw std::runtime_error("HDF5 error while opening " + what);
    }
}

std::vector<double> read_hdf5_vector(hid_t file_id, const std::string& dataset_name) {
    hid_t dset = H5Dopen2(file_id, dataset_name.c_str(), H5P_DEFAULT);
    check_hdf5_id(dset, "dataset '" + dataset_name + "'");

    hid_t space = H5Dget_space(dset);
    check_hdf5_id(space, "dataspace for '" + dataset_name + "'");

    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 1) {
        H5Sclose(space);
        H5Dclose(dset);
        throw std::runtime_error("Dataset is not 1D: " + dataset_name);
    }

    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);

    std::vector<double> data(static_cast<std::size_t>(dim));
    if (!data.empty()) {
        const herr_t status = H5Dread(
            dset,
            H5T_NATIVE_DOUBLE,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            data.data()
        );

        if (status < 0) {
            H5Sclose(space);
            H5Dclose(dset);
            throw std::runtime_error("Failed to read dataset: " + dataset_name);
        }
    }

    H5Sclose(space);
    H5Dclose(dset);

    return data;
}

struct RawMatrix {
    hsize_t rows = 0;
    hsize_t cols = 0;
    std::vector<double> values;
};

RawMatrix read_hdf5_matrix(hid_t file_id, const std::string& dataset_name) {
    hid_t dset = H5Dopen2(file_id, dataset_name.c_str(), H5P_DEFAULT);
    check_hdf5_id(dset, "dataset '" + dataset_name + "'");

    hid_t space = H5Dget_space(dset);
    check_hdf5_id(space, "dataspace for '" + dataset_name + "'");

    const int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims != 2) {
        H5Sclose(space);
        H5Dclose(dset);
        throw std::runtime_error("Dataset is not 2D: " + dataset_name);
    }

    hsize_t dims[2] = {0, 0};
    H5Sget_simple_extent_dims(space, dims, nullptr);

    RawMatrix matrix;
    matrix.rows = dims[0];
    matrix.cols = dims[1];
    matrix.values.resize(static_cast<std::size_t>(dims[0] * dims[1]));

    if (!matrix.values.empty()) {
        const herr_t status = H5Dread(
            dset,
            H5T_NATIVE_DOUBLE,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            matrix.values.data()
        );

        if (status < 0) {
            H5Sclose(space);
            H5Dclose(dset);
            throw std::runtime_error("Failed to read dataset: " + dataset_name);
        }
    }

    H5Sclose(space);
    H5Dclose(dset);

    return matrix;
}

Grid2D read_grid2d(
    const std::string& file_path,
    const std::string& weight_dataset,
    const std::string& x_dataset,
    const std::string& y_dataset
) {
    hid_t file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    check_hdf5_id(file_id, "likelihood file '" + file_path + "'");

    try {
        Grid2D grid;
        grid.x = read_hdf5_vector(file_id, x_dataset);
        grid.y = read_hdf5_vector(file_id, y_dataset);

        if (grid.x.size() < 2 || grid.y.size() < 2) {
            throw std::runtime_error("Likelihood grid axes are too short in file: " + file_path);
        }

        grid.nx = static_cast<int>(grid.x.size());
        grid.ny = static_cast<int>(grid.y.size());

        const RawMatrix raw = read_hdf5_matrix(file_id, weight_dataset);
        const auto nx = static_cast<hsize_t>(grid.nx);
        const auto ny = static_cast<hsize_t>(grid.ny);

        if (raw.rows * raw.cols != nx * ny) {
            throw std::runtime_error("Likelihood grid size mismatch in file: " + file_path);
        }

        grid.w.assign(static_cast<std::size_t>(nx * ny), 0.0);

        if (raw.rows == ny && raw.cols == nx) {
            // HDF5 array is in row-major [y][x]. Convert to x-major [x][y].
            for (hsize_t iy = 0; iy < ny; ++iy) {
                for (hsize_t ix = 0; ix < nx; ++ix) {
                    grid.w[static_cast<std::size_t>(ix * ny + iy)] =
                        raw.values[static_cast<std::size_t>(iy * nx + ix)];
                }
            }
        } else if (raw.rows == nx && raw.cols == ny) {
            // Already x-major [x][y] in the file.
            grid.w = raw.values;
        } else {
            // Square grids are common. If the orientation cannot be inferred
            // from dimensions, use the old-code-compatible transpose.
            for (hsize_t iy = 0; iy < ny; ++iy) {
                for (hsize_t ix = 0; ix < nx; ++ix) {
                    const std::size_t src = static_cast<std::size_t>(iy * nx + ix);
                    const std::size_t dst = static_cast<std::size_t>(ix * ny + iy);
                    if (src < raw.values.size() && dst < grid.w.size()) {
                        grid.w[dst] = raw.values[src];
                    }
                }
            }
        }

        H5Fclose(file_id);
        return grid;
    } catch (...) {
        H5Fclose(file_id);
        throw;
    }
}

std::vector<double> linspace(double a, double b, int n) {
    if (n <= 0) return {};

    std::vector<double> out(static_cast<std::size_t>(n));
    if (n == 1 || a == b) {
        std::fill(out.begin(), out.end(), a);
        return out;
    }

    for (int i = 0; i < n; ++i) {
        out[static_cast<std::size_t>(i)] =
            a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return out;
}

double trapezoid(const std::vector<double>& f, const std::vector<double>& x) {
    if (f.size() != x.size() || f.size() < 2) return 0.0;

    double result = 0.0;
    for (std::size_t i = 1; i < f.size(); ++i) {
        result += 0.5 * (f[i - 1] + f[i]) * (x[i] - x[i - 1]);
    }
    return result;
}

double gaussian(double x, double mu, double sigma) {
    if (sigma <= 0.0) return 0.0;
    const double z = (x - mu) / sigma;
    return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * kPi));
}

double linear_interp_zero_outside(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double xq
) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    if (xq < x.front() || xq > x.back()) return 0.0;

    auto upper = std::upper_bound(x.begin(), x.end(), xq);
    std::size_t i = 0;
    if (upper == x.begin()) {
        i = 0;
    } else if (upper == x.end()) {
        i = x.size() - 2;
    } else {
        i = static_cast<std::size_t>(std::distance(x.begin(), upper) - 1);
    }

    const double dx = x[i + 1] - x[i];
    if (dx == 0.0) return 0.0;

    const double t = (xq - x[i]) / dx;
    return (1.0 - t) * y[i] + t * y[i + 1];
}

double bilinear_interp_zero_outside(const Grid2D& grid, double xq, double yq) {
    if (grid.nx < 2 || grid.ny < 2) return 0.0;
    if (grid.x.size() != static_cast<std::size_t>(grid.nx)) return 0.0;
    if (grid.y.size() != static_cast<std::size_t>(grid.ny)) return 0.0;
    if (grid.w.size() != static_cast<std::size_t>(grid.nx * grid.ny)) return 0.0;

    if (xq < grid.x.front() || xq > grid.x.back()) return 0.0;
    if (yq < grid.y.front() || yq > grid.y.back()) return 0.0;

    auto x_upper = std::upper_bound(grid.x.begin(), grid.x.end(), xq);
    auto y_upper = std::upper_bound(grid.y.begin(), grid.y.end(), yq);

    int ix = 0;
    int iy = 0;

    if (x_upper == grid.x.end()) {
        ix = grid.nx - 2;
    } else {
        ix = std::max(0, static_cast<int>(std::distance(grid.x.begin(), x_upper)) - 1);
    }

    if (y_upper == grid.y.end()) {
        iy = grid.ny - 2;
    } else {
        iy = std::max(0, static_cast<int>(std::distance(grid.y.begin(), y_upper)) - 1);
    }

    if (ix < 0 || iy < 0 || ix + 1 >= grid.nx || iy + 1 >= grid.ny) return 0.0;

    const double x1 = grid.x[static_cast<std::size_t>(ix)];
    const double x2 = grid.x[static_cast<std::size_t>(ix + 1)];
    const double y1 = grid.y[static_cast<std::size_t>(iy)];
    const double y2 = grid.y[static_cast<std::size_t>(iy + 1)];

    const double dx = x2 - x1;
    const double dy = y2 - y1;
    if (dx == 0.0 || dy == 0.0) return 0.0;

    const double tx = (xq - x1) / dx;
    const double ty = (yq - y1) / dy;

    auto w = [&](int x_index, int y_index) -> double {
        return grid.w[static_cast<std::size_t>(x_index * grid.ny + y_index)];
    };

    const double w11 = w(ix, iy);
    const double w21 = w(ix + 1, iy);
    const double w12 = w(ix, iy + 1);
    const double w22 = w(ix + 1, iy + 1);

    const double value =
        (1.0 - tx) * (1.0 - ty) * w11
        + tx * (1.0 - ty) * w21
        + (1.0 - tx) * ty * w12
        + tx * ty * w22;

    return is_finite(value) ? value : 0.0;
}

double m2_from_chirp(double mc, double m1) {
    if (mc <= 0.0 || m1 <= 0.0) return 0.0;

    const double disc = 81.0 * std::pow(m1, 5.0) - 12.0 * std::pow(mc, 5.0);
    if (disc <= 0.0) return 0.0;

    const double a = 9.0 * std::pow(m1, 2.5) + std::sqrt(disc);
    if (a <= 0.0) return 0.0;

    const double value =
        std::pow(mc, 5.0 / 3.0)
        * (
            2.0 * std::pow(3.0, 1.0 / 3.0) * std::pow(mc, 5.0 / 3.0)
            + std::pow(2.0, 1.0 / 3.0) * std::pow(a, 2.0 / 3.0)
        )
        / (
            std::pow(6.0, 2.0 / 3.0)
            * std::pow(m1, 1.5)
            * std::pow(a, 1.0 / 3.0)
        );

    return is_finite(value) && value <= m1 ? value : 0.0;
}

std::vector<double> branch_values(
    const std::vector<double>& values,
    int start,
    int stop
) {
    std::vector<double> out;
    if (values.empty() || start < 0 || stop < start) return out;

    const int last = std::min(stop, static_cast<int>(values.size()) - 1);
    out.reserve(static_cast<std::size_t>(last - start + 1));

    for (int i = start; i <= last; ++i) {
        out.push_back(values[static_cast<std::size_t>(i)]);
    }
    return out;
}

double mass_likelihood(const Config& cfg, const BranchInfo& branches) {
    if (!cfg.enable_mass_likelihood) return 1.0;
    if (cfg.mass_likelihood_points < 2) return 0.0;

    double result = 0.0;

    for (std::size_t b = 0; b < branches.start_mass.size(); ++b) {
        if (branches.stop_mass[b] <= branches.start_mass[b]) continue;

        const std::vector<double> masses = linspace(
            branches.start_mass[b],
            branches.stop_mass[b],
            cfg.mass_likelihood_points
        );

        std::vector<double> values(masses.size(), 0.0);
        for (std::size_t i = 0; i < masses.size(); ++i) {
            values[i] = gaussian(masses[i], cfg.pulsar_mass_mu, cfg.pulsar_mass_sigma);
        }

        result += trapezoid(values, masses);
    }

    return is_finite(result) ? result : 0.0;
}

double gw_likelihood(
    const Config& cfg,
    const TOVSequence& tov,
    const BranchInfo& branches,
    const Grid2D& gw_grid
) {
    if (!cfg.enable_gw_likelihood) return 1.0;
    if (cfg.mass_likelihood_points < 2) return 0.0;

    double result = 0.0;

    for (std::size_t b1 = 0; b1 < branches.idx_start.size(); ++b1) {
        if (branches.stop_mass[b1] <= branches.start_mass[b1]) continue;

        const std::vector<double> m1_grid = linspace(
            branches.start_mass[b1],
            branches.stop_mass[b1],
            cfg.mass_likelihood_points
        );

        const std::vector<double> branch_m1 = branch_values(tov.mass, branches.idx_start[b1], branches.idx_stop[b1]);
        const std::vector<double> branch_l1 = branch_values(tov.lambda, branches.idx_start[b1], branches.idx_stop[b1]);

        if (branch_m1.size() < 2 || branch_l1.size() < 2) continue;

        std::vector<double> prob_total(m1_grid.size(), 0.0);

        for (std::size_t b2 = 0; b2 < branches.idx_start.size(); ++b2) {
            const std::vector<double> branch_m2 = branch_values(tov.mass, branches.idx_start[b2], branches.idx_stop[b2]);
            const std::vector<double> branch_l2 = branch_values(tov.lambda, branches.idx_start[b2], branches.idx_stop[b2]);

            if (branch_m2.size() < 2 || branch_l2.size() < 2) continue;

            for (std::size_t i = 0; i < m1_grid.size(); ++i) {
                const double m1 = m1_grid[i];
                const double m2 = m2_from_chirp(cfg.chirp_mass, m1);
                if (m2 <= 0.0) continue;

                const double lambda1 = linear_interp_zero_outside(branch_m1, branch_l1, m1);
                const double lambda2 = linear_interp_zero_outside(branch_m2, branch_l2, m2);

                if (lambda1 <= 0.0 || lambda2 <= 0.0) continue;
                if (lambda1 > 1600.0 || lambda2 > 1600.0) continue;

                prob_total[i] += bilinear_interp_zero_outside(gw_grid, lambda1, lambda2);
            }
        }

        result += std::abs(trapezoid(prob_total, m1_grid));
    }

    return is_finite(result) ? result : 0.0;
}

double xray_likelihood_one_grid(
    const Config& cfg,
    const TOVSequence& tov,
    const BranchInfo& branches,
    const Grid2D& grid
) {
    if (cfg.mass_likelihood_points < 2) return 0.0;

    double result = 0.0;

    for (std::size_t b = 0; b < branches.idx_start.size(); ++b) {
        if (branches.stop_mass[b] <= branches.start_mass[b]) continue;

        const std::vector<double> m_grid = linspace(
            branches.start_mass[b],
            branches.stop_mass[b],
            cfg.mass_likelihood_points
        );

        const std::vector<double> branch_m = branch_values(tov.mass, branches.idx_start[b], branches.idx_stop[b]);
        const std::vector<double> branch_r = branch_values(tov.radius, branches.idx_start[b], branches.idx_stop[b]);

        if (branch_m.size() < 2 || branch_r.size() < 2) continue;

        std::vector<double> prob(m_grid.size(), 0.0);
        for (std::size_t i = 0; i < m_grid.size(); ++i) {
            const double m = m_grid[i];
            const double r = linear_interp_zero_outside(branch_m, branch_r, m);
            if (r <= 0.0) continue;

            // All MR likelihood grids are loaded with x = radius and y = mass.
            prob[i] = bilinear_interp_zero_outside(grid, r, m);
        }

        result += trapezoid(prob, m_grid);
    }

    return is_finite(result) ? result : 0.0;
}

double xray_likelihood(
    const Config& cfg,
    const TOVSequence& tov,
    const BranchInfo& branches,
    const std::vector<Grid2D>& grids,
    std::vector<double>& each
) {
    each.clear();

    if (!cfg.enable_xray_likelihoods) {
        each.assign(kXrayOutputLength, 0.0);
        return 1.0;
    }

    double total = 1.0;
    each.reserve(std::max<std::size_t>(kXrayOutputLength, grids.size()));

    for (const Grid2D& grid : grids) {
        const double value = xray_likelihood_one_grid(cfg, tov, branches, grid);
        each.push_back(value);
        total *= value;
    }

    while (each.size() < kXrayOutputLength) {
        each.push_back(0.0);
    }

    if (each.size() > kXrayOutputLength) {
        each.resize(kXrayOutputLength);
    }

    return is_finite(total) ? total : 0.0;
}

LikelihoodData load_likelihood_data(const Config& cfg) {
    LikelihoodData data;

    if (cfg.enable_gw_likelihood) {
        data.gw = read_grid2d(
            cfg.gw170817_file,
            "data",
            "x",
            "y"
        );
    }

    if (cfg.enable_xray_likelihoods) {
        data.xray.reserve(kXrayOutputLength);

        // NKS 15, two sources.
        data.xray.push_back(read_grid2d(
            cfg.nks1724b_file,
            "mcarlo/data/weights",
            "mcarlo/xval",
            "mcarlo/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.nks1810b_file,
            "mcarlo/data/weights",
            "mcarlo/xval",
            "mcarlo/yval"
        ));

        // SHB 18, seven sources.
        data.xray.push_back(read_grid2d(
            cfg.shb6304_file,
            "rescaled/data/like",
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shb6397_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shbM28_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shbM30_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shbX7_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shbX5_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        data.xray.push_back(read_grid2d(
            cfg.shbwCen_file,
            cfg.shb_dataset,
            "rescaled/xval",
            "rescaled/yval"
        ));

        // NAT 17.
        data.xray.push_back(read_grid2d(
            cfg.nat17_file,
            "hist2_table/data/avgs",
            "hist2_table/xval",
            "hist2_table/yval"
        ));

        // SHS 18.
        data.xray.push_back(read_grid2d(
            cfg.shsM13_file,
            "rescaled_0/data/like",
            "rescaled_0/xval",
            "rescaled_0/yval"
        ));

        // NICER: x = radius bins, y = mass bins.
        data.xray.push_back(read_grid2d(
            cfg.nicer0437_file,
            "weights",
            "r_bins",
            "m_bins"
        ));

        data.xray.push_back(read_grid2d(
            cfg.nicer0740_file,
            "weights",
            "r_bins",
            "m_bins"
        ));

        data.xray.push_back(read_grid2d(
            cfg.nicer0030_file,
            "weights",
            "r_bins",
            "m_bins"
        ));

        data.xray.push_back(read_grid2d(
            cfg.nicer0614_file,
            "weights",
            "r_bins",
            "m_bins"
        ));
    }

    return data;
}

std::vector<double> computed_stability(const TOVSequence& tov) {
    if (tov.stability.size() == tov.mass.size() && !tov.stability.empty()) {
        return tov.stability;
    }

    std::vector<double> stability(tov.mass.size(), 0.0);
    if (tov.mass.empty()) return stability;

    stability[0] = 1.0;
    for (std::size_t k = 1; k + 1 < tov.mass.size(); ++k) {
        const bool mass_increases = tov.mass[k + 1] >= tov.mass[k];

        if (mass_increases && stability[k - 1] > 0.0) {
            stability[k] = stability[k - 1];
        } else if (mass_increases && stability[k - 1] < 0.0) {
            stability[k] = -stability[k - 1] + 1.0;
        } else if (!mass_increases && stability[k - 1] > 0.0) {
            stability[k] = -stability[k - 1];
        } else {
            stability[k] = stability[k - 1];
        }
    }

    if (tov.mass.size() > 1) {
        stability.back() = stability[tov.mass.size() - 2];
    }

    return stability;
}

} // namespace

BranchInfo find_branches(const TOVSequence& tov) {
    BranchInfo out;
    if (tov.mass.size() < 2) return out;

    const std::vector<double> stability = computed_stability(tov);

    int nbranches = 0;
    for (double s : stability) {
        if (s > 0.0) {
            nbranches = std::max(nbranches, static_cast<int>(std::round(s)));
        }
    }

    // Preserve the original code's convention: use at most two stable branches.
    nbranches = std::min(nbranches, 2);

    for (int branch = 1; branch <= nbranches; ++branch) {
        int start = -1;
        int stop = -1;

        for (std::size_t i = 0; i < stability.size(); ++i) {
            if (static_cast<int>(std::round(stability[i])) == branch) {
                start = static_cast<int>(i);
                break;
            }
        }

        for (int i = static_cast<int>(stability.size()) - 1; i >= 0; --i) {
            if (static_cast<int>(std::round(stability[static_cast<std::size_t>(i)])) == branch) {
                stop = i;
                break;
            }
        }

        if (start < 0 || stop <= start) continue;

        // For the first branch, start near 0.5 solar masses as in the old code.
        if (branch == 1) {
            for (int i = start; i <= stop; ++i) {
                if (tov.mass[static_cast<std::size_t>(i)] >= 0.5) {
                    start = std::max(0, i - 1);
                    break;
                }
            }
        }

        if (stop <= start) continue;

        out.idx_start.push_back(start);
        out.idx_stop.push_back(stop);
        out.start_mass.push_back(tov.mass[static_cast<std::size_t>(start)]);
        out.stop_mass.push_back(tov.mass[static_cast<std::size_t>(stop)]);
    }

    return out;
}

double tov_mass(const TOVSequence& tov) {
    double mtov = 0.0;

    if (tov.mass.empty()) return mtov;

    if (tov.stability.size() == tov.mass.size()) {
        for (std::size_t i = 0; i < tov.mass.size(); ++i) {
            if (tov.stability[i] > 0.0) {
                mtov = std::max(mtov, tov.mass[i]);
            }
        }
    } else {
        for (double mass : tov.mass) {
            mtov = std::max(mtov, mass);
        }
    }

    return mtov;
}

LikelihoodEvaluator::LikelihoodEvaluator(const Config& cfg)
    : cfg_(cfg) {}

LikelihoodResult LikelihoodEvaluator::evaluate(
    const TOVSequence& tov,
    const BranchInfo& branches
) const {
    // Loaded once per process. Construct one LikelihoodEvaluator per run/config.
    static const LikelihoodData data = load_likelihood_data(cfg_);

    LikelihoodResult out;

    out.p_mass = mass_likelihood(cfg_, branches);

    if (cfg_.enable_gw_likelihood) {
        out.p_gw = gw_likelihood(cfg_, tov, branches, data.gw);
    } else {
        out.p_gw = 1.0;
    }

    if (cfg_.enable_xray_likelihoods) {
        out.p_xray = xray_likelihood(cfg_, tov, branches, data.xray, out.p_xray_each);
    } else {
        out.p_xray = 1.0;
        out.p_xray_each.assign(kXrayOutputLength, 0.0);
    }

    return out;
}

} // namespace eos
