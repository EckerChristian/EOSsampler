#pragma once

#include "config.hpp"

#include <vector>

namespace eos {

struct TOVSequence {
    std::vector<double> mass;
    std::vector<double> radius;
    std::vector<double> lambda;
    std::vector<double> stability;
    std::vector<double> e_cent;
    std::vector<double> mu_cent;
    std::vector<double> n_cent;
    std::vector<double> cs2_cent;
    std::vector<double> r_qm;
};

struct BranchInfo {
    std::vector<int> idx_start;
    std::vector<int> idx_stop;
    std::vector<double> start_mass;
    std::vector<double> stop_mass;
};

struct LikelihoodResult {
    double p_mass = 1.0;
    double p_gw = 1.0;
    double p_xray = 1.0;
    std::vector<double> p_xray_each;
};

class LikelihoodEvaluator {
public:
    explicit LikelihoodEvaluator(const Config& cfg);

    LikelihoodResult evaluate(
        const TOVSequence& tov,
        const BranchInfo& branches
    ) const;

private:
    Config cfg_;
};

BranchInfo find_branches(const TOVSequence& tov);

double tov_mass(const TOVSequence& tov);

} // namespace eos
