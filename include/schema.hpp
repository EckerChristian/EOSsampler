#pragma once

#include <array>
#include <string_view>

namespace eos::schema {

inline constexpr std::string_view eos_group = "EOS";
inline constexpr std::string_view tov_group = "TOV";
inline constexpr std::string_view params_group = "params";

inline constexpr std::array<std::string_view, 3> subgroups = {
    eos_group,
    tov_group,
    params_group
};

inline constexpr std::array<std::string_view, 5> eos_datasets = {
    "e",
    "p",
    "n",
    "mu",
    "cs2"
};

inline constexpr std::array<std::string_view, 9> tov_datasets = {
    "M",
    "R",
    "Lambda",
    "stab",
    "e_cent",
    "mu_cent",
    "n_cent",
    "cs2_cent",
    "Rqm"
};

inline constexpr std::array<std::string_view, 14> params_datasets = {
    "dn",
    "ePTl",
    "pPT",
    "nPTl",
    "muPT",
    "cs2PTr",
    "cs2PTl",
    "pCET",
    "pPQCD",
    "pM",
    "pGW",
    "ptot",
    "nbranches",
    "pXray"
};

} // namespace eos::schema
