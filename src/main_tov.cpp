#include "config.hpp"
#include "tov.hpp"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    try {
        const std::string config_path =
            argc >= 2 ? argv[1] : "config/default.yaml";

        eos::Config cfg = eos::load_config(config_path);
        return eos::run_tov_solver(cfg);
    } catch (const std::exception& e) {
        std::cerr << "TOV failed: " << e.what() << '\n';
        return 1;
    }
}
