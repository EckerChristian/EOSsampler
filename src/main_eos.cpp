#include "config.hpp"
#include "eos.hpp"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    try {
        const std::string config_path =
            argc >= 2 ? argv[1] : "config/default.yaml";

        eos::Config cfg = eos::load_config(config_path);

        unsigned int seed = 12345;
        if (argc >= 3) {
            seed = static_cast<unsigned int>(std::stoul(argv[2]));
        }

        return eos::run_eos_sampler(cfg, seed);
    } catch (const std::exception& e) {
        std::cerr << "EOS failed: " << e.what() << '\n';
        return 1;
    }
}
