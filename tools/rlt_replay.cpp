#include "lp_backend.h"

#include <chrono>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc == 1) {
            Highs highs;
            std::cout << "{\"highs_runtime\":\"" << highs.version()
                      << "\",\"compiler\":\"" << __VERSION__ << "\"}\n";
            return 0;
        }
        if (argc != 4) {
            throw std::runtime_error("usage: rlt_replay model.mps options basis");
        }

        for (int repetition = 0; repetition < 3; ++repetition) {
            for (int order = 0; order < 2; ++order) {
                const bool warm = (order + repetition) % 2 != 0;
                Highs highs;
                configure_highs(highs);
                highs_require(highs.readOptions(argv[2]), "options");
                highs_require(highs.readModel(argv[1]), "model");
                if (warm) {
                    highs_require(highs.readBasis(argv[3]), "basis");
                    if (!highs.getBasis().valid) {
                        throw std::runtime_error("invalid replay basis");
                    }
                }

                const auto start = std::chrono::steady_clock::now();
                highs_require(highs.run(), "replay");
                const double seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start).count();
                if (!highs_solution_valid(highs)) {
                    throw std::runtime_error("invalid replay solution");
                }
                std::cout << std::setprecision(17) << "REPLAY {\"repetition\":"
                          << repetition << ",\"warm\":" << (warm ? "true" : "false")
                          << ",\"objective\":" << highs.getObjectiveValue()
                          << ",\"wall_s\":" << seconds
                          << ",\"iterations\":" << highs.getInfo().simplex_iteration_count
                          << "}\n";
            }
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
