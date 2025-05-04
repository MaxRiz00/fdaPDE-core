#include "isogeometric.h"

#include <unsupported/Eigen/SparseExtra>
#include <fstream>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>

using namespace fdapde;

int main() {

    std::array<std::vector<double>, 2> knots = {{
        {0, 0, 0,0.5, 1, 1, 1},
        {0, 0, 0,0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1}
    }};

    std::array<int,2> order = {2, 3};
    std::array<bool, 2> periodic_dims_ = {false, true};

    MdArray<double, full_dynamic_extent_t<2>> weights(4, 8);
    weights.set_constant(1.0);

    auto basis = NurbsBasis<2>(knots, weights, order, periodic_dims_);

    std::cout << "Basis size: " << basis.size() << std::endl;

    const int steps = 6;
    std::cout << "\nChecking periodicity along second dimension (dim = 1):\n";
    for (int i = 0; i < steps; ++i) {
        double u = double(i) / (steps - 1); // from 0 to 1

        Eigen::Vector2d param_at_start(u, 0.0);
        Eigen::Vector2d param_at_end(u, 1.0);

        std::cout << "u = " << std::fixed << std::setprecision(2) << u << '\n';
        for (int j = 0; j < basis.size(); ++j) {
            auto bfun = basis[j];
            auto hess_start = bfun.hessian(param_at_start);
            auto hess_end = bfun.hessian(param_at_end);
            for(int k = 0; k < 2; ++k) {
                for(int l = 0; l < 2; ++l) {
                    if(std::abs(hess_start(k, l) - hess_end(k, l)) > 1e-10) {
                        std::cout << "  Hessian[" << j << "] at v=0.0: " << hess_start(k, l)
                                  << " | at v=1.0: " << hess_end(k, l)
                                  << " | diff: " << std::abs(hess_start(k, l) - hess_end(k, l)) << '\n';
                    }
                }
            }
        }
    }
    std::cout << "\nPeriodic check completed.\n";
    /*
    std::cout << "\nChecking first derivatives along periodic direction (dim = 1):\n";
    for (int i = 0; i < steps; ++i) {
        double u = double(i) / (steps - 1);

        Eigen::Vector2d param_at_start(u, 0.0);
        Eigen::Vector2d param_at_end(u, 1.0);

        std::cout << "u = " << std::fixed << std::setprecision(2) << u << '\n';
        for (int j = 0; j < basis.size(); ++j) {
            auto bfun = basis[j];
            double deriv_start = bfun.derive(1)(param_at_start);
            double deriv_end = bfun.derive(1)(param_at_end);
            double diff = std::abs(deriv_start - deriv_end);
            std::cout << "  Basis[" << j << "]' at v=0.0: " << deriv_start
                      << " | at v=1.0: " << deriv_end
                      << " | diff: " << diff << '\n';
        }
    }
    std::cout << "\nPeriodic derivative check completed.\n";
    */
    

    return 0;
}