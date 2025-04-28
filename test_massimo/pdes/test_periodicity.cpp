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
        {0, 0, 1, 1},
        {0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1}
    }};

    std::array<int,2> order = {1, 2};
    std::array<bool, 2> periodic_dims_ = {false, true};

    MdArray<double, full_dynamic_extent_t<2>> weights(2, 6);
    weights.set_constant(1.0);

    auto basis = NurbsBasis<2>(knots, weights, order, periodic_dims_);

    std::cout << "Basis size: " << basis.size() << std::endl;

    const int steps = 4;
    std::cout << "\nChecking periodicity along second dimension (dim = 1):\n";
    for (int i = 0; i < steps; ++i) {
        double u = double(i) / (steps - 1); // from 0 to 1

        Eigen::Vector2d param_at_start(u, 0.0);
        Eigen::Vector2d param_at_end(u, 1.0);

        std::cout << "u = " << std::fixed << std::setprecision(2) << u << '\n';
        for (int j = 0; j < basis.size(); ++j) {
            auto bfun = basis[j];
            double val_start = bfun.derive(1)(param_at_start);
            double val_end = bfun.derive(1)(param_at_end);
            double diff = std::abs(val_start - val_end);
            if(diff !=0)
            std::cout << "  Basis[" << j << "] at v=0.0: " << val_start
                      << " | at v=1.0: " << val_end
                      << " | diff: " << diff << '\n';
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