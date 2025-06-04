#include "isogeometric.h"
#include "exact_solution.h"
#include "../helpers.h"
#include "utils/utils.h"
int main() {
    using namespace fdapde;

    // Degree 3 with richer knot vector
    std::vector<double> knots = {
        0.0, 0.0, 0.0, 0.0,       // degree+1 left
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 7.0, 7.0, 7.0        // degree+1 right
    };

    int degree = 3;
    bool periodic = true;

    BSplineBasis basis(knots, degree, periodic);
    int n_basis = basis.n_basis();

    double x_start = 0.0;
    double x_end   = 7.0;
    int max_order = 2;

    auto derivs_start = basis.evaluate_der_basis(x_start, max_order, true);
    auto derivs_end   = basis.evaluate_der_basis(x_end,   max_order, true);

    // Print all derivatives at x_start
    std::cout << "Derivatives at x_start = " << x_start << ":\n";
    for (int k = 0; k <= max_order; ++k) {
        std::cout << "D^" << k << " B_i(x_start): ";
        for (int i = 0; i < n_basis; ++i) {
            std::cout << derivs_start[k][i] << " ";
        }
        std::cout << "\n";
    }

    // Print all derivatives at x_end
    std::cout << "\nDerivatives at x_end = " << x_end << ":\n";
    for (int k = 0; k <= max_order; ++k) {
        std::cout << "D^" << k << " B_i(x_end): ";
        for (int i = 0; i < n_basis; ++i) {
            std::cout << derivs_end[k][i] << " ";
        }
        std::cout << "\n";
    }

    // Periodicity check
    auto almost_equal = [](double a, double b, double tol = 1e-10) {
        return std::abs(a - b) < tol;
    };

    std::cout << "\nChecking periodic mapping for degree 3 (up to 2nd derivative):\n";
    std::cout << "Compare D^k B_i(x_start) ≈ D^k B_{(i + n_basis - 3) % n_basis}(x_end)\n";
    std::cout << "n_basis = " << n_basis << "\n\n";

    for (int k = 0; k <= max_order; ++k) {
        std::cout << "Derivative order " << k << ":\n";
        for (int i = 0; i < n_basis; ++i) {
            int j = (i + n_basis - degree) % n_basis;
            double v_start = derivs_start[k][i];
            double v_end   = derivs_end[k][j];
            bool match = almost_equal(v_start, v_end);

            std::cout << "D^" << k << " B_" << i << "(" << x_start << ") = " << v_start
                      << " ~ D^" << k << " B_" << j << "(" << x_end << ") = " << v_end
                      << " --> " << (match ? "OK" : "Mismatch") << "\n";
        }
        std::cout << "\n";
    }

    return 0;
}