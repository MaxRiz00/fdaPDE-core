#include "fields.h"


namespace advdiff_ring {

constexpr int M = 2;
using Vec = Eigen::Matrix<double, M, 1>;
using Fun = std::function<double(const Vec&)>;

// Exact solution
inline fdapde::ScalarField<M> make_u_exact() {
    return fdapde::ScalarField<M>([](const Vec& p) {
        double x = p(0), y = p(1);
        return std::sin(M_PI * x) * std::cos(M_PI * y);
    });
}

// RHS of the Poisson equation
inline fdapde::ScalarField<M> make_rhs() {
    return fdapde::ScalarField<M>([](const Vec& p) {
        double x = p(0), y = p(1);
        return 2 * M_PI * M_PI *std::sin(M_PI * x) * std::cos(M_PI * y);
    });
}

// Gradient of the exact solution
inline fdapde::VectorField<M, M, Fun> make_grad_u_exact() {
    fdapde::VectorField<M, M, Fun> df;

    df(0, 0) = [](const Vec& p) {
        double x = p(0), y = p(1);
        return M_PI * std::cos(M_PI * x) * std::cos(M_PI * y);
    };

    df(1, 0) = [](const Vec& p) {
        double x = p(0), y = p(1);
        return - M_PI * std::sin(M_PI * y) * std::sin(M_PI * x);
    };

    return df;
}


inline fdapde::ScalarField<M> make_g_neumann() {
    return fdapde::ScalarField<M>([](const Vec& p) {
        double x = p(0), y = p(1);
        constexpr double eps = 1e-8;
        if (std::abs(x - 0.0) < eps)
            return - M_PI *  std::cos(M_PI * y);
        else if (std::abs(x - 1.0) < eps)
            return - M_PI *  std::cos(M_PI * y);
        else if (std::abs(y - 0.0) < eps)
            return 0.0;
        else if (std::abs(y - 1.0) < eps)
            return 0.0;
        else
            return 0.0; // interior
    });
}

/*

inline fdapde::ScalarField<M> make_g_neumann() {
    return fdapde::ScalarField<M>([](const Vec& p) {
        double x = p(0), y = p(1);
        if(x == 0 || y == 0) {
            return 0.0; // Neumann condition at the origin
        } else {
            return 2.0; // Some constant value for other points
        }
    });
}
    */

} // namespace ring