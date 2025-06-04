#include "fields.h"

namespace adv_diff_sphere {

constexpr int M = 3;
using Vec = Eigen::Matrix<double, M, 1>;
using Fun = std::function<double(const Vec&)>;


// RHS: simple forcing f(p) = x
inline fdapde::ScalarField<M> make_rhs() {
    return fdapde::ScalarField<M>([](const Vec& p) {
        return 1;
    });
}

// Advection field: dual-rotation plus shear
inline fdapde::VectorField<M, M, Fun> make_b_field() {
    fdapde::VectorField<M, M, Fun> b;
    // First rotation: Ω1 = (0,0,1)
    // Second rotation: Ω2 = (1,0,0)
    // Plus meridional shear strength
    constexpr double osc = 4;
    b(0,0) = [](const Vec& p) {
        double x = p(0), y = p(1), z = p(2);
        double phi = std::atan2(y, x);
        double osc_factor = std::sin(osc * phi);
        return (-y) * osc_factor;
    };
    b(1,0) = [](const Vec& p) {
        double x = p(0), y = p(1), z = p(2);
        double phi = std::atan2(y, x);
        double osc_factor = std::sin(osc * phi);
        return (x) * osc_factor;
    };
    b(2,0) = [](const Vec& p) {
        double x = p(0), y = p(1), z = p(2);
        double phi = std::atan2(y, x);
        double osc_factor = std::sin(osc * phi);
        return (0.0) * osc_factor;
    };
    return b;
}

} // namespace adv_diff_sphere
