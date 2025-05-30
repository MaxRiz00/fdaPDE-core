
#include "nurbs.h"

// Finite difference approximation of Hessian
template<int M>
Eigen::Matrix<double, M, M> finite_difference_hessian(fdapde::Nurbs<M> nurbs, const Eigen::Matrix<double, M, 1>& p, double eps = 1e-5) {
    Eigen::Matrix<double, M, M> H;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            Eigen::Matrix<double, M, 1> p_ijp = p, p_ijm = p, p_jp = p, p_jm = p;

            p_ijp(i) += eps; p_ijp(j) += eps;
            p_ijm(i) += eps; p_ijm(j) -= eps;
            p_jp(i) -= eps;  p_jp(j) += eps;
            p_jm(i) -= eps;  p_jm(j) -= eps;

            double f1 = nurbs(p_ijp);
            double f2 = nurbs(p_ijm);
            double f3 = nurbs(p_jp);
            double f4 = nurbs(p_jm);

            H(i,j) = (f1 - f2 - f3 + f4) / (4 * eps * eps);
        }
    }
    return H;
}


// Finite difference approximation of the gradient
template<int M>
Eigen::Matrix<double, M, 1> finite_difference_gradient(const fdapde::Nurbs<M>& nurbs, const Eigen::Matrix<double, M, 1>& p, double eps = 1e-6) {
    Eigen::Matrix<double, M, 1> grad;
    for (int i = 0; i < M; ++i) {
        Eigen::Matrix<double, M, 1> p_plus = p;
        Eigen::Matrix<double, M, 1> p_minus = p;
        p_plus(i) += eps;
        p_minus(i) -= eps;
        grad(i) = (nurbs(p_plus) - nurbs(p_minus)) / (2.0 * eps);
    }
    return grad;
}

void test_nurbs_hessian() {
    constexpr int M = 2;
    using namespace fdapde;

    // Example: uniform quadratic B-spline over [0, 1]
    std::vector<double> knots = {0, 0, 0, 0.5, 1, 1, 1};  // Clamped uniform B-spline
    MdArray<double, MdExtents<Dynamic,Dynamic>> weights(4,4);
    weights.set_constant(1.0); // Simple weight function (no rationality)
    weights(0,0) = 7.0;
    
    // Choose an internal basis function
    std::array<int, M> index = {1, 1};
    std::array<int, M> degree = {2, 2};
    std::array<std::vector<double>, M> knot_vectors = {knots, knots};

    NurbsBasis<M> nurbs(knot_vectors, weights, degree);

    Eigen::Matrix<double, M, 1> p;
    p << 0.1, 0.1;

    auto analytical_hessian = nurbs[0].hessian(p);
    auto numerical_hessian = finite_difference_hessian(nurbs[0], p);

    auto analytical_grad = nurbs[0].gradient(p);
    auto numerical_grad = finite_difference_gradient(nurbs[0], p);

    std::cout << "Analytical Hessian:\n" << analytical_hessian << "\n";
    std::cout << "Numerical Hessian:\n" << numerical_hessian << "\n";
    std::cout << "Error:\n" << (analytical_hessian - numerical_hessian).norm() << "\n";

    if ((analytical_hessian - numerical_hessian).norm() < 1e-5)
        std::cout << "PASS: Hessian is correct.\n";
    else
        std::cout << "FAIL: Hessian is not matching numerical approximation.\n";


    if((analytical_grad - numerical_grad).norm() < 1e-5)
        std::cout << "PASS: Gradient is correct.\n";
    else
        std::cout << "FAIL: Gradient is not matching numerical approximation.\n";
}

int main() {
    test_nurbs_hessian();
    return 0;
}