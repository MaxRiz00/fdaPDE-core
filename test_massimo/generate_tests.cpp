#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "isogeometric.h"
#include "utils/utils.h"

using SpMatrix = Eigen::SparseMatrix<double>;
using namespace fdapde;
/*
int main() {
    std::vector<double> knots = {0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0};
    int degree = 2;
    std::vector<double> eval_points = {0, 0.25, 0.75, 1.5, 2.25, 3.0, 4.0};

    // Non-periodic basis
    {
        BSplineBasis basis(knots, degree, false);
        int n_basis = basis.n_basis();
        int n_points = eval_points.size();
        SpMatrix basis_out(n_basis + 1, n_points);  // +1 for x values

        for (int j = 0; j < n_points; ++j) {
            basis_out.insert(0, j) = eval_points[j];
            auto values = basis.evaluate_basis(eval_points[j]);
            for (int i = 0; i < values.size(); ++i) {
                basis_out.insert(i + 1, j) = values[i];
            }
        }
        Eigen::saveMarket(basis_out, "../data/sp_basis_nonperiodic.mtx");

        SpMatrix deriv_out(n_basis + 1, n_points);
        for (int j = 0; j < n_points; ++j) {
            deriv_out.insert(0, j) = eval_points[j];
            auto ders = basis.evaluate_der_basis(eval_points[j], 1);
            for (int i = 0; i < ders[1].size(); ++i) {
                deriv_out.insert(i + 1, j) = ders[1][i];
            }
        }
        Eigen::saveMarket(deriv_out, "../data/sp_basis_first_der_nonperiodic.mtx");
    }

    // Periodic basis
    {
        BSplineBasis basis(knots, degree, true);
        int n_basis = basis.n_basis() + 2;
        int n_points = eval_points.size();
        SpMatrix basis_out(n_basis + 1, n_points);  // +1 for x values

        for (int j = 0; j < n_points; ++j) {
            basis_out.insert(0, j) = eval_points[j];
            auto values = basis.evaluate_basis(eval_points[j]);
            // print values
            for (int i = 0; i < values.size(); ++i) {
                std::cout << values[i] << " ";
            }
            std::cout << std::endl;
            for (int i = 0; i < values.size(); ++i) {
                basis_out.insert(i + 1, j) = values[i];
            }
        }
        Eigen::saveMarket(basis_out, "../data/sp_basis_periodic.mtx");
    }

    return 0;
}

int main() {

    // 2D NURBS test setup
    std::array<std::vector<double>, 2> knots;
    knots[0] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0};
    knots[1] = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0};

    std::array<int, 2> degree = {2, 2};
    std::array<bool, 2> periodicity = {false, false};

    MdArray<double, full_dynamic_extent_t<2>> weights(4, 4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            weights(i, j) = 1.0 + 0.1 * i * j;

    NurbsBasis<2> basis(knots, weights, degree, periodicity);

    std::vector<std::array<double, 2>> eval_points = {
        {0,0},{0.25, 0.25}, {0.5, 0.5}, {0.75, 0.75}, {1.25, 1.25}, {1.75, 1.75},{2,2}
    };

    int n_basis = basis.size();
    int n_points = eval_points.size();

    SpMatrix out_eval(n_basis + 2, n_points);      // 2 for x/y
    SpMatrix out_grad(n_basis * 2 + 2, n_points);  // 2 gradients per basis
    SpMatrix out_hess(n_basis * 4 + 2, n_points);  // 4 hessian values per basis

    for (int j = 0; j < n_points; ++j) {
        const auto& pt = eval_points[j];
        Eigen::Vector2d p(pt[0], pt[1]);

        out_eval.insert(0, j) = pt[0];
        out_eval.insert(1, j) = pt[1];
        out_grad.insert(0, j) = pt[0];
        out_grad.insert(1, j) = pt[1];
        out_hess.insert(0, j) = pt[0];
        out_hess.insert(1, j) = pt[1];

        for (int i = 0; i < n_basis; ++i) {
            out_eval.insert(i + 2, j) = basis[i](p);

            Eigen::Vector2d grad = basis[i].gradient(p);
            out_grad.insert(2 + 2 * i, j)     = grad(0);
            out_grad.insert(2 + 2 * i + 1, j) = grad(1);

            Eigen::Matrix2d hess = basis[i].hessian(p);
            out_hess.insert(2 + 4 * i + 0, j) = hess(0, 0);
            out_hess.insert(2 + 4 * i + 1, j) = hess(0, 1);
            out_hess.insert(2 + 4 * i + 2, j) = hess(1, 0);
            out_hess.insert(2 + 4 * i + 3, j) = hess(1, 1);
        }
    }

    Eigen::saveMarket(out_eval, "../data/nurbs2d_eval.mtx");
    Eigen::saveMarket(out_grad, "../data/nurbs2d_grad.mtx");
    Eigen::saveMarket(out_hess, "../data/nurbs2d_hess.mtx");

    return 0;
}
*/


int main() {
    IsoMesh<2, 3> mesh = IsoMesh<2, 3>::torus();

    const int grid = 20;  // increase this for more points
    std::vector<Eigen::Vector2d> u_points;
    for (int i = 0; i <= grid; ++i) {
        for (int j = 0; j <= grid; ++j) {
            double u = static_cast<double>(i) / grid;
            double v = static_cast<double>(j) / grid;
            u_points.emplace_back(u, v);
        }
    }

    int n_points = u_points.size();
    SpMatrix eval_out(4, n_points);         // t, x, y, z
    SpMatrix jacobian_out(7, n_points);     // t, ∂x/∂u, ..., ∂z/∂v
    SpMatrix hessian_out(13, n_points);     // t, H[3][2][2]
    SpMatrix invert_in(4, n_points);        // t, x, y, z
    SpMatrix invert_out(3, n_points);       // t, u, v

    for (int j = 0; j < n_points; ++j) {
        double t = j;
        auto u = u_points[j];
        auto x = mesh.eval_param(u);
        auto derivs = mesh.eval_param_derivatives(u, true);
        double t1 = 0, t2 = 0;
        auto u_back = mesh.invert_point(x, t1, t2);

        eval_out.insert(0, j) = t;
        eval_out.insert(1, j) = x(0);
        eval_out.insert(2, j) = x(1);
        eval_out.insert(3, j) = x(2);

        jacobian_out.insert(0, j) = t;
        for (int d = 0; d < 3; ++d)
            for (int i = 0; i < 2; ++i)
                jacobian_out.insert(1 + d * 2 + i, j) = derivs.first_derivative(d, i);

        hessian_out.insert(0, j) = t;
        const auto& H = *(derivs.second_derivative);
        for (int d = 0; d < 3; ++d)
            for (int i = 0; i < 2; ++i)
                for (int k = 0; k < 2; ++k) {
                    int row = 1 + d * 4 + i * 2 + k;
                    hessian_out.insert(row, j) = H(d, i, k);
                }

        invert_in.insert(0, j) = t;
        invert_in.insert(1, j) = x(0);
        invert_in.insert(2, j) = x(1);
        invert_in.insert(3, j) = x(2);

        invert_out.insert(0, j) = t;
        invert_out.insert(1, j) = u_back(0);
        invert_out.insert(2, j) = u_back(1);
    }

    Eigen::saveMarket(eval_out, "../data/torus_eval_param.mtx");
    Eigen::saveMarket(jacobian_out, "../data/torus_eval_param_jacobian.mtx");
    Eigen::saveMarket(hessian_out, "../data/torus_eval_param_hessian.mtx");
    Eigen::saveMarket(invert_in, "../data/torus_invert_point_input.mtx");
    Eigen::saveMarket(invert_out, "../data/torus_invert_point_output.mtx");

    return 0;
}





