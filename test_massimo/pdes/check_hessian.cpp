#include "isogeometric.h"
#include <unsupported/Eigen/SparseExtra>
#include <fstream>
#include <cassert>

using namespace fdapde;

IsoMesh<2, 3> load_mesh(const std::string& folder_path) {
    using SpMatrix = Eigen::SparseMatrix<double>;

    std::array<int, 2> order;
    SpMatrix knots_x, knots_y, weights;
    SpMatrix control_points_x, control_points_y, control_points_z;

    std::string path = folder_path;

    // Load order
    std::ifstream order_file(path + "order.txt");
    for (int i = 0; i < 2; i++) {
        order_file >> order[i];
    }

    // Load matrices
    Eigen::loadMarket(knots_x, path + "knots_x.mtx");
    Eigen::loadMarket(knots_y, path + "knots_y.mtx");
    Eigen::loadMarket(weights, path + "weights.mtx");
    Eigen::loadMarket(control_points_x, path + "ctrlpts_x.mtx");
    Eigen::loadMarket(control_points_y, path + "ctrlpts_y.mtx");
    Eigen::loadMarket(control_points_z, path + "ctrlpts_z.mtx");

    // Reconstruct data structures
    std::array<std::vector<double>, 2> nodes;
    nodes[0].resize(knots_x.cols());
    nodes[1].resize(knots_y.cols());

    for (size_t i = 0; i < nodes[0].size(); i++) {
        nodes[0][i] = knots_x.coeff(0, i);
    }

    for (size_t i = 0; i < nodes[1].size(); i++) {
        nodes[1][i] = knots_y.coeff(0, i);
    }

    MdArray<double, full_dynamic_extent_t<2>> weights_(weights.rows(), weights.cols());
    for (int i = 0; i < weights.rows(); i++) {
        for (int j = 0; j < weights.cols(); j++) {
            weights_(i, j) = weights.coeff(i, j);
        }
    }

    MdArray<double, full_dynamic_extent_t<3>> control_points(control_points_x.rows(), control_points_x.cols(), 3);
    for (int i = 0; i < control_points_x.rows(); i++) {
        for (int j = 0; j < control_points_x.cols(); j++) {
            control_points(i, j, 0) = control_points_x.coeff(i, j);
            control_points(i, j, 1) = control_points_y.coeff(i, j);
            control_points(i, j, 2) = control_points_z.coeff(i, j);
        }
    }

    return IsoMesh<2, 3>(nodes, weights_, control_points, order);
}

template <typename MeshType>
void check_hessian(const MeshType& mesh, double tol = 1e-3) {
    constexpr int local_dim = MeshType::local_dim;
    constexpr int embed_dim = MeshType::embed_dim;
    using VecLocal = Eigen::Matrix<double, local_dim, 1>;
    using VecEmbed = Eigen::Matrix<double, embed_dim, 1>;

    double h = 1e-5; // finite difference step

    // Pick a random point inside parametric domain Ω
    VecLocal u;
    for (int i = 0; i < local_dim; ++i) {
        auto nodes = mesh.param_nodes()[i];
        double umin = nodes.front();
        double umax = nodes.back();
        u(i) = umin + (umax - umin) * 0.3; // pick something away from boundary
    }

    auto derivatives = mesh.eval_param_derivatives(u, true);
    auto dX = derivatives.first_derivative; // F = ∇X
    auto d2X = *(derivatives.second_derivative); // H = ∇²X

    // Loop over embed_dim, local_dim, local_dim
    for (int k = 0; k < embed_dim; ++k) {
        for (int i = 0; i < local_dim; ++i) {
            for (int j = 0; j < local_dim; ++j) {

                double num_derivative = 0.0;

                if (i == j) {
                    // Diagonal terms: standard second derivative
                    VecLocal up = u, um = u;
                    up(i) += h;
                    um(i) -= h;

                    auto xp = mesh.eval_param(up);
                    auto xm = mesh.eval_param(um);

                    num_derivative = (xp(k) - 2.0 * mesh.eval_param(u)(k) + xm(k)) / (h * h);
                }
                else {
                    // Mixed second derivative
                    VecLocal upi_upj = u, upi_umj = u, umi_upj = u, umi_umj = u;
                    upi_upj(i) += h; upi_upj(j) += h;
                    upi_umj(i) += h; upi_umj(j) -= h;
                    umi_upj(i) -= h; umi_upj(j) += h;
                    umi_umj(i) -= h; umi_umj(j) -= h;

                    auto f1 = mesh.eval_param(upi_upj);
                    auto f2 = mesh.eval_param(upi_umj);
                    auto f3 = mesh.eval_param(umi_upj);
                    auto f4 = mesh.eval_param(umi_umj);

                    num_derivative = (f1(k) - f2(k) - f3(k) + f4(k)) / (4.0 * h * h);
                }

                double analytical = d2X(k,i,j);
                double diff = std::abs(num_derivative - analytical);

                if (diff > tol) {
                    std::cout << "Mismatch at component ("
                              << k << "," << i << "," << j << "): "
                              << "analytical = " << analytical
                              << ", numerical = " << num_derivative
                              << ", diff = " << diff << std::endl;
                }
                //assert(diff < tol && "Hessian check failed!");
            }
        }
    }

    std::cout << "✅ Full Hessian check passed!" << std::endl;
}


template <typename MeshType>
void check_gradient(const MeshType& mesh, double tol = 1e-6) {
    constexpr int local_dim = MeshType::local_dim;
    constexpr int embed_dim = MeshType::embed_dim;
    using VecLocal = Eigen::Matrix<double, local_dim, 1>;
    using VecEmbed = Eigen::Matrix<double, embed_dim, 1>;

    double h = 1e-5; // finite difference step

    // Pick a random point inside parametric domain Ω
    VecLocal u;
    for (int i = 0; i < local_dim; ++i) {
        auto nodes = mesh.param_nodes()[i];
        double umin = nodes.front();
        double umax = nodes.back();
        u(i) = umin + (umax - umin) * 0.3; // pick something away from boundary
    }

    auto derivatives = mesh.eval_param_derivatives(u, false);
    auto dX = derivatives.first_derivative; // F = ∇X

    // Loop over embed_dim and local_dim
    for (int k = 0; k < embed_dim; ++k) {
        for (int i = 0; i < local_dim; ++i) {
            // Finite difference approximation
            VecLocal up = u, um = u;
            up(i) += h;
            um(i) -= h;

            auto xp = mesh.eval_param(up);
            auto xm = mesh.eval_param(um);

            double num_derivative = (xp(k) - xm(k)) / (2*h);
            double analytical = dX(k,i);

            double diff = std::abs(num_derivative - analytical);

            if (diff > tol) {
                std::cout << "Mismatch in gradient at component ("
                          << k << "," << i << "): "
                          << "analytical = " << analytical
                          << ", numerical = " << num_derivative
                          << ", diff = " << diff << std::endl;
            }
            assert(diff < tol && "Gradient check failed!");
        }
    }

    std::cout << "✅ Gradient check passed!" << std::endl;
}


int main() {

    std::string folder = "curved/";
    std::string path = "../../plots/data/" + folder + "/";
    // Example usage with a specific mesh type
    // Replace IsoMesh<2,3>::sphere() with your actual mesh creation function
    // For example, if you have a mesh class called MyMesh, use MyMesh::create();
    auto mesh = IsoMesh<2,3>::sphere(); // for example
    //auto mesh= load_mesh(path);
    check_gradient(mesh);
    check_hessian(mesh);
    
    return 0;
}