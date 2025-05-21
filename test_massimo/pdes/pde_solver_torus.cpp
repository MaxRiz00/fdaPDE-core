#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using matrix_t = Eigen::SparseMatrix<double>;

int main() {
    std::string folder = "torus";
    std::string result_folder = "../results/" + folder + "/";
    std::ofstream file(result_folder + "L2_error.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << result_folder + "L2_error.csv" << std::endl;
        return 1;
    }
    file << "h_max,L2_error,H1_error\n";

    constexpr int M = 3;
    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;

    constexpr int m = 2, n = 2;

    std::vector<int> ref_levels = {0, 1, 2, 3, 4, 5};  // Customize levels

    for (int r : ref_levels) {
        std::cout << "Refinement level: " << r << std::endl;

        auto mesh = IsoMesh<2,3>::torus();
        if (r > 0) mesh.refine_knots({r, r});

        IsoSpace Vh(mesh);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        // Exact solution
        ScalarField<3, decltype([](const Vec& p) {
            double R = 2, r = 1;
            double x = p(0), y = p(1), z = p(2);
            double phi = std::atan2(y,x);
            double theta = std::atan2(z, std::sqrt(x*x + y*y) - R);
            return std::sin(m * theta) * std::sin(n * phi);
        })> u_exact;

        // Forcing term (Laplacian)
        ScalarField<3, decltype([](const Vec& p) {
            double R = 2, r = 1;
            double x = p(0), y = p(1), z = p(2);
            double phi = std::atan2(y,x);
            double theta = std::atan2(z, std::sqrt(x*x + y*y) - R);

            return (m * std::sin(n * phi) * std::cos(m * theta) * std::sin(theta) / (r * (R + r * std::cos(theta)))) +
                   std::sin(n * phi) * std::sin(m * theta) * (m*m / (r*r) + n*n / std::pow(R + r * std::cos(theta), 2));
        })> u;

        // Assemble system
        auto a = integral(mesh, QGL2DP9)(dot(grad(f), grad(v)));
        auto mass = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh, QGL2DP9)(u * v);

        auto& dof_handler = Vh.dof_handler();
        Eigen::SparseMatrix<double> A = a.assemble();
        Eigen::VectorXd b = F.assemble();
        Eigen::VectorXd c = mass.assemble();

        const auto& dof_map = dof_handler.dof_map();
        std::unordered_map<int, int> reduced_indices;
        std::vector<int> keep_dofs;
        int counter = 0;

        for (int i = 0; i < dof_map.size(); ++i) {
            if (dof_map[i] == i) {
                reduced_indices[i] = counter++;
                keep_dofs.push_back(i);
            }
        }

        Eigen::SparseMatrix<double> A_reduced(counter, counter);
        Eigen::VectorXd b_reduced = Eigen::VectorXd::Zero(counter);
        Eigen::VectorXd c_reduced = Eigen::VectorXd::Zero(counter);

        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
                int i = dof_map[it.row()], j = dof_map[it.col()];
                if (reduced_indices.count(i) && reduced_indices.count(j)) {
                    A_reduced.coeffRef(reduced_indices[i], reduced_indices[j]) += it.value();
                }
            }
        }

        for (int i = 0; i < b.size(); ++i) {
            int mapped = dof_map[i];
            if (reduced_indices.count(mapped)) {
                b_reduced[reduced_indices[mapped]] += b[i];
                c_reduced[reduced_indices[mapped]] += c[i];
            }
        }

        Eigen::SparseMatrix<double> Zero(1, 1);
        SparseBlockMatrix<double, 2, 2> D(A_reduced, c_reduced.sparseView(), c_reduced.transpose().sparseView(), Zero);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(D);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(counter + 1);
        rhs.head(counter) = b_reduced;
        Eigen::VectorXd uh_reduced = solver.solve(rhs).head(counter);

        // Expand full solution
        Eigen::VectorXd uh_full(dof_map.size());
        for (int i = 0; i < dof_map.size(); ++i) {
            int mapped = dof_map[i];
            uh_full[i] = uh_reduced[reduced_indices[mapped]];
        }

        IsoFunction solution(Vh);
        solution = uh_full.topRows(A.rows());

        // Approximate gradient field
        VectorField<M, M, Fun> df_appx;
        for (int i = 0; i < M; ++i) {
            df_appx(i, 0) = [&, i](const Vec& p) {
                double t1, t2;
                auto u = mesh.invert_point(p, t1, t2, 5);
                return solution.phys_grad(u)(i);
            };
        }

        // Exact gradient field (replace with your own known exact gradient)
        VectorField<M, M, Fun> df_exact;
        df_exact(0, 0) = [](const Vec& p) {
            double R = 2, r = 1;
            double x = p(0), y = p(1), z = p(2);

            double phi = std::atan2(y, x);
            double rho = std::sqrt(x * x + y * y);
            double theta = std::atan2(z, rho - R);

            double sin_phi = std::sin(phi), cos_phi = std::cos(phi);
            double sin_theta = std::sin(theta), cos_theta = std::cos(theta);

            double dphi = n * std::cos(n * phi) * std::sin(m * theta);
            double dtheta = m * std::cos(m * theta) * std::sin(n * phi);

            double denom_phi = (R + r * cos_theta) * (R + r * cos_theta);
            double denom_theta = r * r;

            // ∂Φ/∂phi (X component)
            double dphi_x = -(R + r * cos_theta) * sin_phi;

            // ∂Φ/∂theta (X component)
            double dtheta_x = -r * sin_theta * cos_phi;

            return (dphi / denom_phi) * dphi_x + (dtheta / denom_theta) * dtheta_x;
        };

        df_exact(1, 0) = [](const Vec& p) {
            double R = 2, r = 1;
            double x = p(0), y = p(1), z = p(2);

            double phi = std::atan2(y, x);
            double rho = std::sqrt(x * x + y * y);
            double theta = std::atan2(z, rho - R);

            double sin_phi = std::sin(phi), cos_phi = std::cos(phi);
            double sin_theta = std::sin(theta), cos_theta = std::cos(theta);

            double dphi = n * std::cos(n * phi) * std::sin(m * theta);
            double dtheta = m * std::cos(m * theta) * std::sin(n * phi);

            double denom_phi = (R + r * cos_theta) * (R + r * cos_theta);
            double denom_theta = r * r;

            // ∂Φ/∂phi (Y component)
            double dphi_y = (R + r * cos_theta) * cos_phi;

            // ∂Φ/∂theta (Y component)
            double dtheta_y = -r * sin_theta * sin_phi;

            return (dphi / denom_phi) * dphi_y + (dtheta / denom_theta) * dtheta_y;
        };

        df_exact(2, 0) = [](const Vec& p) {
            double R = 2, r = 1;
            double x = p(0), y = p(1), z = p(2);

            double phi = std::atan2(y, x);
            double rho = std::sqrt(x * x + y * y);
            double theta = std::atan2(z, rho - R);

            double sin_phi = std::sin(phi);
            double cos_theta = std::cos(theta);

            double dtheta = m * std::cos(m * theta) * std::sin(n * phi);
            double denom_theta = r * r;

            // ∂Φ/∂theta (Z component)
            double dtheta_z = r * cos_theta;

            return (dtheta / denom_theta) * dtheta_z;
        };

        // Error scalar field
        ScalarField<3> err_field(
            std::function<double(const Vec& p)>(
                [&](const Vec& p) {
                    double t1, t2;
                    auto u = mesh.invert_point(p, t1, t2, 5);
                    return solution(u) - u_exact(p);
                }
            )
        );

        // H1 error term
        VectorField<M, M, Fun> grad_diff;
        for (int i = 0; i < M; ++i) {
            grad_diff(i, 0) = [=](const Vec& p) {
                return df_appx(i, 0)(p) - df_exact(i, 0)(p);
            };
        }

        double L2_error = std::sqrt(integral(mesh, QGL2DP9)(err_field * err_field));
        double H1_error = std::sqrt(L2_error * L2_error + integral(mesh, QGL2DP9)(dot(grad_diff, grad_diff)));

        double h_max = mesh.h_max();

        std::cout << "h_max: " << h_max << ", L2 error: " << L2_error << ", H1 error: " << H1_error << std::endl;
        file << h_max << "," << L2_error << "," << H1_error << "\n";

        std::string export_path = "../results/" + folder + std::to_string(r) + "/";
        // create directory if it doesn't exist
        std::cout << "Exporting mesh to: " << export_path << std::endl;
        export_results(mesh, solution, export_path, std::make_optional(u_exact), 10);
    }

    return 0;
}