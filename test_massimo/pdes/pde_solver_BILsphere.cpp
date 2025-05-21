#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using matrix_t = Eigen::SparseMatrix<double>;

int main() {
    std::string folder = "sphere/";
    std::string save_path = "../sphere/"; // or wherever you want
    std::string result_folder = "../results/" + folder + "/";

    std::ofstream file(result_folder + "L2_error.csv");
    file << "h_max,L2_error\n";

    constexpr int M = 3;

    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;

    std::vector<int> ref_levels = {0};

    for (int r : ref_levels) {
        std::cout << "Refinement level: " << r << std::endl;

        auto mesh = IsoMesh<2,3>::sphere();
        if(r > 0) mesh.refine_knots({r, r});
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;

        IsoSpace Vh(mesh);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        constexpr int alpha = 3;
        constexpr int beta = 5;
        constexpr int lambda = 3;

        // Forcing term (Laplacian of u_exact)
        ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
            double x = p(0), y = p(1), z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r);
            return lambda * lambda *(lambda + 1) * (lambda + 1) * std::sin(lambda * phi) * std::pow(std::sin(theta),lambda) ;
        })> u;

        // Exact solution
        ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
            double x = p(0), y = p(1), z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r);
            return std::sin(lambda * phi) * std::pow(std::sin(theta), lambda);
        })> u_exact;



        VectorField<M, M, Fun> df_exact;

        df_exact(0, 0) = [=](const Vec& p) {
            double x = p(0), y = p(1), z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r);
        
            double df_dtheta = lambda *std::cos(theta) * std::sin(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
            double df_dphi = lambda * std::cos(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
        
            return (x * z * df_dtheta - y * df_dphi)/std::sqrt(x*x + y*y);
        };
        
        df_exact(1, 0) = [=](const Vec& p) {
            double x = p(0), y = p(1), z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r);
        
            double df_dtheta = lambda *std::cos(theta) * std::sin(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
            double df_dphi = lambda * std::cos(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
        
            return (y * z * df_dtheta + x * df_dphi)/std::sqrt(x*x + y*y);
        };
        
        df_exact(2, 0) = [=](const Vec& p) {
            double x = p(0), y = p(1), z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r);
        
            double df_dtheta = lambda *std::cos(theta) * std::sin(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
            double df_dphi = lambda * std::cos(lambda * phi) * std::pow(std::sin(theta), lambda - 1);
            return -df_dtheta * std::sqrt(x*x + y*y) ;
        };

        // Assemble system
        auto a = integral(mesh, QGL2DP9)(laplacian(f) *laplacian(v));
        auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh, QGL2DP9)(u * v);

        auto& dof_handler = Vh.dof_handler();
        Eigen::SparseMatrix<double> A = a.assemble();
        Eigen::VectorXd b = F.assemble();
        Eigen::VectorXd c = m.assemble();

        // Apply periodic BC reduction
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
                int i = dof_map[it.row()];
                int j = dof_map[it.col()];
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

        std::cout<<"b reduced: " << A_reduced.toDense() << std::endl;

        // Solve system with constraint (e.g., for unique solution on closed surface)
        Eigen::SparseMatrix<double> Zero(1, 1);
        SparseBlockMatrix<double, 2, 2> D(A_reduced, c_reduced.sparseView(), c_reduced.transpose().sparseView(), Zero);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(D);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(counter + 1);
        rhs.head(counter) = b_reduced;
        Eigen::VectorXd uh_reduced = solver.solve(rhs).head(counter);

        // Expand to full solution
        Eigen::VectorXd uh_full(dof_map.size());
        for (int i = 0; i < dof_map.size(); ++i) {
            int mapped = dof_map[i];
            uh_full[i] = uh_reduced[reduced_indices[mapped]];
        }

        // Create IsoFunction
        IsoFunction solution(Vh);
        solution = uh_full.topRows(A.rows());

        // Define error field
        ScalarField<3> err_field(
            std::function<double(const Eigen::Matrix<double, 3, 1>& p)>(
                [&](const Eigen::Matrix<double, 3, 1>& p) {
                    double t1, t2;
                    auto u = mesh.invert_point(p, t1, t2, 5);
                    return solution(u) - u_exact(p);
                }
            )
        );


        VectorField<M, M, Fun> df_appx;

        df_appx(0, 0) = [&](const Vec& p) {
            double t1,t2;
            auto u = mesh.invert_point(p,t1,t2,5);
            return solution.phys_grad(u)(0);
        };
        df_appx(1, 0) = [&](const Vec& p) {
            double t1,t2;
            auto u = mesh.invert_point(p,t1,t2,5);
            return solution.phys_grad(u)(1);
        };

        df_appx(2, 0) = [&](const Vec& p) {
            double t1,t2;
            auto u = mesh.invert_point(p,t1,t2,5);
            return solution.phys_grad(u)(2);
        };

        double error_L2 = std::sqrt(integral(mesh, QGL2DP9)(err_field * err_field));
        auto error_H1 = std::sqrt(integral(mesh, QGL2DP9)(err_field * err_field + dot(df_appx - df_exact ,df_appx - df_exact )) ) ;
        double h_max = mesh.h_max();

        std::cout << "h_max: " << h_max << ", L2 error: " << error_L2 << ", H1 error: "<<error_H1<<std::endl;
        file << h_max << "," << error_L2 << ","<<error_H1<< std::endl;

        // Optional: export for visualization
        std::string result_folder = "../results/sphere" + std::to_string(r) + "/";
        export_results(mesh, solution, result_folder, std::make_optional(u_exact), 20);
    }

    return 0;
}