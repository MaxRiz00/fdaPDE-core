
#include "isogeometric.h"
#include "exact_solution.h"
#include "../helpers.h"

using namespace fdapde;

int main() {

    constexpr int M = 3;
    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;

    std::string folder = "results/";


    std::string save_path = "../" + folder  ;
    std::filesystem::create_directories(save_path);
    std::ofstream file(save_path + "L2_error.csv");
    file << "h_max,L2_error,H1_error\n";

    std::vector<int> ref = {0, 1, 2, 3, 4, 5};

    auto f_exact = diff_torus::make_u_exact();
    auto u = diff_torus::make_rhs();
    auto df_exact = diff_torus::make_grad_u_exact();

    for (const auto& r : ref) {
        std::cout << "\n=== Refinement level: " << r << " ===\n";

        auto mesh = IsoMesh<2, 3>::torus();
        if (r > 0) mesh.refine_knots({r, r});
        double h_max = mesh.h_max();

        std::cout << "Number of cells: " << mesh.n_cells() << "\n";
        std::cout << "h_max: " << h_max << "\n";

        // Set a periodic Spline basis

        std::array<std::vector<double>,2> open_uniform_knots;
        std::array<int,2> basis_dims;
        std::array<int,2> new_degree;
        for(int i = 0; i < 2; i++){
            new_degree[i] = mesh.degree()[i] ;
        }
        for(int i = 0; i < 2; i++){
            open_uniform_knots[i] = pad_knots(mesh.param_nodes()[i], new_degree[i]);
            basis_dims[i] = open_uniform_knots[i].size() - new_degree[i] - 1;
        }
        MdArray<double, full_dynamic_extent_t<2>> unitary_weights;
        unitary_weights.resize(basis_dims);
        unitary_weights.set_constant(1.0);
        auto basis_pde = NurbsBasis<2>(open_uniform_knots, unitary_weights, new_degree, mesh.is_periodic()); //bas
        //
        

        IsoSpace Vh(mesh, basis_pde);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        auto a = integral(mesh, QGL2DP9)(dot(grad(f), grad(v)));
        auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh, QGL2DP9)(u * v);

        auto& dof_handler = Vh.dof_handler();
        Eigen::SparseMatrix<double> A = a.assemble();
        auto b = F.assemble();
        auto c = m.assemble();

                
        dof_handler.enforce_periodic_constraints(A,b);
        dof_handler.enforce_periodic_constraints(c);

        // Solve system with constraint (e.g., for unique solution on closed surface)
        int counter = b.size();
        Eigen::SparseMatrix<double> Zero(1, 1);
        SparseBlockMatrix<double, 2, 2> D(A, c.sparseView(), c.transpose().sparseView(), Zero);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(D);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(counter + 1);
        rhs.head(counter) = b;
        Eigen::VectorXd uh_reduced = solver.solve(rhs).head(counter);

        Eigen::VectorXd uh_full = dof_handler.expand_solution(uh_reduced);

        // Create IsoFunction
        IsoFunction solution(Vh);
        solution = uh_full;

        ScalarField<M> err_physical(
            [&](const Vec& p) {
                double t1, t2;
                auto u = mesh.invert_point(p, t1, t2, 2);
                return solution(u) - f_exact(p);
            });

        VectorField<M, M, Fun> df_appx;
        df_appx(0, 0) = [&](const Vec& p) {
            double t1, t2;
            auto u = mesh.invert_point(p, t1, t2, 2);
            return solution.phys_grad(u)(0);
        };
        df_appx(1, 0) = [&](const Vec& p) {
            double t1, t2;
            auto u = mesh.invert_point(p, t1, t2, 2);
            return solution.phys_grad(u)(1);
        };
        df_appx(2, 0) = [&](const Vec& p) {
            double t1, t2;
            auto u = mesh.invert_point(p, t1, t2, 2);
            return solution.phys_grad(u)(2);
        };

        auto errorL2 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical));
        auto errorH1 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical + dot(df_appx - df_exact, df_appx - df_exact)));

        std::cout << "L2 error: " << errorL2 << "\n";
        std::cout << "H1 error: " << errorH1 << "\n";
        file << h_max << "," << errorL2 << "," << errorH1 << "\n";

        // Export mesh and solution
        std::string level_path = save_path + "ref" + std::to_string(r) + "/mesh/";
        helpers::export_mesh(mesh, level_path);
        std::string solution_path = save_path + "ref" + std::to_string(r) + "/solution/";
        helpers::export_results(mesh, solution, solution_path, 10);
    }

    return 0;
}