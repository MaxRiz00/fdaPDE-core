
#include "isogeometric.h"
#include "exact_solution.h"
#include "../helpers.h"

using namespace fdapde;

int main() {

    constexpr int M = 2;
    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;

    std::string folder = "results/";


    std::string save_path = "../" + folder  ;
    std::filesystem::create_directories(save_path);
    std::ofstream file(save_path + "L2_error.csv");
    file << "h_max,L2_error,H1_error\n";

    std::vector<int> ref = {1, 2, 3, 4, 5, 6};

    auto f_exact = diff_ring::make_u_exact();
    auto u = diff_ring::make_rhs();
    auto df_exact = diff_ring::make_grad_u_exact();

    for (const auto& r : ref) {
        std::cout << "\n=== Refinement level: " << r << " ===\n";

        auto mesh = IsoMesh<2, 2>::quarter_ring();
        mesh.refine_knots({r, r});
        double h_max = mesh.h_max();

        std::cout << "Number of cells: " << mesh.n_cells() << "\n";
        std::cout << "h_max: " << h_max << "\n";

        std::string level_path = save_path + "ref" + std::to_string(r) + "/mesh/";
        helpers::export_mesh(mesh, level_path);

        IsoSpace Vh(mesh);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        auto a = integral(mesh, QGL2DP9)(dot(grad(f), grad(v)));
        auto F = integral(mesh, QGL2DP9)(u * v);

        auto& dof_handler = Vh.dof_handler();
        Eigen::SparseMatrix<double> A = a.assemble();
        auto b = F.assemble();

        dof_handler.set_hom_dirichlet_constraint();
        dof_handler.enforce_constraints(A, b);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        auto uh_full = solver.solve(b);

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

        auto errorL2 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical));
        auto errorH1 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical + dot(df_appx - df_exact, df_appx - df_exact)));

        std::cout << "L2 error: " << errorL2 << "\n";
        std::cout << "H1 error: " << errorH1 << "\n";
        file << h_max << "," << errorL2 << "," << errorH1 << "\n";

        std::string solution_path = save_path + "ref" + std::to_string(r) + "/solution/";
        helpers::export_results(mesh, solution, solution_path, 10);
    }

    return 0;
}