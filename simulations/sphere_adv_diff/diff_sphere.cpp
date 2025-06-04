
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

    std::vector<int> ref = {20};

    auto u = adv_diff_sphere::make_rhs();
    auto b = adv_diff_sphere::make_b_field();

    for (const auto& r : ref) {
        auto mesh = IsoMesh<2, 3>::sphere_patch(1.0, 20.0, 60.0, -40.0, 50.0); 

        if (r > 0) mesh.refine_knots({r, r});
        
        double h_max = mesh.h_max();

        mesh.mark_boundary(0);
        mesh.mark_boundary(1, [](const auto& edge) {
            return edge.node(0)[1] == 0 && edge.node(1)[1] == 0;
        });


        // Set a periodic Spline basi
        //
        

        IsoSpace Vh(mesh);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        auto a = integral(mesh, QGL2DP9)(dot(grad(f), grad(v)) + dot(grad(f),b) );
        auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh, QGL2DP9)(u * v);

        auto& dof_handler = Vh.dof_handler();
        Eigen::SparseMatrix<double> A = a.assemble();
        auto b = F.assemble();
        auto c = m.assemble();

                
        dof_handler.set_hom_dirichlet_constraint(1);
        dof_handler.enforce_constraints(A,b);
    

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd uh_full = solver.solve(b);

        // Create IsoFunction
        IsoFunction solution(Vh);
        solution = uh_full;

        //std::cout<<"uhfull"<<uh_full.transpose()<<"\n";

        // Export mesh and solution
        std::string level_path = save_path + "ref" + std::to_string(r) + "/mesh/";
        helpers::export_mesh(mesh, level_path);
        std::string solution_path = save_path + "ref" + std::to_string(r) + "/solution/";
        helpers::export_results(mesh, solution, solution_path, 10);

        std::cout << "\n===========================================\n";
        std::cout << "Refinement level: " << r << "\n";
        std::cout << "Number of cells : " << mesh.n_cells() << "\n";
        std::cout << "h_max           : " << h_max << "\n";
        std::cout << "Mesh exported to: " << level_path << "\n";
        std::cout << "PDE results to  : " << solution_path << "\n";
        std::cout << "===========================================\n";
    }

    return 0;
}