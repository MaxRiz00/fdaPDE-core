
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

    std::vector<int> ref = {0};

    auto f_exact = advdiff_ring::make_u_exact();
    auto u = advdiff_ring::make_rhs();
    auto df_exact = advdiff_ring::make_grad_u_exact();
    auto g_N = advdiff_ring::make_g_neumann();

    for (const auto& r : ref) {
        std::cout << "\n=== Refinement level: " << r << " ===\n";

        Eigen::Matrix<double, 2, 2> K;
        K << 1, 0,
             0, 1;



        auto mesh = IsoMesh<2, 2>::square();
        if(r >0 ) mesh.refine_knots({r, r});
        double h_max = mesh.h_max();

        std::cout << "Number of cells: " << mesh.n_cells() << "\n";
        std::cout << "h_max: " << h_max << "\n";

        

        IsoSpace Vh(mesh);
        TrialFunction f(Vh);
        TestFunction v(Vh);

        auto a = integral(mesh,QGL2DP9)(dot(K * grad(f),grad(v)));
        Eigen::SparseMatrix<double> A = a.assemble();
        
        std::cout << "A: " << Eigen::MatrixXd(A) << "\n";
        
        auto F = integral(mesh,QGL2DP9)(u*v) ; //integral(mesh.boundary(1),QGL1DP3)(g_N * v)  
        //std::cout<<"Ecco F: " << F << "\n";
        auto c = integral(mesh,QGL2DP9)(v).assemble(); // 
        auto& dof_handler = Vh.dof_handler();
        //auto b = F.assemble();
        //std::cout<< "b: "<< b.transpose() << "\n";
        //make A full
        //std::cout<<"A: "<<Eigen::MatrixXd(A)<< "\n";
        
        /*
        
        int counter = b.size();
        Eigen::SparseMatrix<double> Zero(1, 1);
        SparseBlockMatrix<double, 2, 2> D(A, c.sparseView(), c.transpose().sparseView(), Zero);

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(D);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(counter + 1);
        rhs.head(counter) = b;
        Eigen::VectorXd uh_full = solver.solve(rhs).head(counter);
        

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

        // Export mesh and solution
        std::string level_path = save_path + "ref" + std::to_string(r) + "/mesh/";
        helpers::export_mesh(mesh, level_path);
        std::string solution_path = save_path + "ref" + std::to_string(r) + "/solution/";
        helpers::export_results(mesh, solution, solution_path, 10);
        */
        
   
    }

    return 0;
}