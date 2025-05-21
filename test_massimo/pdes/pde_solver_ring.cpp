#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;

    constexpr int M = 2;

    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;
    
    std::string folder = "ring/";
    std::string save_path = "../ring/"; // or wherever you want
    std::string result_folder = "../results/" + folder + "/";

    std::ofstream file(result_folder + "L2_error.csv");
    std::ofstream file2(result_folder + "assembly_times.csv");



    std::vector<int> ref = {1,2,3,4,5,6};

    /*

    ScalarField<3, decltype([](const Vec& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);
    
        return 2 * x * y * (x*x + y*y - 1) * (4 - x*x - y*y);
    })> f_exact; // exact solution
    */


    ScalarField<M> f_exact(
        std::function<double(const Vec&)>(
            [&](const Vec& p) {
                double x = p(0);
                double y = p(1);
            
                return 2 * x * y * (x*x + y*y - 1) * (4 - x*x - y*y);
            }));



    
    ScalarField<M, decltype([](const Vec& p) {
        double x = p(0);
        double y = p(1);

        return  8 * x * y * (8 * x*x + 8 * y*y - 15 ) ; //poisson
        //return 80 * x * y * (8 * x*x + 8 * y*y - 15 ) +2*x*y*(6*std::pow(x,4) + 12*std::pow(x,2)*std::pow(y,2) - 20*std::pow(x,2) + 6*std::pow(y,4) - 20*std::pow(y,2) + 8);
    })> u; // forcing term

    
    VectorField<M, M, Fun> df_exact;
    
    df_exact(0, 0) = [](const Vec& p) {
        double x = p(0), y = p(1);
        return - 2 * y * (5 * std::pow(x, 4) + 6 * x * x * y * y - 15 * x * x + std::pow(y, 4) - 5 * y * y + 4);
    };
    
    df_exact(1, 0) = [](const Vec& p) {
        double x = p(0), y = p(1);
        return - 2 * x * (std::pow(x, 4) + 6 * x * x * y * y - 5 * x * x + 5 * std::pow(y, 4) - 15 * y * y + 4);
    };

    MatrixField<M, M, 1> b_;

    b_[0] = [](const Vec& p) {
        return -p(0);  // -x
    };
    
    b_[1] = [](const Vec& p) {
        return -p(1);  // -y
    };
    /*
    b_[2] = [](const Vec& p) {
        return 0.0;  // -y
    };
    */

    
    // Derivative in x-direction
    ScalarField<M, decltype([](const Vec& p) {
        const double x = p(0), y = p(1);
        const double x2 = x * x, x4 = x2 * x2;
        const double y2 = y * y, y4 = y2 * y2;

        return -2.0 * y * (5.0 * x4 + 6.0 * x2 * y2 - 15.0 * x2 + y4 - 5.0 * y2 + 4.0);
    })> df_exact_x;

    // Derivative in y-direction
    ScalarField<M, decltype([](const Vec& p) {
        const double x = p(0), y = p(1);
        const double x2 = x * x, x4 = x2 * x2;
        const double y2 = y * y, y4 = y2 * y2;

        return -2.0 * x * (x4 + 6.0 * x2 * y2 - 5.0 * x2 + 5.0 * y4 - 15.0 * y2 + 4.0);
    })> df_exact_y;
    


    for (const auto& r : ref){

        
        std::string path = "../../plots/data/" + folder + "/";
    
        //auto mesh = load_mesh(path);
        auto mesh = IsoMesh<2,2>::quarter_ring();
    
        
        
        mesh.refine_knots({r,r});

        // print the number of cells
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
        //mesh.refine_knots({7,7});
        // print the knots
        /*
        std::cout << "Knots: " << std::endl;
        for (int d = 0; d < 2; d++) {
            std::cout << "Dimension " << d << ": ";
            for (const auto& knot : mesh.knots()[d]) {
                std::cout << knot << " ";
            }
            std::cout << std::endl;
        }
        */
    
        
    
        //std::string save_path = "../ring" + std::to_string(r) + "/"; // or wherever you want
        //auto mesh2 = IsoMesh<2,3>::quarter_ring();
        //mesh2.refine_knots({r,r});
        //export_mesh(mesh2, save_path);
        
        IsoSpace Vh(mesh);
    
        TrialFunction f(Vh);
        TestFunction v(Vh);
        
    
        //auto start = std::chrono::high_resolution_clock::now();
        //auto a = integral(mesh,QGL2DP9)( 10 * dot(grad(f), grad(v)) +dot(b_, grad(f)) * v ); 
        auto a = integral(mesh,QGL2DP4)( dot(grad(f), grad(v))  );// dot(grad(f), grad(v)) laplacian(f)*laplacian(v) + dot(b,grad(f))*v
        //auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh,QGL2DP9)(u*v);
        //auto mm = integral(mesh,QGL2DP9)(f*v);
    
    
        auto& dof_handler = Vh.dof_handler();

        //std::cout << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        Eigen::SparseMatrix<double> A = a.assemble();
        auto end = std::chrono::high_resolution_clock::now();
        //std::cout << "Assembly time for A: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
        //start = std::chrono::high_resolution_clock::now();
        auto b = F.assemble();   
        //end = std::chrono::high_resolution_clock::now();
        //std::cout << "Assembly time for b: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;  
        //start = std::chrono::high_resolution_clock::now();
        dof_handler.set_hom_dirichlet_constraint();
        dof_handler.enforce_constraints(A,b);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Total time: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
        file2 << dof_handler.n_dofs() << "," << (end - start).count() << std::endl;

    
        
        
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        auto uh_full = solver.solve(b);
    
    
        //std::cout << "uh_full: " << uh_full.transpose() << std::endl;
    
        IsoFunction solution(Vh);
        solution =  uh_full;
    
        
    
        ScalarField<M> err_physical(
            std::function<double(const Vec&)>(
                [&](const Vec& p) {
                    double t1,t2;
                    auto u = mesh.invert_point(p,t1,t2,2); // assumes invert_point returns a Vector2d in parametric space
                    return solution(u) - f_exact(p);            // evaluate IsoFunction at u
                }));

        ScalarField<M> sol_physical(
            std::function<double(const Vec&)>(
                [&](const Vec& p) {
                    double t1,t2;
                    //std::cout<<"Ecco: "<<p;
                    auto u = mesh.invert_point(p,t1,t2,2); // assumes invert_point returns a Vector2d in parametric space
                    return solution(u);            // evaluate IsoFunction at u
                }));

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

        //df_appx(2, 0) = [&](const Vec& p) {
        //    return 0.0; // df_appx(1, 0) is not used
        //};

        

        //auto errorH1 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical + dot(df_appx - df_exact ,df_appx - df_exact )) ) ;
        auto errorL2 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical ) ) ;  // / std::sqrt(integral(mesh, QGL2DP9)(f_exact * f_exact)) //+ dot(grad(sol_physical) - df_exact, grad(sol_physical) - df_exact))


        double h_max = mesh.h_max();
    
    
        std::cout<<"L2 error: "<< errorL2 << std::endl;
        //std::cout<<"H1 error: "<< errorH1 << std::endl;
        //file <<h_max<<","<< errorL2 <<","<<errorH1<< std::endl;
        

        //export_results(mesh, solution, folder, std::make_optional(f_exact), /*nn = */10); //std::nullopt if there is no solution 
        

    }
        
        

    return 0;
}