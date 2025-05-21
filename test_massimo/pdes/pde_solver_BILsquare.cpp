#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;

    using Vec2 = Eigen::Matrix<double, 2, 1>;
    using Fun = std::function<double(const Vec2&)>;
    
    std::string folder = "ring/";
    std::string save_path = "../ring/"; // or wherever you want
    std::string result_folder = "../results/" + folder + "/";

    std::ofstream file(result_folder + "L2_error.csv");
    std::ofstream file2(result_folder + "assembly_times.csv");



    std::vector<int> ref = {1,2,3,4,5,6};


    ScalarField<2> u(
        std::function<double(const Vec2&)>(
            [](const Vec2& p) {
                double x = p(0);
                double y = p(1);
    
                return -8 * std::pow(M_PI, 4) * (-4 * std::pow(std::sin(M_PI * x), 2) * std::pow(std::sin(M_PI * y), 2)
                                                 + std::pow(std::sin(M_PI * x), 2)
                                                 + std::pow(std::sin(M_PI * y), 2))
                       + 4 * std::pow(M_PI, 4) * (std::cos(2 * M_PI * (x - y)) + std::cos(2 * M_PI * (x + y)));
            }));

    ScalarField<2, decltype([](const Vec2& p) {
        double x = p(0);
        double y = p(1);
    
        return std::pow(std::sin(M_PI * x), 2) * std::pow(std::sin(M_PI * y), 2);
    })> f_exact;


    
    fdapde::VectorField<2, 2, Fun> df_exact;

    df_exact(0, 0) = [](const Vec2& p) {
        double x = p(0), y = p(1);
        return 2 * M_PI * std::cos(M_PI * x) * std::sin(M_PI * x) * std::pow(std::sin(M_PI * y), 2);
    };
    
    df_exact(1, 0) = [](const Vec2& p) {
        double x = p(0), y = p(1);
        return 2 * M_PI * std::cos(M_PI * y) * std::sin(M_PI * y) * std::pow(std::sin(M_PI * x), 2);
    };

    for (const auto& r : ref){

        
        std::string path = "../../plots/data/" + folder + "/";
    
        //auto mesh = load_mesh(path);
        auto mesh = IsoMesh<2,2>::square();
    

        
        mesh.refine_knots({r,r});
        // print the number of cells
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;

        // print the knot vector
            // print the knots
            /*
        for(int i = 0; i < 2; i++){
            std::cout<<"Knots "<<i<<" : ";
            for(auto k : mesh.knots()[i]){
                std::cout<<k<<" ";
            }
            std::cout<<std::endl;
        }
            */
        
    
        std::string save_path = "../ring" + std::to_string(r) + "/"; // or wherever you want
        auto mesh2 = IsoMesh<2,3>::quarter_ring();
        mesh2.refine_knots({r,r});
        export_mesh(mesh2, save_path);
        
        IsoSpace Vh(mesh);
    
        TrialFunction f(Vh);
        TestFunction v(Vh);
        
    
        //auto start = std::chrono::high_resolution_clock::now();
        auto a = integral(mesh,QGL2DP9)(laplacian(f)*laplacian(v)); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
        //auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh,QGL2DP9)(u*v);
        //auto mm = integral(mesh,QGL2DP9)(f*v);
    
    
        auto& dof_handler = Vh.dof_handler();

        auto start = std::chrono::high_resolution_clock::now();

        Eigen::SparseMatrix<double> A = a.assemble();
        auto b = F.assemble();  
          
        dof_handler.set_clamped_hom_constraint();
        dof_handler.enforce_constraints(A,b);
        auto end = std::chrono::high_resolution_clock::now();
        
        
        
        std::chrono::duration<double> duration = end - start;
        
        std::cout << "Assembly time: " << duration.count() << " seconds" << std::endl;
        file2 << dof_handler.n_dofs() << "," << duration.count() << std::endl;

    
        
        
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        auto uh_full = solver.solve(b);
    
    
        //std::cout << "uh_full: " << uh_full.transpose() << std::endl;
    
        IsoFunction solution(Vh);
        solution =  uh_full;
    
        
    
        ScalarField<2> err_physical(
            std::function<double(const Vec2&)>(
                [&](const Vec2& p) {
                    double t1,t2;
                    auto u = mesh.invert_point(p,t1,t2,5); // assumes invert_point returns a Vector2d in parametric space
                    return solution(u) - f_exact(p);            // evaluate IsoFunction at u
                }));

        ScalarField<2> sol_physical(
            std::function<double(const Vec2&)>(
                [&](const Vec2& p) {
                    double t1,t2;
                    //std::cout<<"Ecco: "<<p;
                    auto u = mesh.invert_point(p,t1,t2,5); // assumes invert_point returns a Vector2d in parametric space
                    return solution(u);            // evaluate IsoFunction at u
                }));

        fdapde::VectorField<2, 2, Fun> df_appx;

        df_appx(0, 0) = [&](const Vec2& p) {
            double t1,t2;
            auto u = mesh.invert_point(p,t1,t2,2);
            return solution.phys_grad(u)(0);
        };
        df_appx(1, 0) = [&](const Vec2& p) {
            double t1,t2;
            auto u = mesh.invert_point(p,t1,t2,2);
            return solution.phys_grad(u)(1);
        };

        

        auto errorH1 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical + dot(df_appx - df_exact ,df_appx - df_exact )) ) ;
        auto errorL2 = std::sqrt(integral(mesh, QGL2DP9)(err_physical * err_physical ) ) ;  // / std::sqrt(integral(mesh, QGL2DP9)(f_exact * f_exact)) //+ dot(grad(sol_physical) - df_exact, grad(sol_physical) - df_exact))


        double h_max = mesh.h_max();

        std::cout<<"h_max: "<< h_max << std::endl;
    
    
        std::cout<<"L2 error: "<< errorL2 << std::endl;
        std::cout<<"H1 error: "<< errorH1 << std::endl;
        file <<h_max<<","<< errorL2 <<","<<errorH1<< std::endl;

        //export_results(mesh, solution, folder, std::make_optional(f_exact), /*nn = */ 30); //std::nullopt if there is no solution 

    }
        
        

    return 0;
}