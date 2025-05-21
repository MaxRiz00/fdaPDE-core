#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;

    constexpr int M = 3;

    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;
    
    std::string folder = "curly_plate_ref0/";
    std::string save_path = "../curly_plate_ref0/"; // or wherever you want
    std::string result_folder = "../results/" + folder + "/";

    std::ofstream file2(result_folder + "assembly_times.csv");



    std::vector<int> ref = {1,2,3,4,5,6};

    
    ScalarField<M, decltype([](const Vec& p) {
        double x = p(0);
        double y = p(1);

        return  8 * x * y * (8 * x*x + 8 * y*y - 15 ) ;
    })> u; // forcing term


    for (const auto& r : ref){

        
        std::string path = "../../plots/data/" + folder + "/";
    
        auto mesh = load_mesh(path);
    
        
        
        mesh.refine_knots({r,r});

        // print the number of cells
        std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
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
    
        
    
        std::string save_path = "../curly_plate_ref" + std::to_string(r) + "/"; // or wherever you want
        export_mesh(mesh, save_path);
        
        IsoSpace Vh(mesh);
    
        TrialFunction f(Vh);
        TestFunction v(Vh);
        
    
        //auto start = std::chrono::high_resolution_clock::now();
        auto a = integral(mesh,QGL2DP4)(dot(grad(f), grad(v))); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
        //auto m = integral(mesh, QGL2DP9)(v);
        auto F = integral(mesh,QGL2DP4)(u*v);
        //auto mm = integral(mesh,QGL2DP9)(f*v);
    
    
        auto& dof_handler = Vh.dof_handler();

        auto start = std::chrono::high_resolution_clock::now();

        Eigen::SparseMatrix<double> A = a.assemble();
        //auto end = std::chrono::high_resolution_clock::now();
        //std::cout << "Assembly time for A: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
        //start = std::chrono::high_resolution_clock::now();
        auto b = F.assemble();   
        //end = std::chrono::high_resolution_clock::now();
        //std::cout << "Assembly time for b: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;  
        //start = std::chrono::high_resolution_clock::now();
        dof_handler.set_hom_dirichlet_constraint();
        dof_handler.enforce_constraints(A,b);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Assembly time: " << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
        file2 << dof_handler.n_dofs() << "," << (end - start).count() << std::endl;

    
        
        
        //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        //solver.compute(A);
        //auto uh_full = solver.solve(b);
    
    
        //std::cout << "uh_full: " << uh_full.transpose() << std::endl;
    
        //IsoFunction solution(Vh);
        //solution =  uh_full;
    
        //export_results(mesh, solution, folder, std::make_optional(f_exact), /*nn = */ 30); //std::nullopt if there is no solution 

    }
        
        

    return 0;
}