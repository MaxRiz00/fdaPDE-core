#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;
using vector_t = Eigen::Matrix<double, Dynamic, 1>;
using matrix_t = Eigen::SparseMatrix<double>;

int main(){

    auto mesh =  IsoMesh<2,3>::sphere();
    mesh.refine_knots({3,3});

    std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
    // print the knots
    std::cout << "Knots: " << std::endl;
    for (int d = 0; d < 2; d++) {
        std::cout << "Dimension " << d << ": ";
        for (const auto& knot : mesh.knots()[d]) {
            std::cout << knot << " ";
        }
        std::cout << std::endl;
    }

    std::string save_path = "../sphere/"; // or wherever you want
    std::string folder = "sphere";
    export_mesh(mesh, save_path);
    
    IsoSpace Vh(mesh);

    TrialFunction f(Vh);
    TestFunction v(Vh);

    
    /*
    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        return 1;
    })>u;
    */

    /*
    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) { return sin(p[0]) * sin(p[1]) * sin(p[2]); })> u;
    */

    //(x-x0) * (y-y0) * (y-y0) - (y-y0) * (z-z0)* (z-z0) + (x-x0) * (x-x0) * (z-z0);
    
    
    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);

        double r = std::sqrt(x * x + y * y + z * z);
        double phi = std::atan2(y, x);
        double theta = std::acos(z / r );  

        return 2 * std::sin(theta) * std::sin( phi);
        })> u;

        // evaluation of u at poles
   

    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);

        double r = std::sqrt(x * x + y * y + z * z);
        double phi = std::atan2(y, x);
        double theta = std::acos(z / r );  

        return std::sin(theta) * std::sin( phi);
        })> u_exact;

        

    

    /*
    constexpr double alpha = 4.; // 3.
    constexpr double beta = 4.0; //
    
    

    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);

        double r = std::sqrt(x * x + y * y + z * z);
        double phi = std::atan2(y, x);
        double theta = std::acos(z / r );  

        return std::sin(alpha * phi) * std::sin(beta * theta) * (
            alpha * alpha / (std::sin(theta) * std::sin(theta) + 1e-10 ) +
            beta * beta -
            beta * (std::cos(theta) * std::cos(beta * theta)) / (std::sin(theta) * std::sin(beta * theta) + 1e-10)
        );
        })> u;

        // evaluation of u at poles
        std::cout<<"ECCO: " << u(Eigen::Matrix<double, 3, 1>(0,0,1))<<std::endl;
        std::cout<<"ECCO: " << u(Eigen::Matrix<double, 3, 1>(0,0,-1))<<std::endl;


        ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
            double x = p(0);
            double y = p(1);
            double z = p(2);
            double r = std::sqrt(x * x + y * y + z * z);
            double phi = std::atan2(y, x);
            double theta = std::acos(z / r );  
            return std::sin(alpha * phi) * std::sin(beta * theta);
            })> u_exact;
    */




    auto start = std::chrono::high_resolution_clock::now();
    auto a = integral(mesh,QGL2DP9)(dot(grad(f), grad(v))); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
    auto m = integral(mesh, QGL2DP9)(v); // for the boundary conditions
    auto F = integral(mesh,QGL2DP9)(u*v);
    //auto mm = integral(mesh,QGL2DP9)(f*v);


    auto& dof_handler = Vh.dof_handler();
    
    // start the timer
    
    Eigen::SparseMatrix<double> A = a.assemble();
    //auto M = mm.assemble();
    auto c = m.assemble();
    auto b = F.assemble();

    // Reducing the system using the periodic BC

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
    

    Eigen::SparseMatrix<double> Zero(1, 1);
    SparseBlockMatrix<double, 2, 2> D(A_reduced, c_reduced.sparseView(), c_reduced.transpose().sparseView(), Zero);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(D);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(counter + 1);
    rhs.head(counter) = b_reduced;
    Eigen::VectorXd uh_reduced = solver.solve(rhs).head(counter);

    std::cout << "uh_reduced: " << uh_reduced.transpose() << std::endl;
    std::cout << "uh_reduced size: " << uh_reduced.size() << std::endl;
    

    /*
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_reduced);
    Eigen::VectorXd uh_reduced = solver.solve(b_reduced);
    */

    Eigen::VectorXd uh_full(dof_map.size());
    for (int i = 0; i < dof_map.size(); ++i) {
        int mapped = dof_map[i];
        uh_full[i] = uh_reduced[reduced_indices[mapped]];
    }
    

    
    std::cout << "uh_full: " << uh_full.transpose() << std::endl;

    IsoFunction solution(Vh);
    solution =  uh_full.topRows(A.rows());

    int nn = 5; // number of evaluation per cell for each dimension (for plot purposes)

    export_results(mesh, solution, folder, std::make_optional(u_exact), nn); //std::nullopt if there is no solution 

    return 0;
}