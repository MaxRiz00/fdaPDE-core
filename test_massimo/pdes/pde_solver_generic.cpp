#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;


    std::string folder = "torus/";
    std::string path = "../../plots/data/" + folder + "/";

    //auto mesh = load_mesh(path);
    auto mesh =  IsoMesh<2,3>::torus();

    // print the number of cells
    std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
    mesh.refine_knots({1,1});
    
    //mesh.refine_knots({0,0});
    // print the knots
    std::cout << "Knots: " << std::endl;
    for (int d = 0; d < 2; d++) {
        std::cout << "Dimension " << d << ": ";
        for (const auto& knot : mesh.knots()[d]) {
            std::cout << knot << " ";
        }
        std::cout << std::endl;
    }

    

    std::string save_path = "../torus/"; // or wherever you want
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
    
    /*
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

        */

    

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

    

    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);

        int m = 2;
        int n = 3;

        double phi = std::atan2(y, x);
        double theta = std::atan2(z, std::sqrt(x*x + y*y) - 2);
    
        return std::sin(m * theta) * std::sin(n * phi);
    })> u_exact;
    
    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);

        double R = 2;
        double r = 1;

        int m = 2;
        int n = 3;

    
        double phi = std::atan2(y, x);
        double theta = std::atan2(z, std::sqrt(x*x + y*y) - 2);

        return (m * std::sin(n * phi) * std::cos(m * theta) * std::sin(theta) / (r * (R + r * std::cos(theta) ))) +
                std::sin(n * phi) * std::sin(m * theta) * (m * m/(r * r) + n * n / ((R + r * std::cos(theta) ) * (R + r * std::cos(theta) )));
        
        //std::sin(phi) * std::sin(theta)* (2  * (1 + std::cos(theta)) / (2 + std::cos(theta)) + 1/((2 + std::cos(theta)) * (2 + std::cos(theta))) );
    })> u;
    


    auto start = std::chrono::high_resolution_clock::now();
    auto a = integral(mesh,QGL2DP9)(dot(grad(f), grad(v))); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
    auto m = integral(mesh, QGL2DP9)(v);
    auto F = integral(mesh,QGL2DP9)(u*v);
    //auto mm = integral(mesh,QGL2DP9)(f*v);


    auto& dof_handler = Vh.dof_handler();
    
    // start the timer
    
    Eigen::SparseMatrix<double> A = a.assemble();
    //auto M = mm.assemble();
    auto c = m.assemble();
    auto b = F.assemble();


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

    uh_reduced.setZero();
    uh_reduced(11) = 1;
    

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
        




    // save A,M,b in a file

    
    //dof_handler.set_hom_dirichlet_constraint();
    //dof_handler.set_periodic_constraint3();
    //dof_handler.enforce_constraints(A,b);
    //dof_handler.enforce_constraints(M);
    //dof_handler.enforce_constraints(c);
    //Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(M.rows());
    //Eigen::VectorXd c = M * one_vec;
    // remove 
    //std::cout<<"C: "<<c<<std::endl;

    

            /*
    // print the determinant of A
    std::cout<<"Determinant of A: "<<A.toDense().determinant()<<std::endl;
    
    std::ofstream A_file( "../A.txt");
    A_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(A) << std::endl;
    std::ofstream M_file("../c.txt");
    M_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(c) << std::endl;
    std::ofstream b_file("../b.txt");
    b_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(b) << std::endl;
    */


    // imposition of integral constraints
    
    /*
    matrix_t Zero(1, 1);
    SparseBlockMatrix<double, 2, 2> D(A, c.sparseView(), c.sparseView().transpose(), Zero);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> invD;
    invD.compute(D);
    vector_t rhs = vector_t::Zero(A.rows()+1);
    rhs.topRows(A.rows()) = b;
    auto uh_full = invD.solve(rhs);

    double constraint_val = c.dot(uh_full.head(A.rows()));
    std::cout << "Constraint value: " << constraint_val << std::endl;

    */
    
    
    


    /*
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    auto uh_full = solver.solve(b);
    */
    
    
    
    

    
    std::cout << "uh_full: " << uh_full.transpose() << std::endl;

    IsoFunction solution(Vh);
    solution =  uh_full.topRows(A.rows());
    /*
    Eigen::Matrix<double,Dynamic,1> sol;
    sol.resize(uh_full.topRows(A.rows()).size());
    sol.setZero();
    sol(13) = 1;
    sol(17) = 1;
    solution = sol;
    */
    int nn = 5;
    

    // create a result folder
    std::string result_folder = "../results/" + folder + "/";
    std::string command = "mkdir -p " + result_folder;
    system(command.c_str());

    // Export physical nodes
    std::ofstream nodes_file(result_folder + "nodes.txt");
    std::ofstream edges_file(result_folder + "edges.txt");

    std::cout<<"# nodes: "<<mesh.n_nodes()<<std::endl;

    for(int i = 0; i < mesh.n_nodes(); i++){
        //std::cout<<"Processing node "<<i<<std::endl;
        auto p = mesh.phys_node(i);
        nodes_file << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }

    auto edges = mesh.edges();

    std::cout<<"# edges: "<<edges.rows()<<std::endl;

    for(int i = 0; i < edges.rows(); i++){
        for(int j = 0; j < edges.cols(); j++){
            edges_file<<edges(i,j)<<" ";
        }
        edges_file<<std::endl;
    }

    // iterate over edges and print the nodes
    std::ofstream boundary_edges_file(result_folder + "boundary_edges.txt");

    for(auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it){
        if(it->on_boundary()) boundary_edges_file<<1<<std::endl;
        else boundary_edges_file<<0<<std::endl;
    }
        

    // save the control points
    std::ofstream control_points_file(result_folder + "control_points.txt");
    auto control_pointsm = mesh.control_points();
    // fist line the dimension of cp
    control_points_file<<control_pointsm.extent(0)<<" "<<control_pointsm.extent(1)<<std::endl;
    for(int i = 0; i < control_pointsm.extent(0); i++){
        for(int j = 0; j < control_pointsm.extent(1); j++){
            control_points_file<<control_pointsm(i, j, 0)<<" "<<control_pointsm(i, j, 1)<<" "<<control_pointsm(i, j, 2)<<std::endl;
        }
    }

    // save the weights
    std::ofstream weights_file(result_folder + "weights.txt");
    auto weightsm = mesh.weights();
    for(int i = 0; i < weightsm.extent(0); i++){
        for(int j = 0; j < weightsm.extent(1); j++){
            weights_file<<weightsm(i, j)<<std::endl;
        }
    }

    std::cout<<"Computing the edge refinement"<<std::endl;
    std::ofstream edge_refinement_file(result_folder + "edge_refinement.txt");
    // iterate over the edges
    for(auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it){
        auto eval = it->evaluation(nn);
        for(int i = 0; i < eval.rows(); i++){
            edge_refinement_file<<eval(i, 0)<<" "<<eval(i, 1)<<" "<<eval(i, 2)<<std::endl;
        }
    }

    std::cout<<"End edge refinement"<<std::endl;

    std::ofstream cells_file(result_folder + "cells.txt");
    Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor> cells = mesh.cells();
    for(int i = 0; i < cells.rows(); i++){
        for(int j = 0; j < cells.cols(); j++){
            cells_file<<cells(i, j)<<" ";
        }
        cells_file<<std::endl;
    }

    // evaluation for plot purposes

    std::ofstream evaluation_file(result_folder + "refined_surface_points.csv");

    // cellid x i x j x dim
    //MdArray<double, full_dynamic_extent_t<4>> evaluation;
    //evaluation.resize(cells.rows(), n, n, 3);



    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it){
        MdArray<double, full_dynamic_extent_t<3>> param_points;
        MdArray<double, full_dynamic_extent_t<3>> eval = it->linspace_evaluation(nn,param_points);
        Eigen::Matrix<double,Dynamic,Dynamic> error;
        error.resize(eval.extent(0), eval.extent(1));
        for(int i = 0; i < eval.extent(0); i++){
            for(int j = 0; j < eval.extent(1); j++){
                    Eigen::Matrix<double, 2, 1> param_point;
                    param_point << param_points(i, j, 0), param_points(i, j, 1);
                    Eigen::Matrix<double, 3, 1> p;
                    p << eval(i,j,0), eval(i,j,1), eval(i,j,2);
                    //std::cout<<"param_point: "<<param_point.transpose()<<std::endl;
                    evaluation_file << it->id() << "," << i << "," << j << ","
                    << eval(i,j,0) << ","
                    << eval(i,j,1) << ","
                    << eval(i,j,2) << ","
                    << solution(param_point)<< "\n";
                    
                    error(i,j) = solution(param_point) - u_exact(p);
                    /*
                    if(solution(param_point) > 2.){
                        std::cout<<"param_point: "<<param_point.transpose()<<std::endl;
                        std::cout<<"eval: "<<eval(i,j,0)<<" "<<eval(i,j,1)<<" "<<eval(i,j,2)<<std::endl;
                    }
                        */


            }
        }
        std::cout<<"Error: "<<error.mean()<<std::endl;
    }
        

    
        
        

    return 0;
}