#include "isogeometric.h"

#include <unsupported/Eigen/SparseExtra>
#include <fstream>
#include <cassert>

using namespace fdapde;

template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;

void export_mesh(const IsoMesh<2,3>& mesh, const std::string& path) {
    std::cout << "Exporting mesh to: " << path << std::endl;

    // Create directory if not exists
    std::string cmd = "mkdir -p " + path;
    system(cmd.c_str());

    // ---- 1. Order ----
    std::ofstream order_file(path + "order.txt");
    auto order = mesh.order();
    for (int d = 0; d < 2; d++) {
        order_file << order[d] << " ";
    }

// ---- 2. Knots ----
auto knots = mesh.knots();

SpMatrix<double> knots_x(1, knots[0].size());
SpMatrix<double> knots_y(1, knots[1].size());

std::vector<Eigen::Triplet<double>> triplets_kx, triplets_ky;

for (int i = 0; i < knots[0].size(); i++) {
    triplets_kx.emplace_back(0, i, knots[0][i]);
}
for (int i = 0; i < knots[1].size(); i++) {
    triplets_ky.emplace_back(0, i, knots[1][i]);
}

knots_x.setFromTriplets(triplets_kx.begin(), triplets_kx.end());
knots_y.setFromTriplets(triplets_ky.begin(), triplets_ky.end());

Eigen::saveMarket(knots_x, path + "knots_x.mtx");
Eigen::saveMarket(knots_y, path + "knots_y.mtx");

    // ---- 3. Weights ----
    auto weights = mesh.weights();
    int rows = weights.extent(0);
    int cols = weights.extent(1);
    SpMatrix<double> weights_sp(rows, cols);
    std::vector<Eigen::Triplet<double>> w_triplets;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            w_triplets.emplace_back(i, j, weights(i, j));
        }
    }
    weights_sp.setFromTriplets(w_triplets.begin(), w_triplets.end());
    Eigen::saveMarket(weights_sp, path + "weights.mtx");

    // ---- 4. Control points ----
    auto ctrlpts = mesh.control_points();
    int r = ctrlpts.extent(0);
    int c = ctrlpts.extent(1);

    SpMatrix<double> ctrl_x(r, c), ctrl_y(r, c), ctrl_z(r, c);
    std::vector<Eigen::Triplet<double>> triplets_x, triplets_y, triplets_z;

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            triplets_x.emplace_back(i, j, ctrlpts(i, j, 0));
            triplets_y.emplace_back(i, j, ctrlpts(i, j, 1));
            triplets_z.emplace_back(i, j, ctrlpts(i, j, 2));
        }
    }

    ctrl_x.setFromTriplets(triplets_x.begin(), triplets_x.end());
    ctrl_y.setFromTriplets(triplets_y.begin(), triplets_y.end());
    ctrl_z.setFromTriplets(triplets_z.begin(), triplets_z.end());

    Eigen::saveMarket(ctrl_x, path + "ctrlpts_x.mtx");
    Eigen::saveMarket(ctrl_y, path + "ctrlpts_y.mtx");
    Eigen::saveMarket(ctrl_z, path + "ctrlpts_z.mtx");

    std::cout << "Export complete!" << std::endl;
}

IsoMesh<2, 3> load_mesh(const std::string& folder_path) {
    using SpMatrix = Eigen::SparseMatrix<double>;

    std::array<int, 2> order;
    SpMatrix knots_x, knots_y, weights;
    SpMatrix control_points_x, control_points_y, control_points_z;

    std::string path = folder_path;

    // Load order
    std::ifstream order_file(path + "order.txt");
    for (int i = 0; i < 2; i++) {
        order_file >> order[i];
    }

    // Load matrices
    Eigen::loadMarket(knots_x, path + "knots_x.mtx");
    Eigen::loadMarket(knots_y, path + "knots_y.mtx");
    Eigen::loadMarket(weights, path + "weights.mtx");
    Eigen::loadMarket(control_points_x, path + "ctrlpts_x.mtx");
    Eigen::loadMarket(control_points_y, path + "ctrlpts_y.mtx");
    Eigen::loadMarket(control_points_z, path + "ctrlpts_z.mtx");

    // Reconstruct data structures
    std::array<std::vector<double>, 2> nodes;
    nodes[0].resize(knots_x.cols());
    nodes[1].resize(knots_y.cols());

    for (size_t i = 0; i < nodes[0].size(); i++) {
        nodes[0][i] = knots_x.coeff(0, i);
    }

    for (size_t i = 0; i < nodes[1].size(); i++) {
        nodes[1][i] = knots_y.coeff(0, i);
    }

    MdArray<double, full_dynamic_extent_t<2>> weights_(weights.rows(), weights.cols());
    for (int i = 0; i < weights.rows(); i++) {
        for (int j = 0; j < weights.cols(); j++) {
            weights_(i, j) = weights.coeff(i, j);
        }
    }

    MdArray<double, full_dynamic_extent_t<3>> control_points(control_points_x.rows(), control_points_x.cols(), 3);
    for (int i = 0; i < control_points_x.rows(); i++) {
        for (int j = 0; j < control_points_x.cols(); j++) {
            control_points(i, j, 0) = control_points_x.coeff(i, j);
            control_points(i, j, 1) = control_points_y.coeff(i, j);
            control_points(i, j, 2) = control_points_z.coeff(i, j);
        }
    }

    return IsoMesh<2, 3>(nodes, weights_, control_points, order);
}

int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;


    std::string folder = "torus/";
    std::string path = "../../plots/data/" + folder + "/";

    //auto mesh = load_mesh(path);
    auto mesh =  IsoMesh<2,3>::torus();

    // print the number of cells
    std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
    mesh.refine_knots({4,4});
    
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

    //ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) { return sin(2*p[0]) * sin(4*p[1]) * sin(5*p[2]); })> u;


    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);
    
        constexpr double R = 2.0;
        constexpr double r = 1.0;
        constexpr int m = 2;
        constexpr int n = 3;
    
        double theta = std::atan2(z, x);  // rotation around y-axis
        double xz_proj = std::sqrt(x * x + z * z);
        double phi = std::atan2(y, xz_proj - R);
    
        double denom = R + r * std::cos(phi);
        double factor = (n * n) / (r * r) + (m * m) / (denom * denom);
        return factor * std::sin(m * theta) * std::cos(n * phi);
    })> u;
    
    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
        double x = p(0);
        double y = p(1);
        double z = p(2);
    
        constexpr double R = 2.0;
        constexpr double r = 1.0;
        constexpr int m = 2;
        constexpr int n = 3;
    
        double theta = std::atan2(z, x);  // rotation around y-axis
        double xz_proj = std::sqrt(x * x + z * z);
        double phi = std::atan2(y, xz_proj - R);
    
        return std::sin(m * theta) * std::cos(n * phi);
    })> u_exact;



    auto start = std::chrono::high_resolution_clock::now();
    auto a = integral(mesh,QGL2DP16)(dot(grad(f), grad(v)) ); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
    auto m = integral(mesh, QGL2DP16)(v);
    auto F = integral(mesh,QGL2DP16)(u*v);
    //auto mm = integral(mesh,QGL2DP9)(f*v);


    auto& dof_handler = Vh.dof_handler();
    
    // start the timer
    
    Eigen::SparseMatrix<double> A = a.assemble();
    //auto M = mm.assemble();
    auto c = m.assemble();
    auto b = F.assemble();

    // save A,M,b in a file

    
    //dof_handler.set_hom_dirichlet_constraint();
    dof_handler.set_periodic_constraint2();
    dof_handler.enforce_constraints(A,b);
    //dof_handler.enforce_constraints(M);
    //dof_handler.enforce_constraints(c);
    //Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(M.rows());
    //Eigen::VectorXd c = M * one_vec;
    // remove 

    


    // print the determinant of A
    std::cout<<"Determinant of A: "<<A.toDense().determinant()<<std::endl;
    
    std::ofstream A_file( "../A.txt");
    A_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(A) << std::endl;
    std::ofstream M_file("../c.txt");
    M_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(c) << std::endl;
    std::ofstream b_file("../b.txt");
    b_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(b) << std::endl;


    // imposition of integral constraints
    
    matrix_t Zero(1, 1);
    SparseBlockMatrix<double, 2, 2> D(A, c.sparseView(), c.sparseView().transpose(), Zero);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> invD;
    invD.compute(D);
    vector_t rhs = vector_t::Zero(A.rows()+1);
    rhs.topRows(A.rows()) = b;
    auto uh_full = invD.solve(rhs);

    double constraint_val = c.dot(uh_full.head(A.rows()));
    std::cout << "Constraint value: " << constraint_val << std::endl;
    


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
    int nn = 10;
    

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
                    error(i,j) = solution(param_point) - u_exact(p);
                    //std::cout<<"param_point: "<<param_point.transpose()<<std::endl;
                    evaluation_file << it->id() << "," << i << "," << j << ","
                    << eval(i,j,0) << ","
                    << eval(i,j,1) << ","
                    << eval(i,j,2) << ","
                    << solution(param_point)<< "\n";
                   
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