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


    std::string folder = "sphere/";
    //std::string path = "../../plots/data/" + folder + "/";

    //auto mesh = load_mesh(path);
    auto mesh =  IsoMesh<2,3>::sphere();

    // print the number of cells
    std::cout << "Number of cells: " << mesh.n_cells() << std::endl;
    //mesh.refine_knots({0,1});
    /*

    // print the knots
    std::cout << "Knots: " << std::endl;
    for (int d = 0; d < 2; d++) {
        std::cout << "Dimension " << d << ": ";
        for (const auto& knot : mesh.knots()[d]) {
            std::cout << knot << " ";
        }
        std::cout << std::endl;
    }
        */

    //std::string save_path = "../sph/"; // or wherever you want
    //export_mesh(mesh, save_path);


    mesh.refine_knots({2,2});
    
    IsoSpace Vh(mesh);

    TrialFunction f(Vh);
    TestFunction v(Vh);

    ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) { return 1; })> u;

/*

    constexpr double alpha = 3.0;
constexpr double beta = 4.0;

ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>& p) {
    double x = p(0);
    double y = p(1);
    double z = p(2);

    double r = std::sqrt(x * x + y * y + z * z);
    double theta = std::atan2(y, x);
    double phi = std::acos(z / (r + 1e-10));  // small epsilon to avoid division by zero
    double sin_phi = std::sin(phi);
    double beta_phi = beta * phi;
    double sin_beta_phi = std::sin(beta_phi);
    double cos_phi = std::cos(phi);
    double cos_beta_phi = std::cos(beta_phi);

    double term1 = alpha * alpha / (sin_phi * sin_phi + 1e-10);
    double term2 = beta * beta;
    double term3 = beta * cos_phi * cos_beta_phi / (sin_phi * sin_beta_phi + 1e-10);

    return std::sin(alpha * theta) * sin_beta_phi * (term1 + term2 - term3);
    })> u;
    */
    
    

    auto start = std::chrono::high_resolution_clock::now();
    auto a = integral(mesh,QGL2DP4)(dot(grad(f),grad(v))); // dot(grad(f), grad(v)) laplacian(f)*laplacian(v)
    auto m = integral(mesh, QGL2DP4)(f * v);
    auto F = integral(mesh,QGL2DP4)(u*v);


    auto& dof_handler = Vh.dof_handler();
    
    // start the timer
    
    Eigen::SparseMatrix<double> A = a.assemble();
    Eigen::SparseMatrix<double> M = m.assemble();
    auto b = F.assemble();

    // save A,M,b in a file
    
    std::ofstream A_file( "../A.txt");
    A_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(A) << std::endl;
    std::ofstream M_file("../M.txt");
    M_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(M) << std::endl;

    std::ofstream b_file("../b.txt");
    b_file <<std::fixed << std::setprecision(10)<< Eigen::MatrixXd(b) << std::endl;
    

    //dof_handler.set_hom_dirichlet_constraint();
    //dof_handler.enforce_constraints(A,b);

    

    
    // periodic BC
    
    

    std::cout<<"A.rows(): "<<A.rows()<<std::endl;
    std::cout<<"A.cols(): "<<A.cols()<<std::endl;

    std::cout<<"Periodic DOF mapping"<<std::endl;

    auto n_dofs = dof_handler.n_dofs();
    std::vector<int> dof_map(n_dofs);
    for(int i = 0; i < n_dofs; i++){
        dof_map[i] = i;
    }

    for(int d = 0; d < 2; d++){
        std::cout<<"Periodic 0: "<<mesh.is_periodic(d)<<std::endl;
        if(mesh.is_periodic(d)){
            std::vector<int> dofs_min, dofs_max;

            dof_handler.get_boundary_dofs_for_dimension(d, true, dofs_min);
            dof_handler.get_boundary_dofs_for_dimension(d, false, dofs_max);

            assert(dofs_min.size() == dofs_max.size());

            for(int k = 0; k<dofs_max.size(); k++){
                dof_map[dofs_max[k]] = dofs_min[k];
            }
        }
    }

    std::cout<<"Applying periodic mapping"<<std::endl;

    auto map_matrix = [&](const SpMatrix<double>& mat) {
        std::vector<Eigen::Triplet<double>> triplets;
        for (int k = 0; k < mat.outerSize(); ++k) {
            for (SpMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                int i = dof_map[it.row()];
                int j = dof_map[it.col()];
                triplets.emplace_back(i, j, it.value());
            }
        }
        SpMatrix<double> mapped(mat.rows(), mat.cols());
        mapped.setFromTriplets(triplets.begin(), triplets.end());
        return mapped;
    };

    Eigen::VectorXd b_mapped = Eigen::VectorXd::Zero(b.size());
    for (int i = 0; i < b.size(); ++i) {
        b_mapped[dof_map[i]] += b[i];
    }
    
    auto A_mapped = map_matrix(A);
    auto M_mapped = map_matrix(M);

    // save A_mapped and M_mapped in a file (make them full matrices)
    std::ofstream A_mapped_file(folder + "A_mapped.txt");
    A_mapped_file << Eigen::MatrixXd(A_mapped) << std::endl;

    std::ofstream M_mapped_file(folder + "M_mapped.txt");
    M_mapped_file << Eigen::MatrixXd(M_mapped) << std::endl;

    std::cout<<"Removing duplicate dofs"<<std::endl;


    std::vector<int> keep_dofs;
    for (size_t i = 0; i < dof_map.size(); ++i) {
        if (dof_map[i] == i) keep_dofs.push_back(i);
    }
    std::cout<<"Keep dofs size: "<<keep_dofs.size()<<std::endl;

    auto filter_matrix = [&](const SpMatrix<double>& mat) {
        std::unordered_map<int, int> global_to_local;
        for (size_t i = 0; i < keep_dofs.size(); ++i) {
            global_to_local[keep_dofs[i]] = static_cast<int>(i);
        }
    
        std::vector<Eigen::Triplet<double>> triplets;
        Eigen::SparseMatrix<double> filtered(keep_dofs.size(), keep_dofs.size());
    
        for (int k = 0; k < mat.outerSize(); ++k) {
            for (SpMatrix<double>::InnerIterator it(mat, k); it; ++it) {
                int row = it.row();
                int col = it.col();
                auto row_it = global_to_local.find(row);
                auto col_it = global_to_local.find(col);
                if (row_it != global_to_local.end() && col_it != global_to_local.end()) {
                    triplets.emplace_back(row_it->second, col_it->second, it.value());
                }
            }
        }
    
        filtered.setFromTriplets(triplets.begin(), triplets.end());
        return filtered;
    };

    auto A_reduced = filter_matrix(A_mapped);
    auto M_reduced = filter_matrix(M_mapped);


    // save A_mapped and M_mapped in a file (make them full matrices)
    std::ofstream A_red_file(folder + "A_reduced.txt");
    A_red_file << Eigen::MatrixXd(A_reduced) << std::endl;

    std::ofstream M_red_file(folder + "M_reduced.txt");
    M_red_file << Eigen::MatrixXd(M_reduced) << std::endl;

    Eigen::VectorXd b_reduced(keep_dofs.size());
    for (size_t i = 0; i < keep_dofs.size(); ++i) {
        b_reduced[i] = b_mapped[keep_dofs[i]];
    }

    std::ofstream b_red_file(folder + "b_reduced.txt");
    b_red_file << b_reduced << std::endl;

    std::cout << "Adding zero mean constraint (mass matrix aware)" << std::endl;

    //std::cout<<"A_reduced.rows(): "<<A_reduced<<std::endl;

    int n = A_reduced.rows();

    // Compute the mass-row: M_reduced * 1
    Eigen::VectorXd one_vec = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd mass_row = M_reduced * one_vec;

    // Augment system matrix
    Eigen::MatrixXd A_aug(n + 1, n);
    A_aug.topRows(n) = Eigen::MatrixXd(A_reduced);             // n x n
    A_aug.row(n) = mass_row.transpose();                       // 1 x n
    
    // Augment RHS
    Eigen::VectorXd b_aug(n + 1);
    b_aug.head(n) = b_reduced;
    b_aug[n] = 0.0;
    
    // Solve least-squares (rectangular system)
    Eigen::VectorXd uh_reduced = A_aug.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_aug);

    Eigen::VectorXd uh_full(dof_map.size());
    uh_full.setZero();

    for (size_t i = 0; i < keep_dofs.size(); ++i) {
        uh_full[keep_dofs[i]] = uh_reduced[i];
    }

    // Apply mapped values to periodic DOFs
    for (size_t i = 0; i < dof_map.size(); ++i) {
        if (dof_map[i] != i) {
            uh_full[i] = uh_full[dof_map[i]];
        }
    }


    std::cout<<"Uh reduced: "<<uh_reduced.transpose()<<std::endl;
    std::cout<<"Uh size: "<<uh_reduced.size()<<std::endl;

    std::cout<<"Uh full: "<<uh_full.transpose()<<std::endl;
    std::cout<<"Uh full size: "<<uh_full.size()<<std::endl;
    

    
    
    //Eigen::SparseLU<Eigen::SparseMatrix<double>> invA(A);
    //auto uh_full = invA.solve(b);

    //stop the timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

    int nn = 30; 
    
    
    IsoFunction solution(Vh);
    solution =  uh_full;
    std::cout<<"U full: "<<uh_full<<std::endl;
    

    Eigen::Matrix<double, 2, 1> param_point;
    param_point << 0.5, 0.5;

    //std::cout << "solution:\n" << solution(param_point) << std::endl;

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
        for(int i = 0; i < eval.extent(0); i++){
            for(int j = 0; j < eval.extent(1); j++){
                    Eigen::Matrix<double, 2, 1> param_point;
                    param_point << param_points(i, j, 0), param_points(i, j, 1);
                    //std::cout<<"param_point: "<<param_point.transpose()<<std::endl;
                    evaluation_file << it->id() << "," << i << "," << j << ","
                    << eval(i,j,0) << ","
                    << eval(i,j,1) << ","
                    << eval(i,j,2) << ","
                    << solution(param_point)<< "\n";
            }
        }
    }
        
        

    return 0;
}