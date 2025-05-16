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
    auto order = mesh.degree();
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


template <typename ExactFieldT>
void export_results(
    IsoMesh<2, 3>& mesh,
    IsoFunction<IsoSpace<IsoMesh<2, 3>>>& solution,
    const std::string& folder,
    const std::optional<ExactFieldT>& u_exact,
    int nn = 10
) {
    std::string result_folder = "../results/" + folder + "/";
    std::string command = "mkdir -p " + result_folder;
    system(command.c_str());

    std::ofstream nodes_file(result_folder + "nodes.txt");
    for (int i = 0; i < mesh.n_nodes(); ++i) {
        auto p = mesh.phys_node(i);
        nodes_file << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    std::ofstream edges_file(result_folder + "edges.txt");
    auto edges = mesh.edges();
    for (int i = 0; i < edges.rows(); ++i) {
        for (int j = 0; j < edges.cols(); ++j) {
            edges_file << edges(i, j) << " ";
        }
        edges_file << "\n";
    }

    std::ofstream boundary_edges_file(result_folder + "boundary_edges.txt");
    for (auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
        boundary_edges_file << (it->on_boundary() ? 1 : 0) << "\n";
    }

    std::ofstream control_points_file(result_folder + "control_points.txt");
    auto cps = mesh.control_points();
    control_points_file << cps.extent(0) << " " << cps.extent(1) << "\n";
    for (int i = 0; i < cps.extent(0); ++i) {
        for (int j = 0; j < cps.extent(1); ++j) {
            control_points_file << cps(i, j, 0) << " " << cps(i, j, 1) << " " << cps(i, j, 2) << "\n";
        }
    }

    std::ofstream weights_file(result_folder + "weights.txt");
    auto weights = mesh.weights();
    for (int i = 0; i < weights.extent(0); ++i) {
        for (int j = 0; j < weights.extent(1); ++j) {
            weights_file << weights(i, j) << "\n";
        }
    }

    std::ofstream edge_ref_file(result_folder + "edge_refinement.txt");
    for (auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
        auto eval = it->evaluation(nn);
        for (int i = 0; i < eval.rows(); ++i) {
            edge_ref_file << eval(i, 0) << " " << eval(i, 1) << " " << eval(i, 2) << "\n";
        }
    }

    std::ofstream cells_file(result_folder + "cells.txt");
    auto cells = mesh.cells();
    for (int i = 0; i < cells.rows(); ++i) {
        for (int j = 0; j < cells.cols(); ++j) {
            cells_file << cells(i, j) << " ";
        }
        cells_file << "\n";
    }

    std::ofstream eval_file(result_folder + "refined_surface_points.csv");
    for (auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {
        MdArray<double, full_dynamic_extent_t<3>> param_points;
        MdArray<double, full_dynamic_extent_t<3>> eval = it->linspace_evaluation(nn, param_points);
    
        std::ofstream eval_file(result_folder + "refined_surface_points.csv", std::ios_base::app);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> error;
    
        if (u_exact.has_value()) {
            error.resize(eval.extent(0), eval.extent(1));
        }
    
        for (int i = 0; i < eval.extent(0); ++i) {
            for (int j = 0; j < eval.extent(1); ++j) {
                Eigen::Matrix<double, 2, 1> param_point;
                param_point << param_points(i, j, 0), param_points(i, j, 1);
    
                Eigen::Matrix<double, 3, 1> phys_point;
                phys_point << eval(i, j, 0), eval(i, j, 1), eval(i, j, 2);
    
                eval_file << it->id() << "," << i << "," << j << ","
                          << phys_point(0) << "," << phys_point(1) << "," << phys_point(2) << ","
                          << (*u_exact)(phys_point) << "\n";
            }
        }
    }
}