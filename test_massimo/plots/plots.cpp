#include "isogeometric.h"
#include <unsupported/Eigen/SparseExtra>
#include <fstream>

template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;

// this function is use to plot the geometry
using namespace fdapde;

int main(){

    // order
    std::string folder = "torus";
    std::string path = "../data/" + folder + "/";

    SpMatrix<double> knots_x, knots_y, weights;
    SpMatrix<double> control_points_x, control_points_y, control_points_z;

    Eigen::loadMarket(knots_x, path + "knots_x.mtx");
    Eigen::loadMarket(knots_y, path + "knots_y.mtx");

    Eigen::loadMarket(weights, path + "weights.mtx");

    Eigen::loadMarket(control_points_x, path + "ctrlpts_x.mtx");
    Eigen::loadMarket(control_points_y, path + "ctrlpts_y.mtx");
    Eigen::loadMarket(control_points_z, path + "ctrlpts_z.mtx");

    std::array<std::vector<double>,2> nodes;
    MdArray<double, full_dynamic_extent_t<2>> weights_(weights.rows(), weights.cols());
    MdArray<double, full_dynamic_extent_t<3>> control_points(control_points_x.rows(), control_points_x.cols(), 3);

    nodes[0].resize(knots_x.cols());
    nodes[1].resize(knots_y.cols());

    for(size_t i = 0; i < knots_x.cols(); i++){
        nodes[0][i] = knots_x.coeff(0, i);
    }

    for(size_t i = 0; i < knots_y.cols(); i++){
        nodes[1][i] = knots_y.coeff(0, i);
    }

    for(size_t i = 0; i < weights.rows(); i++){
        for(size_t j = 0; j < weights.cols(); j++){
            weights_(i, j) = weights.coeff(i, j);
        }
    }

    for(size_t i = 0; i < control_points_x.rows(); i++){
        for(size_t j = 0; j < control_points_x.cols(); j++){
            control_points(i, j, 0) = control_points_x.coeff(i, j);
            control_points(i, j, 1) = control_points_y.coeff(i, j);
            control_points(i, j, 2) = control_points_z.coeff(i, j);
        }
    }

    

    int order = 2;

    IsoMesh<2, 3> mesh(nodes, weights_, control_points, order);

    // size of knots
    std::cout<<"# knots sizes : "<<mesh.knots()[0].size()<<" "<<mesh.knots()[1].size()<<std::endl;

    
    // Export physical nodes
    std::ofstream nodes_file("../nodes.txt");
    std::ofstream edges_file("../edges.txt");

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

    std::ofstream boundary_edges_file("../boundary_edges.txt");

    for(auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it){
        if(it->on_boundary()) boundary_edges_file<<1<<std::endl;
        else boundary_edges_file<<0<<std::endl;
    }

    // save the control points
    std::ofstream control_points_file("../control_points.txt");
    for(int i = 0; i < control_points.extent(0); i++){
        for(int j = 0; j < control_points.extent(1); j++){
            control_points_file<<control_points(i, j, 0)<<" "<<control_points(i, j, 1)<<" "<<control_points(i, j, 2)<<std::endl;
        }
    }

    // compute a linspace in the parametric domain
    int N =100;

    double step_x = (nodes[0].back() - nodes[0].front())/(N);
    double step_y = (nodes[1].back() - nodes[1].front())/(N);

    // Save the points in a file, to be read from R for the plots
    std::ofstream result("physical_points.txt");


    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            double x = nodes[0].front() + i*step_x;
            double y = nodes[1].front() + j*step_y;
            auto p = mesh.eval_param({x, y});
            result<< p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }

}