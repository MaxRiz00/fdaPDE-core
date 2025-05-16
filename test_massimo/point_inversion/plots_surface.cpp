#include "isogeometric.h"

#include <unsupported/Eigen/SparseExtra>
#include <fstream>

template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;

// this function is use to plot the geometry
using namespace fdapde;

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> array;
    double step = (b - a) / (n - 1);

    while(a <= b) {
        array.push_back(a);
        a += step;
    }

    return array;
}


// Utility: get a sorted pair for edge key
std::pair<int, int> sorted_edge(int a, int b) {
    return std::minmax(a, b);
}

int main(){

    // order
    std::array<int,2> order;

    std::string folder = "curved";
    std::string path = "../../plots/data/" + folder + "/";

    SpMatrix<double> knots_x, knots_y, weights;
    SpMatrix<double> control_points_x, control_points_y, control_points_z;

    // load the order
    std::ifstream order_file(path + "order.txt");
    for(int i = 0; i < 2; i++){
        order_file>>order[i];
        std::cout<<"Order "<<i<<" : "<<order[i]<<std::endl;
    }

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


    IsoMesh<2, 3> mesh(nodes, weights_, control_points, order);
    //auto mesh = IsoMesh<2, 3>::sphere(1.0);
    std::cout<<"# knots sizes : "<<mesh.knots()[0].size()<<" "<<mesh.knots()[1].size()<<std::endl;

    // print the knots
    for(int i = 0; i < 2; i++){
        std::cout<<"Knots "<<i<<" : ";
        for(auto k : mesh.knots()[i]){
            std::cout<<k<<" ";
        }
        std::cout<<std::endl;
    }

    /*
    
    // linspace of the parametric domain
    std::array<std::vector<double>, 3> test_u;
    test_u[0] = linspace(mesh.knots()[0][0] + 1e-5, mesh.knots()[0][mesh.knots()[0].size()-1]-1e-5, 80);
    test_u[1] = linspace(mesh.knots()[1][0]+ 1e-5, mesh.knots()[1][mesh.knots()[1].size()-1]-1e-5, 80);

    std::vector<Eigen::Matrix<double, 3, 1>> test_P;
    std::vector<Eigen::Matrix<double, 2, 1>> test_n;
    

    for(auto u : test_u[0]){
        for(auto v : test_u[1]){
                test_n.push_back({u, v});
                test_P.push_back(mesh.eval_param({u, v}));
        }
    }

    std::cout<<"# test points: "<<test_P.size()<<std::endl;

    // tic
    
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::Matrix<double, 3, 1> P ;
    double t1_sum = 0, t2_sum = 0;
    for(int h=0;h<test_P.size();h++){
        double t1, t2;
        P = test_P[h];
        auto u = test_n[h];
        auto u_ = mesh.invert_point(P,t1,t2);
        t1_sum += t1;
        t2_sum += t2;
        //std::cout<<"P: "<<P.transpose()<<" u: "<<u.transpose()<<" u_: "<<u_.transpose()<<std::endl;
    }

    std::cout<<"t1 mean: "<<t1_sum/test_P.size()<<" t2 mean: "<<t2_sum/test_P.size()<<std::endl;
    */
     
    double t1,t2;
    Eigen::Matrix<double, 3, 1> P(-22.1458,  22.1458,  0.14152);
    Eigen::Matrix<double, 2, 1> u_true(2.9, 2.9);

    std::cout<<"Ecco:"<<mesh.eval_param(u_true).transpose()<<std::endl;

    auto u = mesh.invert_point(P,t1,t2);
    std::cout<<"P: "<<P.transpose()<<" u: "<<u.transpose()<<std::endl;
    
    

    // toc  
    //auto finish = std::chrono::high_resolution_clock::now();

    // print the time
    //std::chrono::duration<double> elapsed = finish - start;
    //std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        
    //ScalarField<3, decltype([](const Eigen::Matrix<double, 3, 1>&) { return 1; })> one;
    //auto a = integral(mesh, QGL2DP9)(one); // physical measure of mesh
    

    //std::cout<<"Ecco a: "<<a<<std::endl;

    // refine
    //mesh.refine_knots({2,0});
    

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
    int n = 30;
    for(auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it){
        auto eval = it->evaluation(n);
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
        MdArray<double, full_dynamic_extent_t<3>> vals;
        MdArray<double, full_dynamic_extent_t<3>> eval = it->linspace_evaluation(n,vals);
        for(int i = 0; i < eval.extent(0); i++){
            for(int j = 0; j < eval.extent(1); j++){
                    evaluation_file << it->id() << "," << i << "," << j << ","
                    << eval(i,j,0) << ","
                    << eval(i,j,1) << ","
                    << eval(i,j,2) << "\n";
            }
        }
    }

    return 0;
}