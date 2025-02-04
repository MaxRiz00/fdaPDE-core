#include "isogeometric.h"
#include <unsupported/Eigen/SparseExtra>
#include <fstream>

template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;

// this function is use to plot the geometry
using namespace fdapde;

int main(){

    SpMatrix<double> knots_x, knots_y, weights;
    SpMatrix<double> control_points_x, control_points_y, control_points_z;

    Eigen::loadMarket(knots_x, "../data/curlyplate_kx_ref3.mtx");
    Eigen::loadMarket(knots_y, "../data/curlyplate_ky_ref3.mtx");

    Eigen::loadMarket(weights, "../data/curlyplate_weights_ref3.mtx");

    Eigen::loadMarket(control_points_x, "../data/curlyplate_cpx_ref3.mtx");
    Eigen::loadMarket(control_points_y, "../data/curlyplate_cpy_ref3.mtx");
    Eigen::loadMarket(control_points_z, "../data/curlyplate_cpz_ref3.mtx");

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

    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it){
        std::cout << it->id() << std::endl;
    }
    
    /*
    
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
            std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
            result<< p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
    }
    */
    
    

}