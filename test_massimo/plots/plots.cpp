#include "isogeometric.h"
#include <unsupported/Eigen/SparseExtra>
#include <fstream>

template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;

// this function is use to plot the geometry
using namespace fdapde;

int main(){
    /*

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
    */

    int order = 1;

    std::array<std::vector<double>,3> knots ;
    MdArray<double,full_dynamic_extent_t<3>> weights(2,3,2);
    MdArray<double,full_dynamic_extent_t<4>> control_points(2,3,2,3);

    knots[0].resize(2);
    knots[1].resize(3);
    knots[2].resize(2);

    for(int i = 0; i<2; i++) knots[0][i] = 1.*i;
    for(int i = 0; i<3; i++) knots[1][i] = 0.5*i;
    for(int i = 0; i<2; i++) knots[2][i] = 1.*i;

    for(int i = 0; i<2; i++) for(int j = 0; j<3; j++) for(int k = 0; k<2; k++) weights(i,j,k) = 1.;

     for(size_t i = 0; i < 3; i++){
        for(size_t j = 0; j < 2; j++){
            control_points(0,i,j,0) = (i<2)?-1.:1.;
            control_points(0,i,j,1) = (i<1)?-1.:1.;
            control_points(0,i,j,2) = (j<1)? 0.:1.;
            control_points(1,i,j,0) = (i<2)? 0.:1.;
            control_points(1,i,j,1) = (i<1)?-1.:0.;
            control_points(1,i,j,2) = (j<1)? 0.:1.;
        }
    }
    IsoMesh<3,3> mesh2(knots, weights, control_points, order);
    /*
    std::cout<<"# faces: "<<mesh2.n_boundary_faces()<<std::endl;
    for(auto it = mesh2.boundary_faces_begin(); it!=mesh2.boundary_faces_end(); ++it){
        std::cout<<"Processing face "<<it->id()<<std::endl;
        std::cout<<"Node indices: "<<std::endl;
        for(int i = 0; i < mesh2.n_nodes_per_face; i++){
            std::cout<<it->node_ids()[i]<<" ";
        }
        std::cout<<std::endl;
        // check if the edge is on the boundary
        if(it->on_boundary()){
            std::cout<<"Face is on the boundary"<<std::endl;
        }
        
    }
    */

    // nuber of cells
    std::cout<<"Number of cells: "<<mesh2.n_cells()<<std::endl;

    // print the neighbors of the cells
    auto neighbors = mesh2.neighbors();
    for(int i = 0; i < neighbors.rows(); i++){
        for(int j = 0; j < neighbors.cols(); j++){
            std::cout<<neighbors(i,j)<<" ";
        }
        std::cout<<std::endl;
    }

    //for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it){
        //std::cout << it->id() << std::endl;
   // }

    // loop over the edgesss
    /*
    for(auto it = mesh.edges_begin(); it != mesh.edges_end(); ++it){
        std::cout<<"Processing edge "<<it->id()<<std::endl;
        std::cout<<"Node indices: "<<std::endl;
        for(int i = 0; i < mesh.n_nodes_per_edge; i++){
            std::cout<<it->node_ids()[i]<<" ";
        }
        std::cout<<std::endl;
        // check if the edge is on the boundary
        if(it->on_boundary()){
            std::cout<<"Edge is on the boundary"<<std::endl;
        }
        
    }
    */
    /*
    for(auto it = mesh.boundary_edges_begin(); it != mesh.boundary_edges_end(); ++it){
        std::cout<<"Processing edge "<<it->id()<<std::endl;
        std::cout<<"Node indices: "<<std::endl;
        for(int i = 0; i < mesh.n_nodes_per_edge; i++){
            std::cout<<it->node_ids()[i]<<" ";
        }
        std::cout<<std::endl;
        // check if the edge is on the boundary
        if(it->on_boundary()){
            std::cout<<"Edge is on the boundary"<<std::endl;
        }
        
    }
    */

    // print the number of boundary edges
    //std::cout<<"Boundary edges: "<<mesh.n_boundary_edges()<<std::endl;
    
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