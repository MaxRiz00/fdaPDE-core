#include "nurbs/nurbs_basis.h"
#include "nurbs/iso_mesh.h"
#include "nurbs/iso_square.h"

#include <gtest/gtest.h>   // testing framework
#include <unsupported/Eigen/SparseExtra>
#include <unistd.h>  // For getcwd()

#include "linear_algebra/mdarray.h"
#include "../utils/utils.h" 


namespace fdapde::testing{

TEST(isogeometric_analysis_test, mesh_structure){

    int order = 2;

    std::array<std::vector<double>,3> knots ;
    MdArray<double,full_dynamic_extent_t<3>> weights(6,6,6);
    MdArray<double,full_dynamic_extent_t<4>> control_points(6,6,6,3);

    knots[0].resize(5);
    knots[1].resize(5);
    knots[2].resize(5);

    for(int i = 0; i<5; i++) knots[0][i] = knots[1][i] = knots[2][i] = 1.*i;

    for(int i = 0; i<6; i++){
        for(int j = 0; j<6; j++){
            for(int k = 0; k<6; k++){
                weights(i,j,k) = 1.;
                control_points(i,j,k,0) = 1.;
                control_points(i,j,k,1) = 1.;
                control_points(i,j,k,2) = 1.;
            }
        }
    }

    IsoMesh<3,3> mesh(knots, weights, order, control_points);

    // iterate over the cells using the class CellIterator for (auto it = mesh.beginCells(); it != end; ++it)

    // Test nodes
    std::vector<std::array<double,3>> nodes = mesh.compute_nodes();
    SpMatrix<double> expected_nodes;
    Eigen::loadMarket(expected_nodes, "../mesh_structure/nodes.mtx");
    for(size_t i = 0; i < nodes.size(); ++i){
        for(size_t j = 0; j < 3; ++j){
            EXPECT_TRUE(almost_equal(expected_nodes.coeff(i,j),nodes[i][j]));
        }
    }

    // use get_neighbors to get the neighbors of a cell, print all neighbors for ll cels

    auto n_elements = mesh.n_elements();

    SpMatrix<size_t> expected_neighbors;
    Eigen::loadMarket(expected_neighbors, "../mesh_structure/neighbors.mtx");

    // compute the neighbors of all elements
    for(size_t i = 0; i < n_elements; ++i){
        auto neighbors = mesh.get_neighbors_ID(i);
        for(size_t j = 0; j < neighbors.size(); ++j){
            //EXPECT_TRUE(neighbors[j] == expected_neighbors.coeff(i,j));
        }
    }

    // check for boundary cells
    auto n_nodes = mesh.n_nodes();
    SpMatrix<size_t> expected_boundary;
    Eigen::loadMarket(expected_boundary, "../mesh_structure/boundary.mtx");

    for(size_t i = 0; i < n_nodes; ++i){
        // print true value and computed value
        //std::cout << "Element " << i << " is boundary: " << mesh.is_node_boundary(i) << " expected: " << expected_boundary.coeff(i,0) << std::endl;
        EXPECT_TRUE(mesh.is_node_boundary(i) == expected_boundary.coeff(i,0));
    }


};


}; // namespace fdapde::testing