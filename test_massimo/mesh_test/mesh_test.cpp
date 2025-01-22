#include "nurbs/nurbs_basis.h"
#include "nurbs/iso_mesh.h"
#include "nurbs/iso_square.h"
#include "nurbs/integrator.h"

#include <gtest/gtest.h>   // testing framework
#include <unsupported/Eigen/SparseExtra>
#include <unistd.h>  // For getcwd()

#include "linear_algebra/mdarray.h"
#include "../utils/utils.h" 

#include <chrono>


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

    // to iterate over the cells using the class CellIterator for (auto it = mesh.beginCells(); it != end; ++it)

    // Test nodes
    std::vector<std::array<double,3>> nodes = mesh.compute_nodes();
    SpMatrix<double> expected_nodes;
    Eigen::loadMarket(expected_nodes, "../mesh_test_data/nodes.mtx");
    for(size_t i = 0; i < nodes.size(); ++i){
        for(size_t j = 0; j < 3; ++j){
            EXPECT_TRUE(almost_equal(expected_nodes.coeff(i,j),nodes[i][j]));
        }
    }

    // use get_neighbors to get the neighbors of a cell, print all neighbors for ll cels

    SpMatrix<size_t> expected_neighbors;
    Eigen::loadMarket(expected_neighbors, "../mesh_test_data/neighbors.mtx");

    auto n_elements = mesh.n_elements();

    // compute the neighbors of all elements
    for(size_t i = 0; i < n_elements; ++i){
        auto neighbors = mesh.get_neighbors_ID(i);
        for(size_t j = 0; j < neighbors.size(); ++j){
            EXPECT_TRUE(neighbors[j] == expected_neighbors.coeff(i,j));
        }
    }

    // check for boundary cells
    auto n_nodes = mesh.n_nodes();
    SpMatrix<size_t> expected_boundary;
    Eigen::loadMarket(expected_boundary, "../mesh_test_data/boundary.mtx");

    bool element = false;

    for(size_t i = 0; i < n_nodes; ++i){
        // print true value and computed value
        EXPECT_TRUE(mesh.is_boundary(i,element) == expected_boundary.coeff(i,0));
    }
};


TEST(isogeometric_analysis_test, mesh_parametrization){
    
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

    SpMatrix<double> expected;
    // expected results from nurbs derivative pointwise evaluations
    Eigen::loadMarket(expected, "../mesh_test_data/nurbs_mesh_test.mtx");

    IsoMesh<3,3> mesh(knots, weights, order, control_points);

    // start measuring the time
    auto start = std::chrono::high_resolution_clock::now();

    for(size_t j = 0; j < expected.cols(); ++j){
        // first three rows of expected contain the x-y-z coordinates of the point at which to evaluate
        std::array<double, 3> x = {expected.coeff(0, j), expected.coeff(1, j), expected.coeff(2, j)};
        for(std::size_t i = 0; i<3 ; ++i){
            EXPECT_TRUE(almost_equal(expected.coeff(3+i,j),mesh.eval_param(x,i)));
            for(std::size_t k = 0; k < 3; ++k){
                EXPECT_TRUE(almost_equal(expected.coeff(6+3*i+k,j),mesh.eval_param_derivative(x,k,i)));
            }
        }
    }

    // stop measuring the time
    auto stop = std::chrono::high_resolution_clock::now();

    //print time elapsed
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

};

TEST(isogeometric_analysis_test, integrator){

    int order =3;

    std::array<std::vector<double>,3>  nodes;
    MdArray<double,full_dynamic_extent_t<3>> weights(4,4,4);
    MdArray<double,full_dynamic_extent_t<4>> controlpoints(4,4,4,3);
    nodes[0].resize(2);
    nodes[1].resize(2);
    nodes[2].resize(2);

    for(size_t i = 0; i < 2; i++)nodes[0][i]=nodes[1][i]=nodes[2][i]=1.*i;

    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 4; j++)
            for(size_t k = 0; k < 4; k++)
                weights(i,j,k) = 1.;

    // control points
    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 4; j++)
            for(size_t k = 0; k < 4; k++){
                controlpoints(i,j,k,0) = 1.*i;
                controlpoints(i,j,k,1) = 1.*j;
                controlpoints(i,j,k,2) = 1.*k;
            }
    
    IsoMesh<3,3> mesh(nodes, weights,order, controlpoints);
    IsoIntegrator<3, 27> itg;

    // exact integral of t is equal to 1
    auto f = [] (const std::array<double,3> & x) -> double 
        {return x[0]*x[0]*x[0] + 2*x[0]*x[0]*x[1] - 2*x[0]*x[1]*x[1] + 4*x[0]*x[1]*x[2] + x[2]*x[2]*x[2];};

    EXPECT_TRUE(almost_equal(1., itg.integrate(*(mesh.beginCells()), f)));
    // print the integral of the function f over the mesh
    EXPECT_TRUE(almost_equal(1., itg.integrate(mesh, f)));

    // try the physical domain integration

    auto f_physical = [] (const std::array<double,3> & x) -> double 
        {return 3*x[0]*x[0] + 2*x[0]*x[1] - 2*x[0]*x[1] + 4*x[0]*x[1] + x[2]*x[2];};

    // print the integral over the physical domain
    std::cout<<itg.integrate_physical(*(mesh.beginCells()), f_physical)<<std::endl;
}

}; // namespace fdapde::testing