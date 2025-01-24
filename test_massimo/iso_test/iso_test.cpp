// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <gtest/gtest.h>   // testing framework
#include <unsupported/Eigen/SparseExtra>
#include <unistd.h>  // For getcwd()

#include "isogeometric.h"

#include "../utils/utils.h" 


template<typename T> using SpMatrix = Eigen::SparseMatrix<T>;


namespace fdapde::testing{
// test 1D nurbs basis (functions compute the correct value)
TEST(nurbs_test, nurbs_basis_1D) {

    std::array<std::vector<double>,1> nodes;
    MdArray<double,MdExtents<Dynamic>> weights(7);
    nodes[0].resize(5);

    int order = 3;

    // open uniform knot vector
    for(size_t i = 0; i < 5; i++)nodes[0][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 7; i++)weights(i)=std::abs(std::sin(i+1));

    SpMatrix<double> expected;
    // expected results from nurbs pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_1.mtx");

    NurbsBasis<1> basis(nodes, weights, order);

    for(size_t i = 0; i < basis.size() ; ++i){ 
        for(size_t j = 0; j < expected.cols(); ++j){ 
            std::array<double, 1> input = {expected.coeff(0, j)};
            EXPECT_TRUE(almost_equal(expected.coeff(i+1,j),basis[i](input)));
        }
    } 

}

// test 2D nurbs basis (functions are accessibile and callable)
TEST(nurbs_test, nurbs_basis_2D) {
    std::array<std::vector<double>, 2> nodes;
    MdArray<double, full_dynamic_extent_t<2>> weights(4,5);
    nodes[0].resize(2);
    nodes[1].resize(3);

    int order=3;

    // open uniform knot vector
    for(size_t i = 0; i < 2; i++)nodes[0][i]=1.*i;
    for(size_t i = 0; i < 3; i++)nodes[1][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 4; i++)for(size_t j = 0; j < 5; j++)weights(i,j)=std::abs(std::sin(i+1))*std::abs(std::cos(j+1));
    

    SpMatrix<double> expected;
    // expected results from nurbs pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_2.mtx");

    NurbsBasis<2> basis(nodes, weights, order);
    
    for(size_t i = 0; i < basis.size(); ++i){
        for(size_t j = 0; j < expected.cols(); ++j){
            // compare values with data from file
            EXPECT_TRUE(almost_equal(expected.coeff(i+2,j),basis[i](Vector<double,2>(expected.coeff(0,j), expected.coeff(1,j)))));
        }
    }
}


// test 1D nurbs basis derivative (functions compute the correct value)
TEST(nurbs_test, nurbs_basis_derivative_1D){
    std::array<std::vector<double>, 1> nodes;
    MdArray<double, full_dynamic_extent_t<1>> weights(7);
    nodes[0].resize(5);

    int order = 3;

    using InputType = std::array<double, 1>;
    
    // open uniform knot vector
    for(size_t i = 0; i < 5; i++)nodes[0][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 7; i++)weights(i)=std::abs(std::sin(i+1));

    SpMatrix<double> expected;
    // expected results from nurbs derivative pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_3.mtx");

    NurbsBasis<1> basis(nodes, weights, order);
    
    for(size_t i = 0; i < basis.size(); ++i){
        for(size_t j = 0; j < expected.cols() ; ++j){
            EXPECT_TRUE(almost_equal(expected.coeff(i+1,j),basis[i].derive()(InputType{expected.coeff(0,j)})));
        }
    }
}

// test 2D nurbs basis derivative (functions are accessibile and callable)
TEST(nurbs_test, nurbs_basis_derivative_2D) {
    std::array<std::vector<double>, 2> nodes;
    MdArray<double, full_dynamic_extent_t<2>>weights(4,5);
    nodes[0].resize(2);
    nodes[1].resize(3);

    using InputType = std::array<double, 2>;

    int order = 3;

    // open uniform knot vector
    for(size_t i = 0; i < 2; i++)nodes[0][i]=1.*i;
    for(size_t i = 0; i < 3; i++)nodes[1][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 4; i++)for(size_t j = 0; j < 5; j++)weights(i,j)=std::abs(std::sin(i+1))*std::abs(std::cos(j+1));
    

    SpMatrix<double> expected;
    // expected results from nurbs derivative pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_4.mtx");

    NurbsBasis<2> basis(nodes, weights, order);
    
    for(size_t i = 0; i < basis.size(); ++i){
        for(size_t j = 0; j < expected.cols(); ++j){
            // compare values with data from file
            EXPECT_TRUE(almost_equal(expected.coeff(2*i+2,j),basis[i].derive(0)(InputType{expected.coeff(0,j), expected.coeff(1,j)})));
            EXPECT_TRUE(almost_equal(expected.coeff(2*i+3,j),basis[i].derive(1)(InputType{expected.coeff(0,j), expected.coeff(1,j)})));
        }
    }
}

// test 1D nurbs basis second derivative (functions compute the correct value)
TEST(nurbs_test, nurbs_basis_second_derivative_1D){

    using InputType = Vector<double, 1>;
    
    std::array<std::vector<double>,1> nodes;
    MdArray<double,MdExtents<Dynamic>> weights(7);
    nodes[0].resize(5);

    int order = 3;
    
    // open uniform knot vector
    for(size_t i = 0; i < 5; i++)nodes[0][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 7; i++)weights(i)=std::abs(std::sin(i+1));

    SpMatrix<double> expected;
    // expected results from nurbs derivative pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_5.mtx");

    NurbsBasis<1> basis(nodes, weights, order);
    
    for(size_t i = 0; i < basis.size(); ++i){
        for(size_t j = 0; j < expected.cols(); ++j){
            EXPECT_TRUE(almost_equal(expected.coeff(i+1,j),basis[i].deriveTwice()((InputType{expected.coeff(0,j)}))));
        }
    }
}

// test 2D nurbs basis second derivative (functions are accessibile and callable)
TEST(nurbs_test, nurbs_basis_second_derivative_2D) {

    std::array<std::vector<double>,2> nodes;
    MdArray<double,MdExtents<Dynamic,Dynamic>>weights(4,5);
    nodes[0].resize(2);
    nodes[1].resize(3);

    using InputType = Vector<double, 2>;

    int order=3;

    // open uniform knot vector
    for(size_t i = 0; i < 2; i++)nodes[0][i]=1.*i;
    for(size_t i = 0; i < 3; i++)nodes[1][i]=1.*i;
    // easily replicable, non trivial weights
    for(size_t i = 0; i < 4; i++)for(size_t j = 0; j < 5; j++)weights(i,j)=std::abs(std::sin(i+1))*std::abs(std::cos(j+1));
    

    SpMatrix<double> expected;
    // expected results from nurbs derivative pointwise evaluations
    Eigen::loadMarket(expected, "../data/mtx/nurbs_test_6.mtx");

    NurbsBasis<2> basis(nodes, weights, order);
    
    for(size_t i = 0; i < basis.size(); ++i){
        for(size_t j = 0; j < expected.cols(); ++j){
            // compare values with data from file
            EXPECT_TRUE(almost_equal(expected.coeff(4*i+2,j),basis[i].deriveTwice(0,0)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
            EXPECT_TRUE(almost_equal(expected.coeff(4*i+3,j),basis[i].deriveTwice(1,0)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
            EXPECT_TRUE(almost_equal(expected.coeff(4*i+4,j),basis[i].deriveTwice(0,1)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
            EXPECT_TRUE(almost_equal(expected.coeff(4*i+5,j),basis[i].deriveTwice(1,1)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
        }
    }
}


TEST(mesh_test, mesh_structure){

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


TEST(mesh_test, mesh_parametrization){
    
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

};

/*

TEST(integration_test, integrator){

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
*/

} // namespace fdapde::testing

