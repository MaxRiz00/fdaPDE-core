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


#include "linear_algebra/mdarray.h"
#include "nurbs/nurbs_basis.h"
#include "../utils/utils.h" 


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
            EXPECT_TRUE(almost_equal(expected.coeff(i+2,j),basis[i](SVector<2>(expected.coeff(0,j), expected.coeff(1,j)))));
        }
    }
}


// test 1D nurbs basis derivative (functions compute the correct value)
TEST(nurbs_test, nurbs_basis_derivative_1D){
    std::array<std::vector<double>, 1> nodes;
    MdArray<double, full_dynamic_extent_t<1>> weights(7);
    nodes[0].resize(5);

    int order = 3;

    using InputType = cexpr::Vector<double, 1>;
    
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
            EXPECT_TRUE(almost_equal(expected.coeff(i+1,j),basis[i].derive()(InputType(expected.coeff(0,j)))));
        }
    }
}

// test 2D nurbs basis derivative (functions are accessibile and callable)
TEST(nurbs_test, nurbs_basis_derivative_2D) {
    std::array<std::vector<double>, 2> nodes;
    MdArray<double, full_dynamic_extent_t<2>>weights(4,5);
    nodes[0].resize(2);
    nodes[1].resize(3);

    using InputType = cexpr::Vector<double, 2>;

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
            EXPECT_TRUE(almost_equal(expected.coeff(2*i+2,j),basis[i].derive(0)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
            EXPECT_TRUE(almost_equal(expected.coeff(2*i+3,j),basis[i].derive(1)(InputType(expected.coeff(0,j), expected.coeff(1,j)))));
        }
    }
}

// test 1D nurbs basis second derivative (functions compute the correct value)
TEST(nurbs_test, nurbs_basis_second_derivative_1D){

    using InputType = cexpr::Vector<double, 1>;
    
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
            EXPECT_TRUE(almost_equal(expected.coeff(i+1,j),basis[i].deriveTwice()((InputType(expected.coeff(0,j))))));
        }
    }
}

// test 2D nurbs basis second derivative (functions are accessibile and callable)
TEST(nurbs_test, nurbs_basis_second_derivative_2D) {

    std::array<std::vector<double>,2> nodes;
    MdArray<double,MdExtents<Dynamic,Dynamic>>weights(4,5);
    nodes[0].resize(2);
    nodes[1].resize(3);

    using InputType = cexpr::Vector<double, 2>;

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
} // namespace fdapde::testing