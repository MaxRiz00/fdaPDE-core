#include <iostream>
#include "linear_algebra/mdarray.h"
#include "fields/nurbs.h"
#include "nurbs/nurbs_basis.h"

using namespace fdapde;

int main() {

    // create a Nurbs object
    std::array<std::vector<double>, 2> knots;
    knots[0] = {0.0, 0.25, 0.5, 0.75, 1.0};
    knots[1] = {0.0, 0.25, 0.5, 0.75};

    MdArray<double, full_dynamic_extent_t<2>> weights(7,6);
    // initialize the weights
    for (int i = 0; i < weights.extent(0); ++i) {
        for (int j = 0; j < weights.extent(1); ++j) {
            weights(i, j) = i + j;
        }
    }

    std::array<int, 2> index = {2,2};
    int order = 3;


    Nurbs<2> nurbs(knots, weights, index, order);

    // evaluate the nurbs at a point
    cexpr::Vector<double, 2> p{0.125,0.125};
    std::cout << "Hello, here the value of the nurbs at point p: " << nurbs(p) << std::endl;


    // create a NurbsBasis object
    NurbsBasis<2> nurbs_basis(knots, weights, order);

/*

    // Declare a 5x4x3 MdArray of doubles
    MdArray<double, MdExtents<Dynamic, Dynamic>> md;

    // Resize the MdArray
    md.resize(5, 4);

    // Fill the MdArray
    for (int i = 0; i < md.extent(0); ++i) {
        for (int j = 0; j < md.extent(1); ++j) {
            md(i, j) = i + j;
        }
    }


    auto sub_block = md.template slice<0>(1);

*/

    return 0;
}