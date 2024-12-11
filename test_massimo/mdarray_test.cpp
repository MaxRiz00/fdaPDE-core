#include <iostream>
#include "linear_algebra/mdarray.h"
#include "fields/nurbs.h"

using namespace fdapde;

int main() {


    /*


    // Declare a 5x4x3 MdArray of doubles
    MdArray<double, MdExtents<Dynamic, Dynamic, Dynamic>> md;
    md.resize(5, 4, 3);

    // Fill the MdArray
    for (int i = 0; i < md.extent(0); ++i) {
        for (int j = 0; j < md.extent(1); ++j) {
            for (int k = 0; k < md.extent(2); ++k) {
                md(i, j, k) = 1.0;
            }
        }
    }

    // Print a value
    std::cout << "Hello, here the value of md(0,0,0): " << md(0, 0, 0) << std::endl;


    auto blk = md.block(full_extent, std::pair{0, 1}, 0); // block of a matrix of the same order

    std::cout << "Hello, here the value of the block: " << blk(0,0,0) << std::endl;

    auto slc = md.template slice<1>(1);

    std::cout << "Hello, here the value of the block: " << slc(0,0) << std::endl;
    */



    return 0;
}