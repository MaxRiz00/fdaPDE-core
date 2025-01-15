#ifndef __ISO_SQUARE_H__

#include "iso_mesh.h"

namespace fdapde {

template<int M, int N>
class IsoSquare{

    IsoMesh<M,N>* mesh_; // mesh pointer

    std::size_t ID_; // element ID

    std::array<int, M> left_coords_; // coordinates of the left corner of the element
    std::array<int, M> right_coords_; // coordinates of the right corner of the element

    bool boundary_; // true if the element has at least one vertex on the boundary

    // should implement "parametrization gradient", "metric tensor" e "metric determinant", gometric properties of the element


    public:

    //construct an element a mesh pointer
    IsoSquare(std::array<size::t,M> eMultiIndex, const IsoMesh<M,N>* mesh): mesh_(mesh) {

        // Compute the ID from the multi-index using the lexicographic order
        ID_ = 0;
        for(std::size_t i = 0; i < M; ++i){
            if(i==0) 
                ID_ = eMultiIndex[i];
            else
                ID_ += (mesh_->knots_[i-1].size() - 1) + eMultiIndex[i];
        }

        // compute the coordinates of the left and right corners of the element, basing on the multi-index
        for(std::size_t i = 0; i < M; ++i){
            left_coords_[i] = eMultiIndex[i];
            right_coords_[i] = eMultiIndex[i] + 1;
        }

    };


};


}; // namespace fdapde

#endif // __ISO_SQUARE_H__
