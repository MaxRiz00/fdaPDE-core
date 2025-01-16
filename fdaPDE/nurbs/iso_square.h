#ifndef __ISO_SQUARE_H__

#include "iso_mesh.h"

namespace fdapde {

template<int M, int N>
class IsoSquare{

    //template<typename T> using DMatrix = MdArray<T, Dynamic, Dynamic>;

    IsoMesh<M,N>* mesh_; // mesh pointer

    std::size_t ID_; // element ID

    std::array<int, M> left_coords_; // coordinates of the left corner of the element
    std::array<int, M> right_coords_; // coordinates of the right corner of the element

    bool boundary_; // true if the element has at least one vertex on the boundary

    // should implement "parametrization gradient", "metric tensor" e "metric determinant", gometric properties of the element


    public:

    //construct an element a mesh pointer
    IsoSquare(std::array<size::t,M> eMultiIndex, const IsoMesh<M,N>* mesh): mesh_(mesh) {

        boundary_=false;

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
            // check if the element is on the boundary
            if(left_coords_[i] == 0 || right_coords_[i] == mesh_->knots_[i].size() - 1){
                boundary_ = true;
            }
        }

    };


    // affine map from the reference element p in [-1,1] to the parametric domain [left_coords, right_coords]
    std::array<double, M> affine_map(const std::array<double, M> & p) const {
            std::array<double, M> x;
            for(std::size_t i = 0; i < M; ++i){
                x[i] = 0.5*(right_coords_[i] + left_coords_[i] + (right_coords_[i] - left_coords_[i]) * p[i]);
            }
            return x;
        }

    
    // All these methods wants a point p in \Omega, [-1,1]^M, and return the corresponding value in the physical domain
    // parametrization gradient F, puo' magari essere resa piu' efficiente

    Eigen::Matrix<double,N,M> parametrization_gradient(const std::array<double, M>& p) const {
        // compute x, with the affine map
        auto x = affine_map(p);
        Eigen::Matrix<double,N,M> F;
        for(int i=0;i<N;i++){
            for(int j=0;j<M;j++){
                F(i,j) = mesh_->eval_param_derivative(x,i,j);
            }
        }
        return F;
    }

    // metric tensor F^T F

    Eigen::Matrix<double,M,M> metric_tensor(const std::array<double, M>& p) const {
        auto x = affine_map(p);
        auto F = parametrization_gradient(x);
        return F.transpose() * F; // non funziona, there are no operations between matrices, guarda per il trasposto
    }

    // metric determinant sqrt(det(F^T F))

    double metric_determinant(const std::array<double, M>& p) const {
        auto x = affine_map(p);
        auto metric_tensor = metric_tensor(x);
        return sqrt(metric_tensor.determinant()); // non funziona, there are no operations between matrices, guarda per il trasposto
    }

};


}; // namespace fdapde

#endif // __ISO_SQUARE_H__
