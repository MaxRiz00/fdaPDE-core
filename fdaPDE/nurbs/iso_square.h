#ifndef __ISO_SQUARE_H__
#define __ISO_SQUARE_H__

#include <array>
#include <cmath>
#include <Eigen/Dense>



namespace fdapde {
        
// Forward declaration of IsoMesh
template<int M, int N> class IsoMesh;

template<int M, int N>
class IsoSquare{

    friend class IsoMesh<M,N>;

    //template<typename T> using DMatrix = MdArray<T, Dynamic, Dynamic>;

    const IsoMesh<M,N>* mesh_; // mesh pointer

    std::size_t ID_; // element ID

    std::array<std::size_t,M> eMultiIndex_; // element multi-index, da capire se tenerlo

    std::array<int, M> left_coords_; // coordinates of the left corner of the element
    std::array<int, M> right_coords_; // coordinates of the right corner of the element

    bool boundary_; // true if the element has at least one vertex on the boundary

    public:

    /// Default constructor
    IsoSquare() = default;

    //construct an element a mesh pointer
    IsoSquare(std::array<std::size_t,M> eMultiIndex,const IsoMesh<M,N>* mesh): mesh_(mesh), eMultiIndex_(eMultiIndex){

        /*
        std::cout << "eMultiIndex della cella in cui siamo: ";
        for(int i=0;i<M;i++){
            std::cout << eMultiIndex[i] << " ";
        }
        std::cout << std::endl;
        */

        // Compute the ID from the multi-index using the lexicographic order
        ID_ = mesh_->compute_el_ID(eMultiIndex);

        //std::cout << "ID della cella in cui siamo: " << ID_ << std::endl;

        left_coords_ = mesh_->compute_param_el_vertices(ID_)[0];
        right_coords_ = mesh_->compute_param_el_vertices(ID_)[1];


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
        return sqrt(metric_tensor(x).determinant()); // non funziona, there are no operations between matrices, guarda per il trasposto
    }

    // da aggiungere misura dell'elemento fisico

    // some getters
    std::size_t ID() const { return ID_; }
    const std::array<int,M> & left_coords() const { return left_coords_; }
    const std::array<int,M> & right_coords() const { return right_coords_; }
    bool is_boundary() const { return boundary_; }

    // get the neighbors of the element
    std::vector<std::size_t> get_neighbors_ID() const {
        return mesh_->get_neighbors_ID(ID_);
    }


};


}; // namespace fdapde

#endif // __ISO_SQUARE_H__
