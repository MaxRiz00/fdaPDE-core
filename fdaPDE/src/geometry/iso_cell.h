#ifndef __FDAPDE_ISO_CELL_H__
#define __FDAPDE_ISO_CELL_H__

#include "header_check.h"


namespace fdapde {

// Md Hypercube
template<int LocalDim_, int EmbedDim_> class IsoCell{
    static_assert(LocalDim_ >= 0 && LocalDim_ <= 3);
    public:
    static constexpr int local_dim = LocalDim_;
    static constexpr int embed_dim = EmbedDim_;
    static constexpr int n_nodes = 1 << LocalDim_;
    static constexpr int n_edges = LocalDim_ * (1 << (LocalDim_ - 1));
    static constexpr int n_faces = LocalDim_ * (LocalDim_ - 1) / 2 * (1 << (LocalDim_ - 2));
    static constexpr int n_nodes_per_face = 1 << (LocalDim_ - 1);
    using BoundaryCellType = std::conditional_t<LocalDim_ == 0, IsoCell<0, EmbedDim_>, IsoCell<LocalDim_ - 1, EmbedDim_>>;
    using NodeType = Eigen::Matrix<double, embed_dim, 1>;

    IsoCell() = default;

    IsoCell(const Eigen::Matrix<double, local_dim, 1> left_coords, const Eigen::Matrix<double, local_dim, 1> right_coords ): 
        left_coords_(left_coords), right_coords_(right_coords) { } 

    //commenta le funzioni
    Eigen::Matrix<double, local_dim,1> affine_map(const Eigen::Matrix<double, local_dim,1> & p) const {
        Eigen::Matrix<double, local_dim,1> x;
            for(std::size_t i = 0; i < LocalDim_; ++i){
                x(i) = 0.5*(right_coords_(i) + left_coords_(i) + (right_coords_(i) - left_coords_(i)) * p(i));
                //std::cout<<"p("<<i<<") = "<<p(i)<<std::endl;
                //std::cout<<"x("<<i<<") = "<<x(i)<<std::endl;
            }
            return x;
        }

    double parametric_measure() const {
        double measure = 1.0;
        for(std::size_t i = 0; i < LocalDim_; i++){
            measure *= right_coords_(i) - left_coords_(i);
        }
        return measure/(1<<LocalDim_);
    }

    // getters
    Eigen::Matrix<double, local_dim, 1> left_coords() const { return left_coords_; }
    Eigen::Matrix<double, local_dim, 1> right_coords() const { return right_coords_; }

    protected:

    Eigen::Matrix<double, local_dim, 1> left_coords_ {} ; // coordinates of the left corner of the element
    Eigen::Matrix<double, local_dim, 1> right_coords_ {} ; // coordinates of the right corner of the element
    // capire se salvare altre quantità

};


};

#endif // __FDAPDE_ISO_CELL_H__