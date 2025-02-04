#ifndef __FDAPDE_ISO_CELL_H__
#define __FDAPDE_ISO_CELL_H__

#include "header_check.h"

namespace fdapde {

// Md Hypercube
template<int Order_, int EmbedDim_> class IsoCell{
    friend class IsoMesh<Order_,EmbedDim_>;
    static_assert(Order_ >= 0 && Order_ <= 3);
    public:
    static constexpr int local_dim = Order_;
    static constexpr int embed_dim = EmbedDim_;
    static constexpr int n_nodes = 1 << Order_;
    static constexpr int n_edges = Order_ * (1 << (Order_ - 1));
    static constexpr int n_faces = Order_ * (Order_ - 1) / 2 * (1 << (Order_ - 2));
    static constexpr int n_nodes_per_face = 1 << (Order_ - 1);
    using BoundaryCellType = std::conditional_t<Order_ == 0, IsoCell<0, EmbedDim_>, IsoCell<Order_ - 1, EmbedDim_>>;
    using NodeType = Eigen::Matrix<double, embed_dim, 1>;

    IsoCell() = default;

    IsoCell(const std::array<int, LocalDim>& left_coords, const std::array<int, LocalDim>& right_coords ): left_coords_(left_coords), right_coords_(right_coords){} 


    std::array<double, Order_> affine_map(const std::array<double, Order_> & p) const {
            std::array<double, Order_> x;
            for(std::size_t i = 0; i < Order_; ++i){
                x[i] = 0.5*(right_coords_[i] + left_coords_[i] + (right_coords_[i] - left_coords_[i]) * p[i]);
            }
            return x;
        }

    // Affine map from reference domain [-1, 1]^M to parametric domain [left_coords, right_coords]^M
    // left_coords
    // map from refernce to parameric domain, map_to_parametric, left_coord e right_coord li prende dalla mesh
    Eigen::Matrix<double, EmbedDim_, 1> parametrization(const std::array<double, Order_>& p) const {
        return mesh_->eval_param(affine_map(p));
    }

    Eigen::Matrix<double, EmbedDim_, Order_, Eigen::RowMajor> parametrization_gradient(const std::array<double, Order_>& p) const {
        return mesh_->eval_param_derivative(affine_map(p));
    }

    // Metric tensor F^T * F
    Eigen::Matrix<double, Order_, Order_, Eigen::RowMajor> metric_tensor(const std::array<double, Order_>& p) const {
        auto F = parametrization_gradient(affine_map(p));
        return F.transpose() * F; 
    }

    // metric determinant sqrt(det(F^T * F)), array diventano matrici eigen
    double metric_determinant(const std::array<double, Order_>& p) const {
        return std::sqrt(metric_tensor(affine_map(p)).determinant()); 
    }

    // Compute the neighbors of the current element, diventa neighbors
    Eigen::Matrix<int, 1, 2 * Order_> neighbors() const {
        return mesh_->neighbors().row(id_);
    }

    double parametric_measure() const {
        double measure = 1.0;
        for(std::size_t i = 0; i < Order_; i++){
            measure *= right_coords_[i] - left_coords_[i];
        }
        return measure/(1<<Order_);
    }

    protected:

    std::array<int, Order_> left_coords_ ; // coordinates of the left corner of the element
    std::array<int, Order_> right_coords_ ; // coordinates of the right corner of the element
    // capire se salvare altre quantità

}


};

#endif // __FDAPDE_ISO_CELL_H__