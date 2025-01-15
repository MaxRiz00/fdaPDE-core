#ifndef __ISO_MESH_H__

#include "nurbs_basis.h"

namespace fdapde {

template<int M, int N>
class IsoMesh {

    template<typename T> using DMatrix = MdArray<T, Dynamic, Dynamic>;

    std::array<std::vector<double>,M> knots_; // knots in each direction
    MdArray<double,full_dynamic_extent_t<M>> weights_; // weights in each direction
    MdArray<double,full_dynamic_extent_t<M+1>> control_points_; // control points in each direction

    
    NurbsBasis<M> basis_; // basis of the mesh


    //std::vector<ElementIGA<M,N>> elements_cache_; // elements of the mesh

    //DMatrix<double> nodes_ {}; // nodes of the mesh, size: r by M, where r is the product of the numbers of unique knots along each dimension
    //std::vector<std::size_t> boundary_ {}; // : r by 1 contains a boolean value for each node, true if the node is on the boundary of Ω
    //DMatrix<std::size_t> elements_ {}; // elements of the mesh
    //DMatrix<std::size_t> boundary_dof_ {}; // dofs of boundary of the mesh

    // implementation of the private methods for the mesh parametrization

    inline double eval_param(const std::array<double, M>& u, std::size_t derivative_index_) const {
        double x = 0.0;
        for (const auto& nurb : this->nurbs_basis_) {
            x += nurb(x) * this->control_points_(nurb.index());
        }
        return x;
    }

    inline double eval_param_derivative(const std::array<double, M>& u) const {
        double x = 0.0;
        for (const auto& nurb : this->nurbs_basis_) {
            x += nurb.derive(derivative_index_)(x) * this->control_points_(nurb.index());
        }
        return x;
    }


    // HERE

    //std::array<Parametrization_type,N> parametrizations_;
    //std::array<Derivative_type,N> gradients_;

    public:

    template<int M, int N>
    IsoMesh<M,N>::IsoMesh(const std::array<std::vector<double>,M> & knots, const MdArray<double,full_dynamic_extent_t<M>> & weights, 
        const MdArray<double,full_dynamic_extent_t<M+1>> & control_points): knots_(knots), weights_(weights), control_points_(control_points){

            nurbs_basis_ = NurbsBasis(knots, weights, order);

        };

    // some getters

    const NurbsBasis<M>& basis() const { return basis_; }
    const MdArray<double, full_dynamic_extent_t<M+1>>& control_points() const { return control_points_; }
    const std::array<std::vector<double>,M>& knots() const { return knots_; }

};






}; // namespace fdapde

#endif // __ISO_MESH_H__