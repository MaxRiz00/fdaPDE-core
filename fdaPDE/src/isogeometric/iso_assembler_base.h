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

#ifndef __FDAPDE_ISO_ASSEMBLER_BASE_H__
#define __FDAPDE_ISO_ASSEMBLER_BASE_H__

#include "header_check.h"

namespace fdapde{

template <typename Derived_> struct IsoMap;

enum class iso_assembler_flags{
    compute_shape_values        = 0x0001,
    compute_shape_grad          = 0x0002,
    compute_shape_hess          = 0x0004,
    compute_shape_div           = 0x0008,
    compute_physical_quad_nodes = 0x0010,
    compute_cell_id             = 0x0020

};

namespace internals {

// informations sent from the assembly loop to the integrated forms
template <int LocalDim> struct iso_assembler_packet {
    static constexpr int local_dim = LocalDim;
    iso_assembler_packet(int n_trial_components, int n_test_components) :
    trial_value(n_trial_components),
    test_value (n_test_components ),
    trial_grad (n_trial_components),
    test_grad  (n_test_components ),
    trial_hess (n_trial_components),
    test_hess  (n_test_components ) { }
    iso_assembler_packet(int n_components) : iso_assembler_packet(n_components, n_components) { }
    iso_assembler_packet() : iso_assembler_packet(1, 1) { }
    iso_assembler_packet(iso_assembler_packet&&) noexcept = default;
    iso_assembler_packet(const iso_assembler_packet&) noexcept = default;

    // geometric informations
    int quad_node_id;       // active physical quadrature node index
    double cell_measure;    // active cell measure
    double cell_id;         // active cell identifier
    double cell_diameter;   // active cell diameter

    // functional informations (Dynamic stands for number of components)
    MdArray<double, MdExtents<Dynamic>> trial_value, test_value;            // \psi_i(q_k), \psi_j(q_k)
    MdArray<double, MdExtents<Dynamic, local_dim>> trial_grad, test_grad;   // \nabla{\psi_i}(q_k), \nabla{\psi_j}(q_k)
    MdArray<double, MdExtents<Dynamic, local_dim, local_dim>> trial_hess, test_hess;
    double trial_div = 0, test_div = 0;
};


// base class for vector finite element assembly loops
template<typename IsoMesh_, typename Form_, int Options_, typename... Quadrature_>
struct iso_assembler_base{
    fdapde_static_assert(sizeof...(Quadrature_) < 2, YOU_CAN_SUPPLY_AT_MOST_ONE_QUADRATURE_RULE_TO_A_ISO_ASSEMBLY_LOOP);
    // detect test space (since a test function is always present in a weak form)
    using TestSpace = test_space_t<Form_>;
    using Form =
      std::decay_t<decltype(xpr_wrap<IsoMap, decltype([]<typename Xpr>() {
	    return !(
	        std::is_invocable_v<Xpr, iso_assembler_packet<Xpr::StaticInputSize>>);
	  })>(std::declval<Form_>()))>; // vector case ???
    using IsoMesh = typename std::decay<IsoMesh_>;
    static constexpr int local_dim = IsoMesh::local_dim;
    static constexpr int embed_dim = IsoMesh::embed_dim;
    static constexpr int Options = Options_;
    using FunctionSpace = TestSpace;
    using DofHandlerType = DofHandler<local_dim, embed_dim, iso_tag>;
    using Quadrature = decltype([]() {
        if constexpr (sizeof...(Quadrature_) == 0) {
            return void();   // quadrature selcted at run-time provided the actual order of spline basis
        } else {
            return std::get<0>(std::tuple<Quadrature_...>());   // user-defined quadrature
        }
    }());
    using geo_iterator = typename IsoMesh::cell_iterator;
    using dof_iterator = typename DofHandlerType::cell_iterator;
    using discretization_category = typename TestSpace::discretization_category;
    fdapde_static_assert(
        std::is_same_v<discretization_category FDAPDE_COMMA iso_tag>, THIS_CLASS_IS_FOR_ISO_DISCRETIZATION_ONLY);

    iso_assembler_base() = default;

    iso_assembler_base( const Form_& form, const geo_iterator& begin, const geo_iterator& end, const Quadrature_&... quadrature )
        requires(sizeof...(quadrature)<=1):
        form_(xpr_wrap<IsoMap, decltype([]<typename Xpr>() {
                    return !(std::is_invocable_v<Xpr, iso_assembler_packet<Xpr::StaticInputSize>>);
                })>(form)),
        dof_handler_(std::addressof(internals::test_space(form_).dof_handler())),
        test_space_ (std::addressof(internals::test_space(form_))),
        begin_(begin),
        end_(end) { 
            fdapde_assert(dof_handler_->n_dofs() > 0);
            // copy quadrature rule
            Eigen::Matrix<double, Dynamic, Dynamic> quad_nodes__;
            if constexpr (sizeof...(quadrature) == 1) {
                auto quad_rule = std::get<0>(std::make_tuple(quadrature...));
                fdapde_assert(local_dim == quad_rule.local_dim);
                quad_nodes__.resize(quad_rule.order, quad_rule.local_dim);
                quad_weights_.resize(quad_rule.order, 1);
                for (int i = 0; i < quad_rule.order; ++i) {
                    quad_weights_(i, 0) = quad_rule.weights[i];
                    for (int j = 0; j < local_dim; ++j) { quad_nodes__(i, j) = quad_rule.nodes(i, j); }
            }
            } else {
                //internals::get_sp_quadrature(test_space_->order(), quad_nodes__, quad_weights_);
            }
            // build grid of quadrature nodes on reference domain
            n_quadrature_nodes_ = quad_nodes__.rows();
            int n_cells = std::distance(begin_, end_), n_src_points = quad_nodes__.rows();

            quad_nodes_.resize(n_cells * n_src_points, local_dim); // global quad nodes
            int i  = 0;
            for(auto it = begin_; it != end_; ++it){
                for(int q_k = 0; q_k < n_quadrature_nodes_; ++q_k){
                    quad_nodes_.row(i) = it->affine_map(quad_nodes__.row(q_k).transpose());
                    i++;
                }
                
            }
            return;

        }

        const TestSpace& test_space() const { return *test_space_; }

        protected:
        
        // evaluation of \psi_i(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
        template<typename BasisType__, typename IteratorType, typename DstMdArray>
        void eval_shape_values(
            BasisType__&& basis, const std::vector<int>& active_dofs, IteratorType cell, DstMdArray& dst) const {

            using BasisType = std::decay_t<BasisType__>;
            int n_basis =active_dofs.size(); // attenzione 1d

            for(int i=0; i < n_basis; ++i){
                // evaluation of \psi_i at q_j, j = 1, ..., n_quadrature_nodes
                for(int j=0; j < n_quadrature_nodes_; ++j){
                    dst(i, j) = basis[active_dofs[i]](quad_nodes_.row(cell->id() * n_quadrature_nodes_ + j).transpose());
                }
            }
            return;
        }
        

        // evaluation of 1-st order derivative of basis function
        template <typename BasisType__, typename IteratorType, typename DstMdArray>
            requires(requires(BasisType__ basis, int i, int k) { basis[i].derive(k); })
            void eval_shape_grads(
            BasisType__&& basis, const std::vector<int>& active_dofs, IteratorType cell, DstMdArray& dst) const {
            using BasisType = std::decay_t<BasisType__>;
            using DerivativeType = decltype(std::declval<BasisType>()[std::declval<int>()].derive(std::declval<int>()));
            int n_basis = active_dofs.size();
            for (int i = 0; i < n_basis; ++i) {
                for(int k = 0; k < local_dim; ++k){
                    DerivativeType der = basis[active_dofs[i]].derive(k);
                    for (int j = 0; j < n_quadrature_nodes_; ++j) {        
                        //evaluation of \nabla{\psi_i}(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
                        dst(i, j, k) = der(quad_nodes_.row(cell->id() * n_quadrature_nodes_ + j).transpose());
                    }
                }
            }
            return;
        }

        //evaluation of hessian matrix of \psi_i(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
        template <typename BasisType__, typename IteratorType, typename DstMdArray>
        void eval_shape_hess(
            BasisType__&& basis, const std::vector<int>& active_dofs, IteratorType cell, DstMdArray& dst) const {
            using BasisType = std::decay_t<BasisType__>;
            int n_basis = active_dofs.size();
            for (int i = 0; i < n_basis; ++i) {
                for(int k = 0; k < local_dim; ++k){
                    for(int l = 0; l < local_dim; ++l){
                        auto hess = basis[active_dofs[i]].deriveTwice(k,l)(quad_nodes_.row(cell->id() * n_quadrature_nodes_).transpose());
                        for (int j = 0; j < n_quadrature_nodes_; ++j) {
                            //evaluation of \nabla{\psi_i}(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
                            dst(i, j, k, l) = hess((quad_nodes_.row(cell->id() * n_quadrature_nodes_).transpose()));
                        }
                    }
                }
            }
            return;
        }
        











    protected:
    Form form_;
    const DofHandlerType* dof_handler_;
    const TestSpace* test_space_;
    geo_iterator begin_, end_;
    // quadrature
    Eigen::Matrix<double, Dynamic, Dynamic> quad_nodes_, quad_weights_;
    int n_quadrature_nodes_;





};




}






}





#endif