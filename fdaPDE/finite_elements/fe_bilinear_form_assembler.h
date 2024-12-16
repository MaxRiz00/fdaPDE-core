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

#ifndef __FE_BILINEAR_FORM_ASSEMBLER_H__
#define __FE_BILINEAR_FORM_ASSEMBLER_H__

#include <unordered_map>

#include "fe_assembler_base.h"

namespace fdapde {
namespace internals {
  
// galerkin and petrov-galerkin vector finite element assembly loop
template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
class fe_bilinear_form_assembly_loop :
    public fe_assembler_base<Triangulation_, Form_, Options_, Quadrature_...>,
    public assembly_xpr_base<fe_bilinear_form_assembly_loop<Triangulation_, Form_, Options_, Quadrature_...>> {
   public:
    // detect trial and test spaces from bilinear form
    using TrialSpace = trial_space_t<Form_>;
    using TestSpace  = test_space_t <Form_>;
    static_assert(TrialSpace::local_dim == TestSpace::local_dim && TrialSpace::embed_dim == TestSpace::embed_dim);
    static constexpr bool is_galerkin = std::is_same_v<TrialSpace, TestSpace>;
    static constexpr bool is_petrov_galerkin = !is_galerkin;
    using Base = fe_assembler_base<Triangulation_, Form_, Options_, Quadrature_...>;
    using Form = typename Base::Form;
    using DofHandlerType = typename Base::DofHandlerType;
    static constexpr int local_dim = Base::local_dim;
    static constexpr int embed_dim = Base::embed_dim;
    using Base::form_;
    // as trial and test spaces could be different, we here need to redefine some properties of Base
    // trial space properties
    using TrialFeType = typename TrialSpace::FeType;
    using trial_fe_traits = std::conditional_t<
      Options_ == CellMajor, fe_cell_assembler_traits<TrialSpace, Quadrature_...>,
      fe_face_assembler_traits<TrialSpace, Quadrature_...>>;
    using trial_dof_descriptor = typename trial_fe_traits::dof_descriptor;
    using TrialBasisType = typename trial_dof_descriptor::BasisType;
    static constexpr int n_trial_basis = TrialBasisType::n_basis;
    static constexpr int n_trial_components = TrialSpace::n_components;
    // test space properties
    using TestFeType = typename TestSpace::FeType;
    using test_fe_traits = std::conditional_t<
      Options_ == CellMajor, fe_cell_assembler_traits<TestSpace, Quadrature_...>,
      fe_face_assembler_traits<TestSpace, Quadrature_...>>;
    using test_dof_descriptor = typename test_fe_traits::dof_descriptor;
    using TestBasisType = typename test_dof_descriptor::BasisType;
    static constexpr int n_test_basis = TestBasisType::n_basis;
    static constexpr int n_test_components = TestSpace::n_components;

    using Quadrature = std::conditional_t<
      (is_galerkin || sizeof...(Quadrature_) > 0), typename Base::Quadrature,
      higher_degree_fe_quadrature_t<
        typename TrialFeType::template cell_quadrature_t<local_dim>,
        typename TestFeType ::template cell_quadrature_t<local_dim>>>;
    static constexpr int n_quadrature_nodes = Quadrature::order;
   private:
    // selected Quadrature could be different than Base::Quadrature, evaluate trial and (re-evaluate) test functions
    static constexpr auto test_shape_values_  = Base::template eval_shape_values<Quadrature, test_fe_traits >();
    static constexpr auto trial_shape_values_ = Base::template eval_shape_values<Quadrature, trial_fe_traits>();
    static constexpr auto test_shape_grads_   = Base::template eval_shape_grads <Quadrature, test_fe_traits >();
    static constexpr auto trial_shape_grads_  = Base::template eval_shape_grads <Quadrature, trial_fe_traits>();
    // private data members
    const DofHandlerType* trial_dof_handler_;
    Quadrature quadrature_ {};
    constexpr const DofHandlerType* test_dof_handler() const { return Base::dof_handler_; }
    constexpr const DofHandlerType* trial_dof_handler() const {
        return is_galerkin ? Base::dof_handler_ : trial_dof_handler_;
    }
    const TrialSpace* trial_space_;
   public:
    fe_bilinear_form_assembly_loop() = default;
    fe_bilinear_form_assembly_loop(
      const Form_& form, typename Base::fe_traits::geo_iterator begin, typename Base::fe_traits::geo_iterator end,
      const Quadrature_&... quadrature)
        requires(sizeof...(quadrature) <= 1)
        : Base(form, begin, end, quadrature...), trial_space_(&internals::trial_space(form_)) {
        if constexpr (is_petrov_galerkin) { trial_dof_handler_ = &internals::trial_space(form_).dof_handler(); }
        fdapde_assert(test_dof_handler()->n_dofs() != 0 && trial_dof_handler()->n_dofs() != 0);
    }

    SpMatrix<double> assemble() const {
        SpMatrix<double> assembled_mat(test_dof_handler()->n_dofs(), trial_dof_handler()->n_dofs());
        std::vector<Eigen::Triplet<double>> triplet_list;
	assemble(triplet_list);
	// linearity of the integral is implicitly used here, as duplicated triplets are summed up (see Eigen docs)
        assembled_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
        assembled_mat.makeCompressed();
        return assembled_mat;
    }
    void assemble(std::vector<Eigen::Triplet<double>>& triplet_list) const {
        using iterator = typename Base::fe_traits::dof_iterator;
        iterator begin(Base::begin_.index(), test_dof_handler(), Base::begin_.marker());
        iterator end  (Base::end_.index(),   test_dof_handler(), Base::end_.marker()  );
        // prepare assembly loop
        DVector<int> test_active_dofs, trial_active_dofs;
        MdArray<double, MdExtents<n_test_basis,  n_quadrature_nodes, local_dim, n_test_components >> test_grads;
        MdArray<double, MdExtents<n_trial_basis, n_quadrature_nodes, local_dim, n_trial_components>> trial_grads;
	cexpr::Matrix<double, n_test_basis,  n_quadrature_nodes> test_divs;
	cexpr::Matrix<double, n_trial_basis, n_quadrature_nodes> trial_divs;
	
        std::unordered_map<const void*, DMatrix<double>> fe_map_buff;
        if constexpr (Form::XprBits & int(fe_assembler_flags::compute_physical_quad_nodes)) {
            Base::distribute_quadrature_nodes(
              fe_map_buff, begin, end);   // distribute quadrature nodes on physical mesh (if required)
        }
        // start assembly loop
        internals::fe_assembler_packet<local_dim> fe_packet(n_trial_components, n_test_components);
	int local_cell_id = 0;
        for (iterator it = begin; it != end; ++it) {
            // update fe_packet content based on form requests
            fe_packet.cell_measure = it->measure();
            if constexpr (Form::XprBits & int(fe_assembler_flags::compute_cell_diameter)) {
                fe_packet.cell_diameter = it->diameter();
            }
            if constexpr (Form::XprBits & int(fe_assembler_flags::compute_shape_grad)) {
                Base::eval_shape_grads_on_cell(it, test_shape_grads_, test_grads);
                if constexpr (is_petrov_galerkin) Base::eval_shape_grads_on_cell(it, trial_shape_grads_, trial_grads);
            }
            if constexpr (Form::XprBits & int(fe_assembler_flags::compute_shape_div)) {
                // divergence is defined only for vector elements, skeep computation for scalar element case
                if constexpr (n_test_components != 1) Base::eval_shape_div_on_cell(it, test_shape_grads_, test_divs);
                if constexpr (is_petrov_galerkin && n_trial_components != 1)
                    Base::eval_shape_div_on_cell(it, trial_shape_grads_, trial_divs);
            }

            // perform integration of weak form for (i, j)-th basis pair
            test_active_dofs = it->dofs();
            if constexpr (is_petrov_galerkin) { trial_active_dofs = trial_dof_handler()->active_dofs(it->id()); }
            for (int i = 0; i < n_trial_basis; ++i) {      // trial function loop
                for (int j = 0; j < n_test_basis; ++j) {   // test function loop
                    double value = 0;
                    for (int q_k = 0; q_k < n_quadrature_nodes; ++q_k) {
                        if constexpr (Form::XprBits & int(fe_assembler_flags::compute_shape_values)) {
                            fe_packet.trial_value.assign_inplace_from(trial_shape_values_.template slice<0, 1>(i, q_k));
                            fe_packet.test_value .assign_inplace_from(test_shape_values_ .template slice<0, 1>(j, q_k));
                        }
                        if constexpr (Form::XprBits & int(fe_assembler_flags::compute_shape_grad)) {
                            fe_packet.trial_grad.assign_inplace_from(is_galerkin ?
                                test_grads.template slice<0, 1>(i, q_k) : trial_grads.template slice<0, 1>(i, q_k));
                            fe_packet.test_grad .assign_inplace_from(test_grads.template slice<0, 1>(j, q_k));

                        }
                        if constexpr (Form::XprBits & int(fe_assembler_flags::compute_shape_div)) {
                            if constexpr (n_trial_components != 1) {
                                fe_packet.trial_div =
                                  (is_galerkin && n_test_components != 1) ? test_divs(i, q_k) : trial_divs(i, q_k);
                            }
                            if constexpr (n_test_components != 1) fe_packet.test_div = test_divs(j, q_k);
                        }
                        if constexpr (Form::XprBits & int(fe_assembler_flags::compute_physical_quad_nodes)) {
                            fe_packet.quad_node_id = local_cell_id * n_quadrature_nodes + q_k;
                        }
                        value += Quadrature::weights[q_k] * form_(fe_packet);
                    }
                    triplet_list.emplace_back(
                      test_active_dofs[j], is_galerkin ? test_active_dofs[i] : trial_active_dofs[i],
                      value * fe_packet.cell_measure);
                }
            }
            local_cell_id++;
        }
        return;
    }
    constexpr int n_dofs() const { return trial_dof_handler()->n_dofs(); }
    constexpr int rows() const { return test_dof_handler()->n_dofs(); }
    constexpr int cols() const { return trial_dof_handler()->n_dofs(); }
    constexpr const TrialSpace& trial_space() const { return *trial_space_; }
};


// optimized computation of discretized laplace operator (\int_D (\grad{\psi_i} * \grad{\psi_j})) for scalar elements
template <typename DofHandler, typename FeType> class scalar_fe_grad_grad_assembly_loop {
    static constexpr int local_dim = DofHandler::local_dim;
    static constexpr int embed_dim = DofHandler::embed_dim;
    using Quadrature = typename FeType::template cell_quadrature_t<local_dim>;
    using cell_dof_descriptor = FeType::template cell_dof_descriptor<local_dim>;
    using BasisType = typename cell_dof_descriptor::BasisType;
    static constexpr int n_quadrature_nodes = Quadrature::order;
    static constexpr int n_basis = BasisType::n_basis;
    static constexpr int n_components = FeType::n_components;
    fdapde_static_assert(n_components == 1, THIS_CLASS_IS_FOR_SCALAR_FINITE_ELEMENTS_ONLY);
    // compile-time evaluation of \nabla{\psi_i}(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
    static constexpr std::array<cexpr::Matrix<double, local_dim, n_quadrature_nodes>, n_basis> shape_grad_ {[]() {
        std::array<cexpr::Matrix<double, local_dim, n_quadrature_nodes>, n_basis> shape_grad_ {};
        BasisType basis {cell_dof_descriptor().dofs_phys_coords()};
        for (int i = 0; i < n_basis; ++i) {
            // evaluation of \nabla{\psi_i} at q_j, j = 1, ..., n_quadrature_nodes
            std::array<double, n_quadrature_nodes * local_dim> grad_eval_ {};
            for (int k = 0; k < n_quadrature_nodes; ++k) {
                auto grad = basis[i].gradient()(Quadrature::nodes.row(k).transpose());
                for (int j = 0; j < local_dim; ++j) { grad_eval_[j * n_quadrature_nodes + k] = grad[j]; }
            }
            shape_grad_[i] = cexpr::Matrix<double, local_dim, n_quadrature_nodes>(grad_eval_);
        }
        return shape_grad_;
    }()};
    DofHandler* dof_handler_;
   public:
    scalar_fe_grad_grad_assembly_loop() = default;
    scalar_fe_grad_grad_assembly_loop(DofHandler& dof_handler) : dof_handler_(&dof_handler) { }
    SpMatrix<double> assemble() {
        if (!dof_handler_) dof_handler_->enumerate(FeType {});
        int n_dofs = dof_handler_->n_dofs();
        SpMatrix<double> assembled_mat(n_dofs, n_dofs);
        // prepare assembly loop
        std::vector<Eigen::Triplet<double>> triplet_list;
        DVector<int> active_dofs;
        std::array<cexpr::Matrix<double, local_dim, n_quadrature_nodes>, n_basis> shape_grad;
        for (typename DofHandler::cell_iterator it = dof_handler_->cells_begin(); it != dof_handler_->cells_end();
             ++it) {
            active_dofs = it->dofs();
            // map reference cell shape functions' gradient on physical cell
            for (int i = 0; i < n_basis; ++i) {
                for (int j = 0; j < n_quadrature_nodes; ++j) {
                    shape_grad[i].col(j) =
                      cexpr::Map<const double, local_dim, embed_dim>(it->invJ().data()).transpose() *
                      shape_grad_[i].col(j);
                }
            }
            for (int i = 0; i < BasisType::n_basis; ++i) {
                for (int j = 0; j < i + 1; ++j) {
                    std::pair<const int&, const int&> minmax(std::minmax(active_dofs[i], active_dofs[j]));
                    double value = 0;
                    for (int k = 0; k < n_quadrature_nodes; ++k) {
                        value += Quadrature::weights[k] * (shape_grad[i].col(k)).dot(shape_grad[j].col(k));
                    }
                    triplet_list.emplace_back(minmax.first, minmax.second, value * it->measure());
                }
            }
        }
        // linearity of the integral is implicitly used here, as duplicated triplets are summed up (see Eigen docs)
        assembled_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
        assembled_mat.makeCompressed();
	return assembled_mat.selfadjointView<Eigen::Upper>();
    }
    constexpr int n_dofs() const { return dof_handler_->n_dofs(); }
    constexpr int rows() const { return n_dofs(); }
    constexpr int cols() const { return n_dofs(); }  
};
  
}   // namespace internals
}   // namespace fdapde

#endif   // __FE_BILINEAR_FORM_ASSEMBLER_H__
