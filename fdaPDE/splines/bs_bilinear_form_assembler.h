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

#ifndef __BS_BILINEAR_FORM_ASSEMBLER_H__
#define __BS_BILINEAR_FORM_ASSEMBLER_H__

#include <unordered_map>

#include "bs_assembler_base.h"

namespace fdapde {
namespace internals {
  
// galerkin and petrov-galerkin spline assembly loop
template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
class bs_bilinear_form_assembly_loop :
    public bs_assembler_base<Triangulation_, Form_, Options_, Quadrature_...>,
    public assembly_xpr_base<bs_bilinear_form_assembly_loop<Triangulation_, Form_, Options_, Quadrature_...>> {
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
    // private data members
    const DofHandlerType* trial_dof_handler_;
    Quadrature quadrature_ {};
    constexpr const DofHandlerType* test_dof_handler() const { return Base::dof_handler_; }
    constexpr const DofHandlerType* trial_dof_handler() const {
        return is_galerkin ? Base::dof_handler_ : trial_dof_handler_;
    }
    const TrialSpace* trial_space_;
    MdArray<double, MdExtents<Dynamic, Dynamic>> shape_values_, shape_dx_, shape_ddx_;
   public:
    bs_bilinear_form_assembly_loop() = default;
    bs_bilinear_form_assembly_loop(
      const Form_& form, typename Base::fe_traits::geo_iterator begin, typename Base::fe_traits::geo_iterator end,
      const Quadrature_&... quadrature)
        requires(sizeof...(quadrature) <= 1)
        : Base(form, begin, end, quadrature...), trial_space_(std::adressof(internals::trial_space(form_))) {
        if constexpr (is_petrov_galerkin) {
            trial_dof_handler_ = std::addressof(internals::trial_space(form_).dof_handler());
        }
        fdapde_assert(test_dof_handler()->n_dofs() != 0 && trial_dof_handler()->n_dofs() != 0);
        if constexpr (sizeof...(Quadrature_) == 0) {
            // default to higher-order quadrature
            internals::get_bs_quadrature(
              test_space_->order() > trial_space_->order() ? test_space_->order() : trial_space_->order(),
              Base::quad_nodes_, Base::quad_weights_);
        }
	// evaluate basis system on reference interval at quadrature nodes according to form requests
        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_values) {
            test_shape_values_  = Base::eval_shape_values(test_space_ ->basis());
            trial_shape_values_ = Base::eval_shape_values(trial_space_->basis());
        }
        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_dx) {
            test_shape_dx_  = Base::eval_shape_dx(test_space_ ->basis());
            trial_shape_dx_ = Base::eval_shape_dx(trial_space_->basis());
        }
        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_ddx) {
            test_shape_ddx_  = Base::eval_shape_ddx(test_space_ ->basis());
            trial_shape_ddx_ = Base::eval_shape_ddx(trial_space_->basis());
        }
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
        using iterator = typename Base::dof_iterator;
        iterator begin(Base::begin_.index(), test_dof_handler(), Base::begin_.marker());
        iterator begin(Base::end_.index(), test_dof_handler(), Base::end_.marker());
	// prepare assembly loop
        Eigen::Matrix<int, Dynamic, 1> test_active_dofs, trial_active_dofs;
        int n_trial_basis = trial_space_->n_shape_functions();
        int n_test_basis = test_space_->n_shape_functions();
	int n_quadrature_nodes = Base::quad_nodes_.rows();
	
        std::unordered_map<const void*, Eigen::Matrix<double, Dynamic, Dynamic>> bs_map_buff;
        if constexpr (Form::XprBits & bs_assembler_flags::compute_physical_quad_nodes) {
            Base::distribute_quadrature_nodes(
              bs_map_buff, begin, end);   // distribute quadrature nodes on physical mesh (if required)
        }
	
	// start assembly loop
        internals::bs_assembler_packet<local_dim> bs_packet();
        int local_cell_id = 0;
        for (iterator it = begin; it != end; ++it) {
            bs_packet.cell_measure = it->measure();
            // perform integration of weak form for (i, j)-th basis pair
            test_active_dofs = it->dofs();
            if constexpr (is_petrov_galerkin) { trial_active_dofs = trial_dof_handler()->active_dofs(it->id()); }
            for (int i = 0; i < n_trial_basis; ++i) {      // trial function loop
                for (int j = 0; j < n_test_basis; ++j) {   // test function loop
                    double value = 0;
                    for (int q_k = 0; q_k < n_quadrature_nodes; ++q_k) {
                        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_values) {
                            bs_packet.trial_value = Base::trial_shape_values_(i, q_k);
                            bs_packet.test_value  = Base::test_shape_values_ (j, q_k);
                        }
                        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_dx) {
                            bs_packet.trial_dx = Base::trial_shape_dx_(i, q_k);
                            bs_packet.test_dx  = Base::test_shape_dx_ (j, q_k);
                        }
                        if constexpr (Form::XprBits & bs_assembler_flags::compute_shape_ddx) {
                            bs_packet.trial_ddx = Base::trial_shape_ddx_(i, q_k);
                            bs_packet.test_ddx  = Base::test_shape_ddx_ (j, q_k);
                        }
                        if constexpr (Form::XprBits & fe_assembler_flags::compute_physical_quad_nodes) {
                            bs_packet.quad_node_id = local_cell_id * n_quadrature_nodes + q_k;
                        }
                        value += Base::quad_weights_[q_k] * form_(bs_packet);
                    }
                    triplet_list.emplace_back(
                      test_active_dofs[j], is_galerkin ? test_active_dofs[i] : trial_active_dofs[i],
                      value * fe_packet.cell_measure);
                }
            }
        }
	return;
    }
    constexpr int n_dofs() const { return trial_dof_handler()->n_dofs(); }
    constexpr int rows() const { return test_dof_handler()->n_dofs(); }
    constexpr int cols() const { return trial_dof_handler()->n_dofs(); }
    constexpr const TrialSpace& trial_space() const { return *trial_space_; }
};

}   // namespace internals
}   // namespace fdapde

#endif   // __BS_BILINEAR_FORM_ASSEMBLER_H__
