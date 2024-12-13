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

#ifndef __BS_SPACE_H__
#define __BS_SPACE_H__

#include "../utils/symbols.h"
#include "bs_assembler_base.h"
#include "dof_handler.h"

namespace fdapde {

// forward declarations
template <typename BsSpace_> class BsFunction;
namespace internals {

template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
class bs_bilinear_form_assembly_loop;
template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
class bs_linear_form_assembly_loop;

}   // namespace internals
  
template <typename Triangulation_> class BsSpace {
    fdapde_static_assert(
      Triangulation_::local_dim == 1 && Triangulation_::embed_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    template <typename T> struct subscript_t_impl {
        using type = std::decay_t<decltype(std::declval<T>().operator[](std::declval<int>()))>;
    };
    template <typename T> using subscript_t = typename subscript_t_impl<T>::type;
   public:
    using Triangulation = std::decay_t<Triangulation_>;
    static constexpr int local_dim = Triangulation::local_dim;
    static constexpr int embed_dim = Triangulation::embed_dim;
    using BasisType = SplineBasis;
    using ShapeFunctionType = subscript_t<BasisType>;
    using DofHandlerType = DofHandler<local_dim, embed_dim, fdapde::bspline>;
    using space_category = fdapde::bspline;
    template <typename Triangulation__, typename Form__, int Options__, typename... Quadrature__>
    using bilinear_form_assembly_loop =
      internals::bs_bilinear_form_assembly_loop<Triangulation__, Form__, Options__, Quadrature__...>;
    template <typename Triangulation__, typename Form__, int Options__, typename... Quadrature__>
    using linear_form_assembler_loop =
      internals::bs_linear_form_assembly_loop  <Triangulation__, Form__, Options__, Quadrature__...>;

    BsSpace() = default;
    BsSpace(const Triangulation_& interval, int order) :
        triangulation_(std::addressof(interval)), dof_handler_(interval), order_(order) {
        dof_handler_.enumerate(BasisType(interval, order));
	// build reference [-1, 1] interval with nodes mapped from physical interval [a, b]
        Eigen::Matrix<double, Dynamic, 1> ref_nodes(triangulation_->n_nodes());
        for (int i = 0; i < triangulation_->n_nodes(); ++i) {
            ref_nodes[i] = map_to_reference(triangulation_->nodes()(i, 0));
	    std::cout << triangulation_->nodes()(i, 0) << " : " << ref_nodes[i] << " -- ";
        }
	std::cout << std::endl;
        // generate basis on reference [-1, 1] interval
        basis_ = BasisType(Triangulation(ref_nodes), order);
    }
    // observers
    const Triangulation& triangulation() const { return *triangulation_; }
    const DofHandlerType& dof_handler() const { return dof_handler_; }
    DofHandlerType& dof_handler() { return dof_handler_; }
    constexpr int n_shape_functions() const { return basis_.size(); }
    constexpr int n_shape_functions_face() const { return 1; }
    int n_dofs() const { return dof_handler_.n_dofs(); }
    const BasisType& basis() const { return basis_; }
    int order() const { return order_; }
    // evaluation
    template <typename InputType>
        requires(std::is_invocable_v<ShapeFunctionType, InputType>)
    constexpr auto eval_shape_value(int i, const InputType& p) const {
        return basis_[i](p);
    }
    template <typename InputType>
        requires(std::is_invocable_v<decltype(std::declval<ShapeFunctionType>().gradient(1)), InputType>)
    constexpr auto eval_shape_dx(int i, const InputType& p) const {
        return basis_[i].gradient(1)(p);
    }
    template <typename InputType>
        requires(std::is_invocable_v<decltype(std::declval<ShapeFunctionType>().gradient(2)), InputType>)
    constexpr auto eval_shape_ddx(int i, const InputType& p) const {
        return basis_[i].gradient(2)(p);
    }
    template <typename InputType>
    constexpr auto eval_face_shape_value(int i, [[maybe_unused]] const InputType& p) const {
        return (i == 0 || i == n_dofs() - 1) ? 1.0 : 0.0;
    }
    template <typename InputType> constexpr auto eval_face_shape_dx (int i, const InputType& p) const {
        return basis_[i].gradient(1)(p);
    }
    template <typename InputType> constexpr auto eval_face_shape_ddx(int i, const InputType& p) const {
        return basis_[i].gradient(2)(p);
    }
    // evaluation on physical domain
    // evaluate value of the i-th shape function defined on physical domain [a, b]
    template <typename InputType> auto eval_cell_value(int i, const InputType& p) const {
        double p_;
        if constexpr (fdapde::is_subscriptable<InputType, int>) { p_ = p[0]; }
	else { p_ = p; }
        if (p_ < triangulation_->range()[0] || p_ > triangulation_->range()[1]) {
            return std::numeric_limits<double>::quiet_NaN();   // return NaN if point lies outside domain
        }
        return eval_shape_value(i, map_to_reference(p_));
    }

    // need to return something which represent a basis function on the whole physical domain

    // generate fe_function bounded to this finite element space
    /* FeFunction<FeSpace<Triangulation_, FeType_>> make_fe_function(const Eigen::Matrix<double, Dynamic, 1>& coeff_vec) { */
    /*     return FeFunction<FeSpace<Triangulation_, FeType_>>(*this, coeff_vec); */
    /* } */
   private:
    // linear mapping from p \in [a, b] to \hat p \in [-1, +1]
    double map_to_reference(double p) {
        double a = triangulation_->range()[0], b = triangulation_->range()[1];
        return (p - (b - a) / 2) * 2 / (b + a);
    };

    const Triangulation* triangulation_;
    DofHandlerType dof_handler_;   // dof_handler over physical domain
    BasisType basis_;              // basis_ over reference interval [-1, +1]
    int order_;                    // spline order
};

}   // namespace fdapde

#endif   // __BS_SPACE_H__
