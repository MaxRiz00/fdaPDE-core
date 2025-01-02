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

#ifndef __ASSEMBLY_H__
#define __ASSEMBLY_H__

#include "fields/meta.h"
#include "finite_elements/fe_integration.h"

namespace fdapde {

[[maybe_unused]] static constexpr int CellMajor = 0;
[[maybe_unused]] static constexpr int FaceMajor = 1;

// test function forward decl
template <typename FunctionSpace_, typename SpaceCategory_> struct TestFunction;
template <typename FunctionSpace_>   // deduction guide
TestFunction(FunctionSpace_ function_space) -> TestFunction<FunctionSpace_, typename FunctionSpace_::space_category>;
// trial function forward decl
template <typename FunctionSpace, typename SpaceCategory_> struct TrialFunction;
template <typename FunctionSpace_>   // deduction guide
TrialFunction(FunctionSpace_ function_space) -> TrialFunction<FunctionSpace_, typename FunctionSpace_::space_category>;

namespace internals {

// set of internal utilities to write weak form assembly loops

// detect trial space from bilinear form
template <typename Xpr> constexpr decltype(auto) trial_space(Xpr&& xpr) {
    constexpr bool found = xpr_any_of<
      decltype([]<typename Xpr_>() { return requires { typename Xpr_::TrialSpace; }; }), std::decay_t<Xpr>>();
    fdapde_static_assert(found, NO_TRIAL_SPACE_FOUND_IN_EXPRESSION);
    return xpr_apply_if<
      decltype([]<typename Xpr_>(Xpr_&& xpr) -> auto& { return xpr.function_space(); }),
      decltype([]<typename Xpr_>() { return requires { typename Xpr_::TrialSpace; }; })>(std::forward<Xpr>(xpr));
}
template <typename Xpr> using trial_space_t = std::decay_t<decltype(trial_space(std::declval<Xpr>()))>;
// detect test space from bilinear form
template <typename Xpr> constexpr decltype(auto)  test_space(Xpr&& xpr) {
    constexpr bool found = xpr_any_of<
      decltype([]<typename Xpr_>() { return requires { typename Xpr_::TestSpace; }; }), std::decay_t<Xpr>>();
    fdapde_static_assert(found, NO_TEST_SPACE_FOUND_IN_EXPRESSION);
    return xpr_apply_if<
      decltype([]<typename Xpr_>(Xpr_&& xpr) -> auto& { return xpr.function_space(); }),
      decltype([]<typename Xpr_>() { return requires { typename Xpr_::TestSpace; }; })>(std::forward<Xpr>(xpr));
}
template <typename Xpr> using test_space_t = std::decay_t<decltype(test_space(std::declval<Xpr>()))>;

template <typename T>
    requires(
      std::is_floating_point_v<T> ||
      requires(T t, int k) {
          { t.size() } -> std::convertible_to<std::size_t>;
          { t.operator[](k) } -> std::convertible_to<double>;
      })
constexpr auto scalar_or_kth_component_of(const T& t, std::size_t k) {
    if constexpr (std::is_floating_point_v<T>) {
        return t;
    } else {
        fdapde_constexpr_assert(k < t.size());
        return t[k];
    }
}
  
// arithmetic between matrix assembly loops
template <typename Derived> struct assembly_xpr_base {
    constexpr const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

template <typename Lhs, typename Rhs>
struct assembly_add_op : public assembly_xpr_base<assembly_add_op<Lhs, Rhs>> {
    fdapde_static_assert(
      std::is_same_v<decltype(std::declval<Lhs>().assemble()) FDAPDE_COMMA decltype(std::declval<Rhs>().assemble())>,
      YOU_ARE_SUMMING_NON_COMPATIBLE_ASSEMBLY_LOOPS);
    fdapde_static_assert(
      std::is_same_v<typename Lhs::space_category FDAPDE_COMMA typename Rhs::space_category>,
      CANNOT_SUM_ASSEMBLY_LOOPS_OF_DIFFERENT_CATEGORY);
    using OutputType = decltype(std::declval<Lhs>().assemble());
    assembly_add_op(const Lhs& lhs, const Rhs& rhs) : lhs_(lhs), rhs_(rhs) {
        fdapde_assert(lhs.rows() == rhs.rows() && lhs.cols() == rhs.cols());
    }
    OutputType assemble() const {
        if constexpr (std::is_same_v<OutputType, Eigen::SparseMatrix<double>>) {
            Eigen::SparseMatrix<double> assembled_mat(lhs_.rows(), lhs_.cols());
            std::vector<Eigen::Triplet<double>> triplet_list;
            assemble(triplet_list);
            // linearity of the integral is implicitly used here, as duplicated triplets are summed up (see Eigen docs)
            assembled_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
            assembled_mat.makeCompressed();
            return assembled_mat;
        } else {
            Eigen::Matrix<double, Dynamic, 1> assembled_vec(lhs_.rows());
            assembled_vec.setZero();
            assemble(assembled_vec);
            return assembled_vec;
        }
    }
    template <typename T> void assemble(T& assembly_buff) const {
        lhs_.assemble(assembly_buff);
        rhs_.assemble(assembly_buff);
        return;
    }
    constexpr int n_dofs() const { return lhs_.n_dofs(); }
    constexpr int rows() const { return lhs_.rows(); }
    constexpr int cols() const { return lhs_.cols(); }
   private:
    Lhs lhs_;
    Rhs rhs_;
};

template <typename Lhs, typename Rhs>
assembly_add_op<Lhs, Rhs> operator+(const assembly_xpr_base<Lhs>& lhs, const assembly_xpr_base<Rhs>& rhs) {
    return assembly_add_op<Lhs, Rhs>(lhs.derived(), rhs.derived());
}

// generic integration loop to integrate scalar expressions over physical domains
template <typename Triangulation_, typename Xpr_, int Options_, typename... Quadrature_> class std_integration_loop {
    fdapde_static_assert(is_scalar_field_v<Xpr_>, THIS_CLASS_IS_FOR_SCALAR_FIELDS_ONLY);
    using Triangulation = std::decay_t<Triangulation_>;
    static constexpr int local_dim = Triangulation::local_dim;
    static constexpr int embed_dim = Triangulation::embed_dim;
    static constexpr int Options = Options_;
    using iterator = std::conditional_t<
      Options == CellMajor, typename Triangulation::cell_iterator, typename Triangulation::boundary_iterator>;

    template <typename... Quad_> struct empty_quadrature {
        static constexpr int order = 0;
    };
    template <typename... Quad_>
    using maybe_empty_quadrature = decltype([]() {
        if constexpr (sizeof...(Quadrature_) == 0) {
            return empty_quadrature<Quadrature_...> {};
        } else {
            return std::tuple_element_t<0, std::tuple<Quadrature_...>> {};
        }
    }());
    using Quadrature = maybe_empty_quadrature<Quadrature_...>;
  
    Xpr_ xpr_;
    iterator begin_, end_;
    Quadrature quadrature_;
   public:
    std_integration_loop() = default;
    std_integration_loop(const Xpr_& xpr, iterator begin, iterator end, const Quadrature_&... quadrature) :
        xpr_(xpr), begin_(begin), end_(end), quadrature_(quadrature...) { }

    double operator()() const {
        double integral_ = 0;
        if constexpr (Quadrature::order == 0) {
            fdapde_static_assert(false, THIS_METHOD_REQUIRES_A_QUADRATURE_RULE);
        } else {
            fdapde_static_assert(
              Triangulation::local_dim == Quadrature::local_dim &&
                (Triangulation::embed_dim == 1 ||
                 internals::is_fe_quadrature_simplex_v<std::tuple_element_t<0, std::tuple<Quadrature_...>>>),
              INVALID_QUADRATURE_RULE);
            constexpr int n_quadrature_nodes = Quadrature::order;
            Eigen::Map<const Eigen::Matrix<double, n_quadrature_nodes, Triangulation::local_dim, Eigen::RowMajor>>
              ref_quad_nodes(quadrature_.nodes.data());
            for (iterator it = begin_; it != end_; ++it) {
                double partial = 0;
                for (int q_k = 0; q_k < n_quadrature_nodes; ++q_k) {
                    partial +=
                      xpr_(it->J() * ref_quad_nodes.row(q_k).transpose() + it->node(0)) * quadrature_.weights[q_k];
                }
                integral_ += partial * it->measure();
            }
        }
        return integral_;
    }
};

// main assembly loop dispatching type. This class instantiates the correctly assembly loop for the form to integrate
template <typename Triangulation, int Options, typename... Quadrature> class integrator_dispatch {
    fdapde_static_assert(sizeof...(Quadrature) < 2, THIS_CLASS_ACCEPTS_AT_MOST_ONE_QUADRATURE_RULE);
    std::tuple<Quadrature...> quadrature_;
    std::conditional_t<
      Options == CellMajor, typename Triangulation::cell_iterator, typename Triangulation::boundary_iterator>
      begin_, end_;
   public:
    integrator_dispatch() = default;
    template <typename Iterator>
    integrator_dispatch(const Iterator& begin, const Iterator& end, const Quadrature&... quadrature) :
        begin_(begin), end_(end), quadrature_(std::make_tuple(quadrature...)) { }

    template <typename Form> auto operator()(const Form& form) const {
        static constexpr bool trial_space_detected = xpr_any_of<
          decltype([]<typename Xpr_>() { return requires { typename Xpr_::TrialSpace; }; }), std::decay_t<Form>>();
        static constexpr bool test_space_detected  = xpr_any_of<
          decltype([]<typename Xpr_>() { return requires { typename Xpr_::TestSpace;  }; }), std::decay_t<Form>>();

        if constexpr (trial_space_detected && test_space_detected) {   // discretizing bilinear form
            using TrialSpace = std::decay_t<decltype(trial_space(form))>;
            using TestSpace  = std::decay_t<decltype(test_space (form))>;
            fdapde_static_assert(
              std::is_same_v<typename TrialSpace::space_category FDAPDE_COMMA typename TestSpace::space_category>,
              TRIAL_AND_TEST_SPACES_ARE_DEFINED_ON_DIFFERENT_SPACE_CATEGORIES);
            if constexpr (sizeof...(Quadrature) == 0) {
                return typename TestSpace::template bilinear_form_assembly_loop<Triangulation, Form, Options> {
                  form, begin_, end_};
            } else {
                return typename TestSpace::template bilinear_form_assembly_loop<
                  Triangulation, Form, Options, std::tuple_element_t<0, std::tuple<Quadrature...>>> {
                  form, begin_, end_, std::get<0>(quadrature_)};
            }
        } else {
            if constexpr (test_space_detected) {   // discretizing linear form
                using TestSpace = std::decay_t<decltype(test_space(form))>;
                if constexpr (sizeof...(Quadrature) == 0) {
                    return typename TestSpace::template linear_form_assembly_loop<Triangulation, Form, Options> {
                      form, begin_, end_};
                } else {
                    return typename TestSpace::template linear_form_assembly_loop<
                      Triangulation, Form, Options, std::tuple_element_t<0, std::tuple<Quadrature...>>> {
                      form, begin_, end_, std::get<0>(quadrature_)};
                }
            } else {
                // integration of callable on domain, directly return the approximated integral
                return internals::std_integration_loop<Triangulation, Form, Options, Quadrature...>(
                  form, begin_, end_, Quadrature {}...)();
            }	    
        }
    }
};

}   // namespace internals

// main entry points for operator discretization
template <typename Triangulation, typename... Quadrature>
auto integral(const Triangulation& triangulation, Quadrature... quadrature) {
    return internals::integrator_dispatch<Triangulation, CellMajor, Quadrature...>(
      triangulation.cells_begin(), triangulation.cells_end(), quadrature...);
}
template <typename Triangulation, typename... Quadrature>
auto integral(
  const CellIterator<Triangulation>& begin, const CellIterator<Triangulation>& end, Quadrature... quadrature) {
    return internals::integrator_dispatch<Triangulation, CellMajor, Quadrature...>(begin, end, quadrature...);
}
// boundary integration
template <typename Triangulation, typename... Quadrature>
auto integral(
  const BoundaryIterator<Triangulation>& begin, const BoundaryIterator<Triangulation>& end, Quadrature... quadrature) {
    return internals::integrator_dispatch<Triangulation, FaceMajor, Quadrature...>(begin, end, quadrature...);
}
template <typename Triangulation, typename... Quadrature>
auto integral(
  const std::pair<BoundaryIterator<Triangulation>, BoundaryIterator<Triangulation>>& range, Quadrature... quadrature) {
    return internals::integrator_dispatch<Triangulation, FaceMajor, Quadrature...>(
      range.first, range.second, quadrature...);
}
  
}   // namespace fdapde

#endif   // __ASSEMBLY_H__
