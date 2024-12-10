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

#ifndef __BS_ASSEMBLER_BASE_H__
#define __BS_ASSEMBLER_BASE_H__
  
namespace fdapde{

template <typename Derived_> struct BsMap;
  
enum bs_assembler_flags {
    compute_shape_values        = 0x0001,
    compute_shape_dx            = 0x0002,
    compute_shape_ddx           = 0x0004,
    compute_physical_quad_nodes = 0x0010,
    compute_cell_diameter       = 0x0020
};
  
namespace internals {

// informations sent from the assembly loop to the integrated forms
template <int LocalDim> struct bs_assembler_packet {
    static constexpr int local_dim = LocalDim;
    bs_assembler_packet() = default;
    bs_assembler_packet(bs_assembler_packet&&) noexcept = default;
    bs_assembler_packet(const bs_assembler_packet&) noexcept = default;

    // geometric informations
    int quad_node_id;       // active physical quadrature node index
    double cell_measure;    // active cell measure
    double cell_diameter;   // active cell diameter

    // functional informations
    double trial_value, test_value;   // \psi_i(q_k), \psi_j(q_k)
    double trial_dx, test_dx;         // d{\psi_i}/dx(q_k), d{\psi_j}/dx(q_k)
    double trial_ddx, test_ddx;       // d^2{\psi_i}/dx^2(q_k), d^2{\psi_j}/dx^2(q_k)
};

// base class for spline-based assembly loops
template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
struct bs_assembler_base {
    fdapde_static_assert(sizeof...(Quadrature_) < 2, YOU_CAN_SUPPLY_AT_MOST_ONE_QUADRATURE_RULE_TO_A_BS_ASSEMBLY_LOOP);
    // detect test space (since a test function is always present in a weak form)
    using TestSpace = test_space_t<Form_>;
    using Form =
      std::decay_t<decltype(meta::xpr_wrap<BsMap, decltype([]<typename Xpr>() {
	    return !(
	        std::is_invocable_v<Xpr, bs_assembler_packet<Xpr::StaticInputSize>>);
	  })>(std::declval<Form_>()))>;
    using Triangulation = typename std::decay_t<Triangulation_>;
    static constexpr int local_dim = Triangulation::local_dim;
    static constexpr int embed_dim = Triangulation::embed_dim;
    static constexpr int Options = Options_;
    using FunctionSpace = TestSpace;
    using DofHandlerType = DofHandler<local_dim, embed_dim, fdapde::bspline>;
    using Quadrature = decltype([]() {
        if constexpr (sizeof...(Quadrature_) == 0) {
            return void();   // quadrature selcted at run-time provided the actual order of spline basis
        } else {
            return std::get<0>(std::tuple<Quadrature_...>());   // user-defined quadrature
        }
    }());
    using geo_iterator = typename Triangulation::cell_iterator;
    using dof_iterator = typename DofHandlerType::cell_iterator;

    bs_assembler_base() = default;
    bs_assembler_base(
      const Form_& form, typename fe_traits::geo_iterator begin, typename fe_traits::geo_iterator end,
      const Quadrature_&... quadrature)
        requires(sizeof...(quadrature) <= 1):
        form_(meta::xpr_wrap<BsMap, decltype([]<typename Xpr>() {
	      return !(
	          std::is_invocable_v<Xpr, fe_assembler_packet<Xpr::StaticInputSize>>);
	    })>(form)),
        dof_handler_(std::addressof(internals::test_space(form_).dof_handler())),
        test_space_ (std::addressof(internals::test_space(form_))),
        begin_(begin),
        end_(end) {
        fdapde_assert(dof_handler_->n_dofs() > 0);
        // copy quadrature rule
        if constexpr (sizeof...(quadrature) == 1) {
            auto quad_rule = std::get<0>(std::make_tuple(quadrature...));
            for (int i = 0; i < quad_rule.order; ++i) {
                quad_nodes_  [i] = quad_rule.nodes  [i];
                quad_weights_[i] = quad_rule.weights[i];
            }
        } else {
            internals::get_bs_quadrature(test_space_->order(), quad_nodes_, quad_weights_);
        }
    }
    const TestSpace& test_space() const { return *test_space_; }
   protected:
    // evaluation of \psi_i(q_j), i = 1, ..., n_basis, j = 1, ..., n_quadrature_nodes
    template <typename BasisType__>
    MdArray<double, MdExtents<Dynamic, Dynamic>> eval_shape_values(BasisType__&& basis) {
        using BasisType = std::decay_t<BasisType__>;
	int n_quadrature_nodes = quad_nodes_.rows();
        int n_basis = basis.size();

        MdArray<double, MdExtents<Dynamic, Dynamic>> shape_values_(n_basis, n_quadrature_nodes);
        for (int i = 0; i < n_basis; ++i) {
            // evaluation of \psi_i at q_j, j = 1, ..., n_quadrature_nodes
            for (int j = 0; j < n_quadrature_nodes; ++j) {
                shape_values_(i, j) = basis[i](quad_nodes_.row(j).transpose());
            }
        }
        return shape_values_;
    }
    // evaluation of k-th order derivative of basis function
    template <typename BasisType__>
        requires(requires(BasisType__ basis, int i, int k) { basis[i].gradient(k); })
    MdArray<double, MdExtents<Dynamic, Dynamic>> eval_shape_derivative(BasisType__&& basis, int order) {
        using BasisType = std::decay_t<BasisType__>;
	using DerivativeType = decltype(std::declval<BasisType>()[std::declval<int>()].gradient(std::declval<int>()));
        int n_quadrature_nodes = quad_nodes_.rows();
        int n_basis = basis.size();
	// instantiate derivative functors once
        std::vector<DerivativeType> der(n_basis);
        for (int i = 0; i < n_basis; ++i) { ders[i] = basis[i].gradient(order); }

        MdArray<double, MdExtents<Dynamic, Dynamic>> shape_derivatives_(n_basis, n_quadrature_nodes);
        for (int i = 0; i < n_basis; ++i) {
            // evaluation of d^k\psi_i/dx^k at q_j, j = 1, ..., n_quadrature_nodes
            for (int j = 0; j < n_quadrature_nodes; ++j) {
                shape_derivatives_(i, j) = der[i](quad_nodes_.row(j).transpose());
            }
        }
        return shape_derivatives_;
    }
    template <typename BasisType__> MdArray<double, MdExtents<Dynamic, Dynamic>> eval_shape_dx (BasisType__&& basis) {
        return eval_shape_derivatives(basis, /* order = */ 1);
    }
    template <typename BasisType__> MdArray<double, MdExtents<Dynamic, Dynamic>> eval_shape_ddx(BasisType__&& basis) {
        return eval_shape_derivatives(basis, /* order = */ 2);
    }
    void distribute_quadrature_nodes(
      std::unordered_map<const void*, Eigen::Matrix<double, Dynamic, Dynamic>>& bs_map_buff, dof_iterator begin,
      dof_iterator end) {
        Eigen::Matrix<double, Dynamic, Dynamic> quad_nodes;
        int n_quadrature_nodes = quad_nodes_.rows();
        quad_nodes.resize(n_quadrature_nodes * (end_.begin() - begin_.index()), embed_dim);
        int local_cell_id = 0;
        for (geo_iterator it = begin_; it != end_; ++it) {
            for (int q_k = 0; q_k < n_quadrature_nodes; ++q_k) {
                quad_nodes.row(local_cell_id * n_quadrature_nodes + q_k) =
                  it->J() * quad_nodes_.row(q_k).transpose() + it->node(0);
            }
            local_cell_id++;
        }
        // evaluate Map nodes at quadrature nodes
        meta::xpr_apply_if<
          decltype([]<typename Xpr_, typename... Args>(Xpr_& xpr, Args&&... args) {
              xpr.init(std::forward<Args>(args)...);
              return;
          }),
          decltype([]<typename Xpr_>() {
              return requires(Xpr_ xpr) { xpr.init(bs_map_buff, quad_nodes, begin, end); };
          })>(form_, bs_map_buff, quad_nodes, begin, end);
        return;
    }

    Form form_;
    Eigen::Matrix<double, Dynamic, local_dim> quad_nodes_, quad_weights_;
    const DofHandlerType* dof_handler_;
    const TestSpace* test_space_;
    typename fe_traits::geo_iterator begin_, end_;
};

}   // namespace internals
}   // namespace fdapde

#endif   // __BS_ASSEMBLER_BASE_H__
