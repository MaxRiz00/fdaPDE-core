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

#ifndef __BS_INTEGRATION_H__
#define __BS_INTEGRATION_H__

#include "../linear_algebra/constexpr_matrix.h"

namespace fdapde {
namespace internals {

struct bs_quadrature_gauss_base { };
template <typename T>
concept is_bs_quadrature_gauss = std::is_base_of_v<bs_quadrature_gauss_base, T>;
template <typename T> [[maybe_unused]] static constexpr bool is_bs_quadrature_gauss_v = is_bs_quadrature_gauss<T>;

template <typename LhsQuadrature, typename RhsQuadrature>
    requires(
      LhsQuadrature::local_dim == RhsQuadrature::local_dim && is_bs_quadrature_gauss_v<LhsQuadrature> &&
      is_bs_quadrature_gauss_v<RhsQuadrature>)
struct higher_degree_bs_quadrature :
    std::type_identity<
      std::conditional_t<(LhsQuadrature::order > RhsQuadrature::order), LhsQuadrature, RhsQuadrature>> { };
template <typename LhsQuadrature, typename RhsQuadrature>
using higher_degree_bs_quadrature_t = higher_degree_bs_quadrature<LhsQuadrature, RhsQuadrature>::type;

// quadrature points and weights at: https://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
template <int LocalDim, int Size> struct bs_quadrature_gauss_legendre;

// degree: highest polynomial degree correctly integrated by the formula
// order : number of quadrature nodes

// 1D 1 point formula
template <> struct bs_quadrature_gauss_legendre<1, 1> : public bs_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 1;
    static constexpr int degree = 2 * order - 1;   // 1

    static constexpr cexpr::Vector<double, order> nodes {
      std::array<double, order> {0.000000000000000}
    };
    static constexpr cexpr::Vector<double, order> weights {
      std::array<double, order> {2.000000000000000}
    };
};

// 1D 2 point formula
template <> struct bs_quadrature_gauss_legendre<1, 2> : public bs_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 2;
    static constexpr int degree = 2 * order - 1;   // 3

    static constexpr cexpr::Vector<double, order> nodes {
      std::array<double, order> {-0.5773502691896257, 0.5773502691896257}
    };
    static constexpr cexpr::Vector<double, order> weights {
      std::array<double, order> { 0.9999999999999998, 0.9999999999999998}
    };
};

// 1D 3 point formula
template <> struct bs_quadrature_gauss_legendre<1, 3> : public bs_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 3;
    static constexpr int degree = 2 * order - 1;   // 5

    static constexpr cexpr::Vector<double, order> nodes {
      std::array<double, order> {-0.7745966692414834, 0.0000000000000000, 0.7745966692414834}
    };
    static constexpr cexpr::Vector<double, order> weights {
      std::array<double, order> { 0.5555555555555555, 0.8888888888888888, 0.5555555555555555}
    };
};

// 1D 4 point formula
template <> struct bs_quadrature_gauss_legendre<1, 4> : public bs_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 4;
    static constexpr int degree = 2 * order - 1;   // 7

    static constexpr cexpr::Vector<double, order> nodes {
      std::array<double, order> {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526}
    };
    static constexpr cexpr::Vector<double, order> weights {
      std::array<double, order> { 0.3478548451374538,  0.6521451548625461, 0.6521451548625461, 0.3478548451374538}
    };
};

// copy weights and nodes of quadrature rule based on run-time polynomial order
template <typename T>
    requires(requires(T t, int i, int j) {
        { t(i, j) } -> std::same_as<double&>;
        { t.resize(i, j) } -> std::same_as<void>;
    })
void get_bs_quadrature(int order, T& quad_nodes, T& quad_weights) {
    fdapde_assert(order >= 0 && order <= 7);
    auto copy_ = []<typename QuadRule>(QuadRule q, T& quad_nodes_, T& quad_weights_) {
        quad_nodes_  .resize(q.order, q.local_dim);
        quad_weights_.resize(q.order, q.local_dim);
        for (int i = 0; i < q.order; ++i) {
            quad_nodes_  (i, 0) = q.nodes  [i];
            quad_weights_(i, 0) = q.weights[i];
        }
    };
    if (order == 0 || order == 1) { copy_(bs_quadrature_gauss_legendre<1, 1> {}, quad_nodes, quad_weights); }
    if (order >  1 && order <= 3) { copy_(bs_quadrature_gauss_legendre<1, 3> {}, quad_nodes, quad_weights); }
    if (order >  3 && order <= 5) { copy_(bs_quadrature_gauss_legendre<1, 3> {}, quad_nodes, quad_weights); }
    if (order >  5 && order <= 7) { copy_(bs_quadrature_gauss_legendre<1, 4> {}, quad_nodes, quad_weights); }
}

}   // namespace internals

// 1D formulas
[[maybe_unused]] static struct QGL1D1P_ : internals::bs_quadrature_gauss_legendre<1, 1> { } QGL1D1P;
[[maybe_unused]] static struct QGL1D2P_ : internals::bs_quadrature_gauss_legendre<1, 2> { } QGL1D2P;
[[maybe_unused]] static struct QGL1D3P_ : internals::bs_quadrature_gauss_legendre<1, 3> { } QGL1D3P;
[[maybe_unused]] static struct QGL1D4P_ : internals::bs_quadrature_gauss_legendre<1, 4> { } QGL1D4P;
  
}   // namespace fdapde

#endif   // __BS_INTEGRATION_H__
