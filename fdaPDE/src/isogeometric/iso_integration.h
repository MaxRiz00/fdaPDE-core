#ifndef __ISO_INTEGRATION_H__
#define __ISO_INTEGRATION_H__

#include "header_check.h"


namespace fdapde {
namespace internals{

struct iso_quadrature_gauss_base { };
template <typename T>
concept is_iso_quadrature_gauss = std::is_base_of_v<iso_quadrature_gauss_base, T>;
template <typename T> [[maybe_unused]] static constexpr bool is_iso_quadrature_gauss_v = is_iso_quadrature_gauss<T>;

// Higher degree quadrature
template <typename LhsQuadrature, typename RhsQuadrature>
requires(
    LhsQuadrature::local_dim == RhsQuadrature::local_dim && is_iso_quadrature_gauss_v<LhsQuadrature> &&
    is_iso_quadrature_gauss_v<RhsQuadrature>)
struct higher_degree_iso_quadrature :
    std::type_identity<
    std::conditional_t<(LhsQuadrature::order > RhsQuadrature::order), LhsQuadrature, RhsQuadrature>> { };

template <typename LhsQuadrature, typename RhsQuadrature>
using higher_degree_iso_quadrature_t = higher_degree_iso_quadrature<LhsQuadrature, RhsQuadrature>::type;

// quadrature points and weights at: https://people.sc.fsu.edu/~jburkardt/datasets/datasets.html
template <int LocalDim, int Size> struct iso_quadrature_gauss_legendre;

// 1D 1 point formula
template <> struct iso_quadrature_gauss_legendre<1, 1> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 1;
    static constexpr int degree = 2 * order - 1;   // 1

    static constexpr Vector<double, order> nodes {
        std::array<double, order> {0.000000000000000}
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {2.000000000000000}
    };
};

// 1D 2 point formula
template <> struct iso_quadrature_gauss_legendre<1, 2> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 2;
    static constexpr int degree = 2 * order - 1;   // 3

    static constexpr Vector<double, order> nodes {
        std::array<double, order> {-0.5773502691896257, 0.5773502691896257}
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> { 0.9999999999999998, 0.9999999999999998}
    };
};

// 1D 3 point formula
template <> struct iso_quadrature_gauss_legendre<1, 3> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 1;
    static constexpr int order  = 3;
    static constexpr int degree = 2 * order - 1;   // 5

    static constexpr Vector<double, order> nodes {
        std::array<double, order> {-0.7745966692414834, 0.000000000000000, 0.7745966692414834}
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556}
    };
};

// 2D 1 point formula
template <> struct iso_quadrature_gauss_legendre<2, 1> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 2;
    static constexpr int order = 1;
    static constexpr int degree = 2 * order - 1;   // 1

    static constexpr Matrix<double, order, local_dim> nodes {
        std::array<double, order * local_dim> {
            0.000000000000000, 0.000000000000000
        }
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {4.000000000000000}
    };
};

// 2D 4 point formula
template <> struct iso_quadrature_gauss_legendre<2, 4> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 2;
    static constexpr int order = 4;
    static constexpr int degree = 2 * order - 1;   // 7

    static constexpr Matrix<double, order, local_dim> nodes {
        std::array<double, order * local_dim> {
            -0.5773502691896257, -0.5773502691896257, 
            -0.5773502691896257, 0.5773502691896257,
            0.5773502691896257, -0.5773502691896257,
            0.5773502691896257, 0.5773502691896257
        }
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {
            0.9999999999999996, 0.9999999999999996, 0.9999999999999996, 0.9999999999999996
        }
    };
};

// 2D 9 point formula
template <> struct iso_quadrature_gauss_legendre<2, 9> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 2;
    static constexpr int order = 9;
    static constexpr int degree = 2 * order - 1;   // 17

    static constexpr Matrix<double, order, local_dim> nodes {
        std::array<double, order * local_dim> {
            -0.7745966692414835, -0.7745966692414835,
            -0.7745966692414835,  0.0000000000000000,
            -0.7745966692414835,  0.7745966692414835,
            0.0000000000000000, -0.7745966692414835,
            0.0000000000000000,  0.0000000000000000,
            0.0000000000000000,  0.7745966692414835,
            0.7745966692414835, -0.7745966692414835,
            0.7745966692414835,  0.0000000000000000,
            0.7745966692414835,  0.7745966692414835
        }
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {
            0.3086419753086420, 0.4938271604938272, 0.3086419753086420, 0.4938271604938272,
            0.7901234567901235, 0.4938271604938272, 0.3086419753086420, 0.4938271604938272,
            0.3086419753086420
        }
    };
};

// 3D 1 point formula
template <> struct iso_quadrature_gauss_legendre<3, 1> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 3;
    static constexpr int order = 1;
    static constexpr int degree = 2 * order - 1;   // 2

    static constexpr Matrix<double, order, local_dim> nodes {
        std::array<double, order * local_dim> {
            0.000000000000000, 0.000000000000000, 0.000000000000000
        }
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {8.000000000000000}
    };
};

// 3D 8 point formula
template <> struct iso_quadrature_gauss_legendre<3, 8> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 3;
    static constexpr int order = 8;
    static constexpr int degree = 2 * order - 1;   // 15

static constexpr Matrix<double, order, local_dim> nodes {
    std::array<double, order * local_dim> {
        -0.5773502691896257, -0.5773502691896257, -0.5773502691896257,
            0.5773502691896257, -0.5773502691896257, -0.5773502691896257,
        -0.5773502691896257,  0.5773502691896257, -0.5773502691896257,
            0.5773502691896257,  0.5773502691896257, -0.5773502691896257,
        -0.5773502691896257, -0.5773502691896257,  0.5773502691896257,
            0.5773502691896257, -0.5773502691896257,  0.5773502691896257,
        -0.5773502691896257,  0.5773502691896257,  0.5773502691896257,
            0.5773502691896257,  0.5773502691896257,  0.5773502691896257
    }
};
static constexpr Vector<double, order> weights {
    std::array<double, order> {
        0.9999999999999996, 0.9999999999999996, 0.9999999999999996, 0.9999999999999996,
        0.9999999999999996, 0.9999999999999996, 0.9999999999999996, 0.9999999999999996
    }
};
};

// 3D 27 point formula
template <> struct iso_quadrature_gauss_legendre<3, 27> : public iso_quadrature_gauss_base {
    static constexpr int local_dim = 3;
    static constexpr int order = 27;
    static constexpr int degree = 2 * order - 1;   

    static constexpr Matrix<double, order, local_dim> nodes {
        std::array<double, order * local_dim> {
            -0.5773502691896257, -0.5773502691896257, -0.5773502691896257,
            0.5773502691896257, -0.5773502691896257, -0.5773502691896257,
            -0.5773502691896257,  0.5773502691896257, -0.5773502691896257,
            0.5773502691896257,  0.5773502691896257, -0.5773502691896257,
            -0.5773502691896257, -0.5773502691896257,  0.5773502691896257,
            0.5773502691896257, -0.5773502691896257,  0.5773502691896257,
            -0.5773502691896257,  0.5773502691896257,  0.5773502691896257,
            0.5773502691896257,  0.5773502691896257,  0.5773502691896257
        }
    };
    static constexpr Vector<double, order> weights {
        std::array<double, order> {
            0.9999999999999996, 0.9999999999999996, 0.9999999999999996, 0.9999999999999996,
            0.9999999999999996, 0.9999999999999996, 0.9999999999999996, 0.9999999999999996
        }
    };
};


}// namespace internals

// 1D formulas ( da chiedere se usare quella di spline, altrimenti ci sono conflitti)
//[[maybe_unusued]] static struct QGL1DP1_ : internals::iso_quadrature_gauss_legendre<1, 1> { } QGL1DP1;
//[[maybe_unusued]] static struct QGL1DP2_ : internals::iso_quadrature_gauss_legendre<1, 2> { } QGL1DP2;
//[[maybe_unusued]] static struct QGL1DP3_ : internals::iso_quadrature_gauss_legendre<1, 3> { } QGL1DP3;
// 2D formulas
[[maybe_unused]] static struct QGL2DP1_ : internals::iso_quadrature_gauss_legendre<2, 1> { } QGL2DP1;
[[maybe_unused]] static struct QGL2DP4_ : internals::iso_quadrature_gauss_legendre<2, 4> { } QGL2DP4;
[[maybe_unused]] static struct QGL2DP9_ : internals::iso_quadrature_gauss_legendre<2, 9> { } QGL2DP9;
// 3D formulas
[[maybe_unused]] static struct QGL3DP1_ : internals::iso_quadrature_gauss_legendre<3, 1> { } QGL3DP1;
[[maybe_unused]] static struct QGL3DP8_ : internals::iso_quadrature_gauss_legendre<3, 8> { } QGL3DP8;
[[maybe_unused]] static struct QGL3DP27_ : internals::iso_quadrature_gauss_legendre<3, 27> { } QGL3DP27;

}


#endif // __ISO_INTEGRATION_H__