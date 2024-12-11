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

#ifndef __OLDSPLINE_H__
#define __OLDSPLINE_H__

#include <limits>
#include "../../fields/scalar_field.h"
//#include "../../fields/scalar_expressions.h"

namespace fdapde {
namespace core {

// Let u_0, u_1, ..., u_n n distinct knots. Call U = [u_0, u_1, ..., u_n] the knot vector.
// By Cox-DeBoor formula the i-th spline basis of order j N_{ij} is recursively defined as
//
// N_i0(x) = 1 if x \in [u_i, u_i+1) 0 otherwise
// N_ij(x) = [(x-u_i)/(u_i+j - u_i)]*N_i,j-1(x) + [(u_i+j+1 - x)/(u_i+j+1 - u_i+1)]*N_i+1,j-1(x)

// A spline of order R centered in knot u_i.
// Template parameter M is used to keep track of the order of the spline while developing the Cox-DeBoor recursion.
template <int R, int M = R> class OldSpline : public ScalarBase<1, OldSpline<R, M>> {
   private:
    DVector<double> knots_ {};
    int i_ = 0;              // knot index where this basis element is centered
    double a_ = 0, b_ = 0;   // constants a_ = 1/(u_i+j - u_i), b_ = 1/(u_i+j+1 - u_i+1)
   public:
        using Base = ScalarBase<1, Spline>;
        using Scalar = double;
        static constexpr int StaticInputSize = 1;
        using InputType = double; //cexpr::Vector<Scalar, StaticInputSize>;
        
    static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporariess
    OldSpline() = default;
    OldSpline(const DVector<double>& knots, int i) : knots_(knots), i_(i) {
        // avoid possible divisions by zero
        a_ = knots_[i_ + R] - knots_[i_] != 0 ? 1.0 / (knots_[i_ + R] - knots_[i_]) : 0.0;
        b_ = knots_[i_ + R + 1] - knots_[i_ + 1] != 0 ? 1.0 / (knots_[i_ + R + 1] - knots_[i_ + 1]) : 0.0;
    };
    // evaluates the spline at a given point by Cox-DeBoor recursion
    inline double operator()(double x) const {
        return a_ * (x  - knots_[i_]) * OldSpline<R - 1, M>(knots_, i_)(x) +
               b_ * (knots_[i_ + R + 1] - x) * OldSpline<R - 1, M>(knots_, i_ + 1)(x);
    }
    // compute derivative of order K as a ScalarExpr
    // d^K/dx^K N_ij(x) = j/(u_i+j - u_i)*[d^{K-1}/dx^{K-1} N_i,j-1(x)]  - j/(u_i+j+1 - u_i+1)*
    // [d^{K-1}/dx^{K-1} N_i+1,j-1(x)]
    template <int K> auto derive() const {
        if constexpr (K == 1) {   // end of recursion
            return (R * a_) * OldSpline<R - 1, M>(knots_, i_) - (R * b_) * OldSpline<R - 1, M>(knots_, i_ + 1);
        } else {   // exploit Cox-DeBoor recursive formula
            return (R * a_) * OldSpline<R - 1, M>(knots_, i_).template derive<K - 1>() -
                   (R * b_) * OldSpline<R - 1, M>(knots_, i_ + 1).template derive<K - 1>();
        }
    }
};

// partial template specialization for order 0 splines (end of recursion)
template <int M> class OldSpline<0, M> : public ScalarBase<1, OldSpline<0, M>> {
   private:
    DVector<double> knots_ {};
    int i_ = 0;   // knot index where this basis is centered
    static constexpr double tol_ = 50 * std::numeric_limits<double>::epsilon();   // approx 10^-14
   public:
    static constexpr int NestAsRef = 0;   // avoid nesting as reference
    OldSpline() = default;
    OldSpline(const DVector<double>& knots, int i) : knots_(knots), i_(i) { }
    // indicator function over the interval [u_i, u_{i+1}). Returns 1 if at the end of the knot span.
    inline double operator()(double x) const {
        return (knots_[i_] <= x && x < knots_[i_ + 1]) ||
               (std::abs(x - knots_[knots_.rows() - 1]) < tol_ && i_ == knots_.rows() - M - 2) ? 1.0 : 0.0;
    }
    // return derivative as a zero callable field
    auto derive() const { return 0; } // Zerofield 1 
};

}   // namespace core
}   // namespace fdapde

#endif   // __OLDSPLINE_H__
