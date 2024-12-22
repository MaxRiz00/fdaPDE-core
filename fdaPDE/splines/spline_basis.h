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

#ifndef __SPLINE_BASIS_H__
#define __SPLINE_BASIS_H__

#include "../geometry/interval.h"
#include "../fields/spline.h"

namespace fdapde {

// given vector of knots u_1, u_2, ..., u_N, this class represents the set of N + order - 1 spline basis functions
// {l_1(x), l_2(x), ..., l_{N + order - 1}(x)} centered at knots u_1, u_2, ..., u_N
class SplineBasis {
   private:
    int order_;
    std::vector<Spline> basis_ {};
    std::vector<double> knots_ {};
   public:
    static constexpr int StaticInputSize = 1;
    static constexpr int Order = Dynamic;
    // constructors
    constexpr SplineBasis() : order_(0) { }
    // constructor from user defined knot vector
    template <typename KnotsVectorType>
        requires(requires(KnotsVectorType knots, int i) {
                    { knots[i] } -> std::convertible_to<double>;
                    { knots.size() } -> std::convertible_to<std::size_t>;
                })
    SplineBasis(KnotsVectorType&& knots, int order) : order_(order) {

        //fdapde_assert(std::is_sorted(knots.begin() FDAPDE_COMMA knots.end(), std::less_equal<double>()));
        int n = knots.size();
        knots_.resize(n);
        for (int i = 0; i < n; ++i) { knots_[i] = knots[i]; }
        // define basis system
        basis_.reserve(n - order_ + 1);
        for (int i = 0; i < n - order_ - 1; ++i) { basis_.emplace_back(knots_, i, order_); }
    }
    // constructor from geometric interval (no repeated knots)
    SplineBasis(const Triangulation<1, 1>& interval, int order, bool already_padded = false) : order_(order) {
        // construct knots vector
        Eigen::Matrix<double, Dynamic, 1> knots = interval.nodes();
        fdapde_assert(std::is_sorted(knots.begin() FDAPDE_COMMA knots.end() FDAPDE_COMMA std::less_equal<double>()));
        int n = knots.size();
        knots_.resize(n + 2 * order_);
        // pad the knot vector to obtain a full basis for the whole knot span [knots[0], knots[n-1]]
        if (not already_padded)
            for (int i = 0; i < n + 2 * order_; ++i) {
                if (i < order_) {
                    knots_[i] = knots[0];
                } else {
                    if (i < n + order_) {
                        knots_[i] = knots[i - order_];
                    } else {
                        knots_[i] = knots[n - 1];
                    }
                }
	}
        // define basis system
        basis_.reserve(knots_.size() - order_ + 1);
        for (int i = 0; i < knots_.size() - order_ - 1; ++i) { basis_.emplace_back(knots_, i, order_); }
    }

    //Algorithm A2. 2 from NURBS book pag. 70 computes all the nonvanishing
    //basis functions and stores them in the array N [0] , ... ,N [p]
    // it does not return zeros for the basis functions that are zero
    std::vector<double> compute_basis(double x) const {

        std::vector<double> N(order_ , 0.0);
        std::vector<double> left(order_, 0.0);
        std::vector<double> right(order_ , 0.0);
        N[0] = 1.0;

        // find the span of the knot vector
        int i = 0;
        while (x >= knots_[i+1] && i < knots_.size() - order_ - 1) i++;

        for (int j=1;j<=order_;j++){
            left[j] = x - knots_[i+1-j];
            right[j] = knots_[i+j] - x;
            double saved = 0.0;
            for (int r=0;r<j;r++){
                double temp = N[r]/(right[r+1]+left[j-r]);
                N[r] = saved + right[r+1]*temp;
                saved = left[j-r]*temp;
            }
            N[j] = saved;
        }

        // create a vector of the same size of the basis functions, copy N in the right position, zeros elsewhere
        std::vector<double> nonvanishing_basis(knots_.size() - order_ + 1, 0.0);
        for (int j = 0; j < order_ + 1; ++j) { nonvanishing_basis[i - order_ + j] = N[j]; }
        
        return nonvanishing_basis;
    }

    // getters
    constexpr const Spline& operator[](int i) const { return basis_[i]; }
    constexpr int size() const { return basis_.size(); }
    constexpr const std::vector<double>& knots_vector() const { return knots_; }
    int n_knots() const { return knots_.size(); }
    int order() const { return order_; }
};


} // namespace fdapde

#endif // __SPLINE_BASIS_H__
