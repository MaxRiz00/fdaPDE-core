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

#ifndef __DIVERGENCE_H__
#define __DIVERGENCE_H__

#include "scalar_field.h"

namespace fdapde {

template <typename Derived_> class Divergence : public ScalarBase<Derived_::StaticInputSize, Divergence<Derived_>> {
    fdapde_static_assert(
      Derived_::Cols == 1 && (Derived_::StaticInputSize == Dynamic || Derived_::StaticInputSize == Derived_::Rows),
      DIVERGENCE_OPERATOR_IS_FOR_VECTOR_FIELDS_ONLY);  
   public:
    using Derived = Derived_;
   private:
    class wrap_t {   // wrapper for subscripted types having just a call operator
        using functor_t = std::decay_t<decltype(std::declval<Derived>().operator[](std::declval<int>()))>;
        functor_t xpr_;
       public:
        using InputType = typename Derived::InputType;
        using Scalar = typename Derived::Scalar;
        static constexpr int StaticInputSize = Derived::StaticInputSize;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0;

        constexpr wrap_t() noexcept : xpr_() { }
        constexpr wrap_t(const functor_t& xpr) noexcept : xpr_(xpr) { }
        // call operator
        template <typename InputType_> constexpr Scalar operator()(InputType_&& p) const { return xpr_(p); }
    };
   public:
    template <typename T> using Meta = Divergence<T>;
    using Base = ScalarBase<Derived::StaticInputSize, Divergence<Derived>>;
    // using FunctorType =
    //   PartialDerivative<std::decay_t<decltype(std::declval<Derived>().operator[](std::declval<int>()))>, 1>;
    using FunctorType = PartialDerivative<wrap_t, 1>;

    using InputType = typename Derived::InputType;
    using Scalar = typename Derived::Scalar;
    static constexpr int StaticInputSize = Derived::StaticInputSize;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = FunctorType::XprBits;

    explicit constexpr Divergence(const Derived_& xpr) : Base(), data_(), xpr_(xpr) {
        if constexpr (StaticInputSize == Dynamic) {
            fdapde_constexpr_assert(xpr_.input_size() == xpr_.rows());
            data_.resize(xpr_.rows());
        }
        for (int i = 0; i < xpr_.rows(); ++i) { data_[i] = FunctorType(xpr_[i], i); }
    }
    constexpr Scalar operator()(const InputType& p) const {
        Scalar div_ = 0;
        for (int i = 0; i < xpr_.rows(); ++i) { div_ += data_[i](p); }
        return div_;
    }
    constexpr int input_size() const { return xpr_.input_size(); }
    constexpr const Derived& derived() const { return xpr_; }
   private:
    using StorageType = typename std::conditional_t<
      Derived::StaticInputSize == Dynamic, std::vector<FunctorType>, std::array<FunctorType, StaticInputSize>>;
    StorageType data_;
    typename internals::ref_select<const Derived>::type xpr_;
};

template <typename XprType> Divergence<XprType> constexpr div(const XprType& xpr) { return Divergence<XprType>(xpr); }

}   // namespace fdapde

#endif // __DIVERGENCE_H__
