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

#ifndef __FDAPDE_SP_OBJECTS_H__
#define __FDAPDE_SP_OBJECTS_H__

#include "header_check.h"

namespace fdapde {  
namespace internals {

template <typename IsoSpace_>
struct iso_scalar_test_function_impl : public ScalarFieldBase<IsoSpace_::local_dim, TestFunction<IsoMesh_, iso_tag>> {
    using TestSpace = std::decay_t<IsoSpace_>;
    using Base = ScalarFieldBase<IsoSpace_::local_dim, TestFunction<IsoMesh_, iso_tag>>;
    using InputType = internals::iso_assembler_packet<TestSpace::local_dim>;
    using Scalar = double; 
    static constexpr int StaticInputSize = TestSpace::local_dim;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_values);

    private:
     template<typename Derived_>
     struct FirstPartialDerivative_ : ScalarFieldBase<TestSpace::local_dim, FirstPartialDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = FirstPartialDerivative_<T>;
        using TestSpace = std::decay_t<IsoSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TestSpace::local_dim, FirstDerivative_<Derived_>>;
        using InputType = internals::iso_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_grad);

        FirstPartialDerivative_() noexcept = default;
        FirstPartialDerivative_(const Derived_& xpr) noexcept : xpr_(xpr), i_(0) { }
        FirstPartialDerivative_(const Derived_& xpr, int i) noexcept : xpr_(xpr), i_(i) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& iso_packet) const { return iso_packet.test_grad(i_); }
        constexpr TestSpace& function_space() { return *(xpr_.iso_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.iso_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        int i_;
        typename internals::ref_select<const Derived>::type xpr_;
     };

     template <typename Derived_>
    struct MixedPartialDerivative_ : ScalarFieldBase<TestSpace::local_dim, MixedPartialDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = MixedPartialDerivative_<T>;      
        using TestSpace = std::decay_t<IsoSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TestSpace::local_dim, MixedPartialDerivative_<Derived_>>;
        using InputType = internals::sp_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_hessian);

        MixedPartialDerivative_() noexcept = default;
        MixedPartialDerivative_(const Derived_& xpr) noexcept : xpr_(xpr), i_(0), j_(0) { }
        MixedPartialDerivative_(const Derived_& xpr, int i, int j) noexcept : xpr_(xpr), i_(i), j_(j) { }

        // assembly evaluation
        constexpr Scalar operator()(const InputType& iso_packet) const { return iso_packet.test_hessian(i_,j_); }
        constexpr TestSpace& function_space() { return *(xpr_.iso_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.iso_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        int i_, j_;
        typename internals::ref_select<const Derived>::type xpr_;
    };

    // expose derivative types
    using FirstPartialDerivative  = FirstPartialDerivative_ <TestFunction<IsoSpace_, iso_tag>>;
    using MixedPartialDerivative = MixedPartialDerivative_<TestFunction<IsoSpace_, iso_tag>>;

    constexpr sp_scalar_test_function_impl() noexcept = default;
    constexpr sp_scalar_test_function_impl(IsoSpace_& iso_space) noexcept : iso_space_(std::addressof(iso_space)) { }  
    // assembly evaluation
    constexpr Scalar operator()(const InputType& iso_packet) const { return iso_packet.test_value; }
    constexpr TestSpace& function_space() { return *iso_space_; }
    constexpr const TestSpace& function_space() const { return *iso_space_; }
    constexpr int input_size() const { return StaticInputSize; }
   private:
    TestSpace* sp_space_;


};

template <typename IsoSpace_>
struct iso_scalar_test_function_impl : public ScalarFieldBase<IsoSpace_::local_dim, TestFunction<IsoMesh_, iso_tag>> {
    
}

} // namespace internals

template<typename IsoSpace_>
    requires(std::is_same_v<typename std::decay_t<IsoSpace_>::discretization_category, iso_tag>)
struct TestFunction<IsoSpace_, iso_tag> : public internals::iso_scalar_test_function_impl<IsoSpace_> {
    using Base = internals::iso_scalar_test_function_impl<IsoSpace_>;
    constexpr TestFunction() = default;
    constexpr TestFunction(IsoSpace_& iso_space) : Base(iso_space) { }
};


// partial derivatives of scalar test function
template <typename IsoSpace_>
struct PartialDerivative<TestFunction<IsoSpace_, spline_tag>, 1> :
    public TestFunction<IsoSpace_, iso_tag>::FirstPartialDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<IsoSpace_, spline_tag>& f, int i) :
        TestFunction<SpSpace_, spline_tag>::FirstPartialDerivative(f,i) { }
};
template <typename IsoSpace_>
struct PartialDerivative<TestFunction<IsoSpace_, spline_tag>, 2> :
    public TestFunction<IsoSpace_, iso_tag>::MixedPartialDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<IsoSpace_, spline_tag>& f, int i, int j) :
        TestFunction<SpSpace_, spline_tag>::MixedPartialDerivative(f,i,j) { }
};





} // namespace fdapde