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

#ifndef __BS_FUNCTION_H__
#define __BS_FUNCTION_H__

#include "../fields/scalar_field.h"
#include "bs_bilinear_form_assembler.h"

namespace fdapde {  
namespace internals {

template <typename BsSpace_>
struct bs_scalar_test_function_impl :
    public fdapde::ScalarBase<BsSpace_::local_dim, TestFunction<BsSpace_, bspline>> {
    fdapde_static_assert(BsSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using TestSpace = std::decay_t<BsSpace_>;
    using Base = fdapde::ScalarBase<BsSpace_::local_dim, TestFunction<BsSpace_, bspline>>;
    using InputType = internals::bs_assembler_packet<TestSpace::local_dim>;
    using Scalar = double;
    static constexpr int StaticInputSize = TestSpace::local_dim;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_values);
   private:
    template <typename Derived_>
    struct FirstDerivative_ : fdapde::ScalarBase<TestSpace::local_dim, FirstDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = FirstDerivative_<T>;
        using TestSpace = std::decay_t<BsSpace_>;   // required from xpr_query<>
        using Base = fdapde::ScalarBase<TestSpace::local_dim, FirstDerivative_<Derived_>>;
        using InputType = internals::bs_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_dx);

        FirstDerivative_() noexcept = default;
        FirstDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.test_dx; }
        constexpr TestSpace& function_space() { return *(xpr_.bs_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.bs_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
    template <typename Derived_>
    struct SecondDerivative_ : fdapde::ScalarBase<TestSpace::local_dim, SecondDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = SecondDerivative_<T>;      
        using TestSpace = std::decay_t<BsSpace_>;   // required from xpr_query<>
        using Base = fdapde::ScalarBase<TestSpace::local_dim, SecondDerivative_<Derived_>>;
        using InputType = internals::bs_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_ddx);

        SecondDerivative_() noexcept = default;
        SecondDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.test_ddx; }
        constexpr TestSpace& function_space() { return *(xpr_.bs_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.bs_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
   public:
    // expose derivative types
    using FirstDerivative  = FirstDerivative_ <TestFunction<BsSpace_, bspline>>;
    using SecondDerivative = SecondDerivative_<TestFunction<BsSpace_, bspline>>;

    constexpr bs_scalar_test_function_impl() noexcept = default;
    constexpr bs_scalar_test_function_impl(BsSpace_& bs_space) noexcept : bs_space_(std::addressof(bs_space)) { }  
    // assembly evaluation
    constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.test_value; }
    constexpr TestSpace& function_space() { return *bs_space_; }
    constexpr const TestSpace& function_space() const { return *bs_space_; }
    constexpr int input_size() const { return StaticInputSize; }
   private:
    TestSpace* bs_space_;
};

template <typename BsSpace_>
struct bs_scalar_trial_function_impl :
    public fdapde::ScalarBase<BsSpace_::local_dim, TrialFunction<BsSpace_, bspline>> {
    fdapde_static_assert(BsSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using TrialSpace = std::decay_t<BsSpace_>;
    using Base = fdapde::ScalarBase<BsSpace_::local_dim, TrialFunction<BsSpace_, bspline>>;
    using InputType = internals::bs_assembler_packet<TrialSpace::local_dim>;
    using Scalar = double;
    static constexpr int StaticInputSize = TrialSpace::local_dim;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_values);
   private:
    // definition of derivative functors
    template <typename Derived_>
    struct FirstDerivative_ : fdapde::ScalarBase<TrialSpace::local_dim, FirstDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = FirstDerivative_<T>;
        using TrialSpace = std::decay_t<BsSpace_>;   // required from xpr_query<>
        using Base = fdapde::ScalarBase<TrialSpace::local_dim, FirstDerivative_<Derived_>>;
        using InputType = internals::bs_assembler_packet<TrialSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TrialSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_dx);

        FirstDerivative_() noexcept = default;
        FirstDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.trial_dx; }
        constexpr TrialSpace& function_space() { return *(xpr_.bs_space_); }
        constexpr const TrialSpace& function_space() const { return *(xpr_.bs_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
    template <typename Derived_>
    struct SecondDerivative_ : fdapde::ScalarBase<TrialSpace::local_dim, SecondDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = SecondDerivative_<T>;
        using TrialSpace = std::decay_t<BsSpace_>;   // required from xpr_query<>
        using Base = fdapde::ScalarBase<TrialSpace::local_dim, SecondDerivative_<Derived_>>;
        using InputType = internals::bs_assembler_packet<TrialSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TrialSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(bs_assembler_flags::compute_shape_ddx);

        SecondDerivative_() noexcept = default;
        SecondDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.trial_ddx; }
        constexpr TrialSpace& function_space() { return *(xpr_.bs_space_); }
        constexpr const TrialSpace& function_space() const { return *(xpr_.bs_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
   public:
    // expose derivative types
    using FirstDerivative  = FirstDerivative_ <TrialFunction<BsSpace_, bspline>>;
    using SecondDerivative = SecondDerivative_<TrialFunction<BsSpace_, bspline>>;
  
    constexpr bs_scalar_trial_function_impl() noexcept = default;
    constexpr bs_scalar_trial_function_impl(BsSpace_& bs_space) noexcept : bs_space_(std::addressof(bs_space)) { }
    // assembly evaluation
    constexpr Scalar operator()(const InputType& bs_packet) const { return bs_packet.trial_value; }
    constexpr TrialSpace& function_space() { return *bs_space_; }
    constexpr const TrialSpace& function_space() const { return *bs_space_; }
    constexpr int input_size() const { return StaticInputSize; }
   private:
    TrialSpace* bs_space_;
};

}   // namespace internals

template <typename BsSpace_>
    requires(std::is_same_v<typename std::decay_t<BsSpace_>::space_category, bspline>)
struct TestFunction<BsSpace_, bspline> : public internals::bs_scalar_test_function_impl<BsSpace_> {
    using Base = internals::bs_scalar_test_function_impl<BsSpace_>;
    constexpr TestFunction() = default;
    constexpr TestFunction(BsSpace_& bs_space) : Base(bs_space) { }
};
// partial derivatives of scalar test function
template <typename BsSpace_>
struct PartialDerivative<TestFunction<BsSpace_, bspline>, 1> : public TestFunction<BsSpace_, bspline>::FirstDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<BsSpace_, bspline>& f, [[maybe_unused]] int i) :
        TestFunction<BsSpace_, bspline>::FirstDerivative(f) { }
};
template <typename BsSpace_>
struct PartialDerivative<TestFunction<BsSpace_, bspline>, 2> :
    public TestFunction<BsSpace_, bspline>::SecondDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<BsSpace_, bspline>& f, [[maybe_unused]] int i, [[maybe_unused]] int j) :
        TestFunction<BsSpace_, bspline>::SecondDerivative(f) { }
};
template <typename BsSpace_> constexpr auto dx (const TestFunction<BsSpace_, bspline>& xpr) {
    return typename TestFunction<BsSpace_, bspline>::FirstDerivative(xpr);
}
template <typename BsSpace_> constexpr auto ddx(const TestFunction<BsSpace_, bspline>& xpr) {
    return typename TestFunction<BsSpace_, bspline>::SecondDerivative(xpr);
}
  
template <typename BsSpace_>
    requires(std::is_same_v<typename std::decay_t<BsSpace_>::space_category, bspline>)
struct TrialFunction<BsSpace_, bspline> : public internals::bs_scalar_trial_function_impl<BsSpace_> {
    using Base = internals::bs_scalar_trial_function_impl<BsSpace_>;
    using TrialSpace = typename Base::TrialSpace;
    static constexpr int local_dim = BsSpace_::local_dim;
    static constexpr int embed_dim = BsSpace_::embed_dim;
    constexpr TrialFunction() = default;
    constexpr TrialFunction(BsSpace_& bs_space) : Base(bs_space) { }
    // norm evaluation
    double l2_squared_norm() {
        TrialFunction u(*Base::bs_space_);
        TestFunction  v(*Base::bs_space_);
        auto assembler = integrate(*Base::bs_space_->triangulation())(u * v);
        return coeff_.dot(assembler.assemble() * coeff_);
    }
    double l2_norm() { return std::sqrt(l2_squared_norm()); }
    const DVector<double>& coeff() const { return coeff_; }
    void set_coeff(const DVector<double>& coeff) { coeff_ = coeff; }
   private:
    Eigen::Matrix<double, Dynamic, 1> coeff_;
};
// partial derivative of scalar trial function
template <typename BsSpace_>
struct PartialDerivative<TrialFunction<BsSpace_, bspline>, 1> :
    public TrialFunction<BsSpace_, bspline>::FirstDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TrialFunction<BsSpace_, bspline>& f, [[maybe_unused]] int i) :
        TrialFunction<BsSpace_, bspline>::FirstDerivative(f, i) { }
};
template <typename BsSpace_>
struct PartialDerivative<TrialFunction<BsSpace_, bspline>, 2> :
    public TrialFunction<BsSpace_, bspline>::SecondDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TrialFunction<BsSpace_, bspline>& f, [[maybe_unused]] int i, [[maybe_unused]] int j) :
        TrialFunction<BsSpace_, bspline>::SecondDerivative(f, i, j) { }
};
template <typename BsSpace_> constexpr auto dx (const TrialFunction<BsSpace_, bspline>& xpr) {
    return typename TrialFunction<BsSpace_, bspline>::FirstDerivative(xpr);
}
template <typename BsSpace_> constexpr auto ddx(const TrialFunction<BsSpace_, bspline>& xpr) {
    return typename TrialFunction<BsSpace_, bspline>::SecondDerivative(xpr);
}

// representation of u(x) = \sum_{i=1}^{n_dofs} u_i \psi_i(x) with \{ \psi_i \}_i a B-Spline basis system
template <typename BsSpace_> class BsFunction : public fdapde::ScalarBase<BsSpace_::local_dim, BsFunction<BsSpace_>> {
    using Triangulation = typename BsSpace_::Triangulation;
    fdapde_static_assert(BsSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);  
   public:
    using BsSpace = std::decay_t<BsSpace_>;
    using Base = fdapde::ScalarBase<BsSpace_::local_dim, BsFunction<BsSpace_>>;
    using DofHandlerType = typename BsSpace::DofHandlerType;
    using InputType = Eigen::Matrix<double, BsSpace::local_dim, 1>;
    using Scalar = double;
    static constexpr int StaticInputSize = BsSpace::local_dim;
    static constexpr int Rows = 1;
    static constexpr int Cols = 1;
    static constexpr int NestAsRef = 1;
    static constexpr int local_dim = Triangulation::local_dim;
    static constexpr int embed_dim = Triangulation::embed_dim;
    static constexpr int XprBits = 0;

    BsFunction() = default;
    explicit BsFunction(BsSpace_& bs_space) : bs_space_(std::addressof(bs_space)) {
        coeff_ = Eigen::Matrix<double, Dynamic, 1>::Zero(bs_space_->n_dofs());
    }
    BsFunction(BsSpace_& bs_space, const Eigen::Matrix<double, Dynamic, 1>& coeff) :
        bs_space_(std::addressof(bs_space)), coeff_(coeff) {
        fdapde_assert(coeff.size() > 0 && coeff.size() == bs_space_->n_dofs());
    }
    Scalar operator()(const InputType& p) const { 
        int e_id = bs_space_->triangulation().locate(p);
        if (e_id == -1) return std::numeric_limits<Scalar>::quiet_NaN();   // return NaN if point lies outside domain
        // map p to reference interval [-1, 1]
        double a = bs_space_->triangulation().range()[0], b = bs_space_->triangulation().range()[1];
        double ref_p;
        if constexpr (fdapde::is_subscriptable<InputType, int>) { ref_p = p[0]; }
	else { ref_p = p; }
        ref_p = ((b - a) / 2) * ref_p + (b + a) / 2;
	// get active dofs
        typename DofHandlerType::CellType cell = bs_space_->dof_handler().cell(e_id);
	std::vector<int> active_dofs = cell.dofs();
        Scalar value = 0;
        for (int i = 0, n = active_dofs.size(); i < n; ++i) {
            value += coeff_[active_dofs[i]] * bs_space_->eval_shape_value(active_dofs[i], ref_p);
        }
        return value;
    }
    // norm evaluation
    double l2_squared_norm() {
        TrialFunction u(*bs_space_);
        TestFunction  v(*bs_space_);
        auto assembler = integrate(*(bs_space_->triangulation()))(u * v);
        return coeff_.dot(assembler.assemble() * coeff_);
    }
    double l2_norm() { return std::sqrt(l2_squared_norm()); }
    // getters
    const Eigen::Matrix<double, Dynamic, 1>& coeff() const { return coeff_; }
    constexpr BsSpace& function_space() { return *bs_space_; }
    constexpr const BsSpace& function_space() const { return *bs_space_; }
    constexpr int rows() const { return Rows; }
    constexpr int cols() const { return Cols; }
    constexpr int input_size() const { return StaticInputSize; }
    void set_coeff(const Eigen::Matrix<double, Dynamic, 1>& coeff) { coeff_ = coeff; }
    // linear algebra between fe functions
    friend constexpr BsFunction<BsSpace_> operator+(BsFunction<BsSpace_>& lhs, BsFunction<BsSpace_>& rhs) {
        return BsFunction<BsSpace_>(lhs.function_space(), lhs.coeff() + rhs.coeff());
    }
    friend constexpr BsFunction<BsSpace_> operator-(BsFunction<BsSpace_>& lhs, BsFunction<BsSpace_>& rhs) {
        return BsFunction<BsSpace_>(lhs.function_space(), lhs.coeff() - rhs.coeff());
    }
    // assignment from expansion coefficient vector
    BsFunction& operator=(const Eigen::Matrix<double, Dynamic, 1>& coeff) {
        fdapde_assert(coeff.size() > 0 && coeff.size() == bs_space_->n_dofs());
        coeff_ = coeff;
        return *this;
    }
   private:
    Eigen::Matrix<double, Dynamic, 1> coeff_;
    BsSpace* bs_space_;
};

}   // namespace fdapde

#endif   // __BS_FUNCTION_H__
