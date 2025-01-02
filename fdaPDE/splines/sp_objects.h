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

#ifndef __SP_OBJECTS_H__
#define __SP_OBJECTS_H__

#include "../fields/scalar_field.h"
#include "sp_bilinear_form_assembler.h"

namespace fdapde {  
namespace internals {

template <typename SpSpace_>
struct sp_scalar_test_function_impl : public ScalarFieldBase<SpSpace_::local_dim, TestFunction<SpSpace_, bspline>> {
    fdapde_static_assert(SpSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using TestSpace = std::decay_t<SpSpace_>;
    using Base = ScalarFieldBase<SpSpace_::local_dim, TestFunction<SpSpace_, bspline>>;
    using InputType = internals::sp_assembler_packet<TestSpace::local_dim>;
    using Scalar = double;
    static constexpr int StaticInputSize = TestSpace::local_dim;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_values);
   private:
    template <typename Derived_>
    struct FirstDerivative_ : ScalarFieldBase<TestSpace::local_dim, FirstDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = FirstDerivative_<T>;
        using TestSpace = std::decay_t<SpSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TestSpace::local_dim, FirstDerivative_<Derived_>>;
        using InputType = internals::sp_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_dx);

        FirstDerivative_() noexcept = default;
        FirstDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.test_dx; }
        constexpr TestSpace& function_space() { return *(xpr_.sp_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.sp_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
    template <typename Derived_>
    struct SecondDerivative_ : ScalarFieldBase<TestSpace::local_dim, SecondDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = SecondDerivative_<T>;      
        using TestSpace = std::decay_t<SpSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TestSpace::local_dim, SecondDerivative_<Derived_>>;
        using InputType = internals::sp_assembler_packet<TestSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TestSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_ddx);

        SecondDerivative_() noexcept = default;
        SecondDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.test_ddx; }
        constexpr TestSpace& function_space() { return *(xpr_.sp_space_); }
        constexpr const TestSpace& function_space() const { return *(xpr_.sp_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
   public:
    // expose derivative types
    using FirstDerivative  = FirstDerivative_ <TestFunction<SpSpace_, bspline>>;
    using SecondDerivative = SecondDerivative_<TestFunction<SpSpace_, bspline>>;

    constexpr sp_scalar_test_function_impl() noexcept = default;
    constexpr sp_scalar_test_function_impl(SpSpace_& sp_space) noexcept : sp_space_(std::addressof(sp_space)) { }  
    // assembly evaluation
    constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.test_value; }
    constexpr TestSpace& function_space() { return *sp_space_; }
    constexpr const TestSpace& function_space() const { return *sp_space_; }
    constexpr int input_size() const { return StaticInputSize; }
   private:
    TestSpace* sp_space_;
};

template <typename SpSpace_>
struct sp_scalar_trial_function_impl : public ScalarFieldBase<SpSpace_::local_dim, TrialFunction<SpSpace_, bspline>> {
    fdapde_static_assert(SpSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using TrialSpace = std::decay_t<SpSpace_>;
    using Base = ScalarFieldBase<SpSpace_::local_dim, TrialFunction<SpSpace_, bspline>>;
    using InputType = internals::sp_assembler_packet<TrialSpace::local_dim>;
    using Scalar = double;
    static constexpr int StaticInputSize = TrialSpace::local_dim;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_values);
   private:
    // definition of derivative functors
    template <typename Derived_>
    struct FirstDerivative_ : ScalarFieldBase<TrialSpace::local_dim, FirstDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = FirstDerivative_<T>;
        using TrialSpace = std::decay_t<SpSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TrialSpace::local_dim, FirstDerivative_<Derived_>>;
        using InputType = internals::sp_assembler_packet<TrialSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TrialSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_dx);

        FirstDerivative_() noexcept = default;
        FirstDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.trial_dx; }
        constexpr TrialSpace& function_space() { return *(xpr_.sp_space_); }
        constexpr const TrialSpace& function_space() const { return *(xpr_.sp_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
    template <typename Derived_>
    struct SecondDerivative_ : ScalarFieldBase<TrialSpace::local_dim, SecondDerivative_<Derived_>> {
        using Derived = Derived_;
        template <typename T> using Meta = SecondDerivative_<T>;
        using TrialSpace = std::decay_t<SpSpace_>;   // required from xpr_query<>
        using Base = ScalarFieldBase<TrialSpace::local_dim, SecondDerivative_<Derived_>>;
        using InputType = internals::sp_assembler_packet<TrialSpace::local_dim>;
        using Scalar = double;
        static constexpr int StaticInputSize = TrialSpace::local_dim;
        static constexpr int NestAsRef = 0;
        static constexpr int XprBits = 0 | int(sp_assembler_flags::compute_shape_ddx);

        SecondDerivative_() noexcept = default;
        SecondDerivative_(const Derived_& xpr) noexcept : xpr_(xpr) { }
        // assembly evaluation
        constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.trial_ddx; }
        constexpr TrialSpace& function_space() { return *(xpr_.sp_space_); }
        constexpr const TrialSpace& function_space() const { return *(xpr_.sp_space_); }
        constexpr int input_size() const { return StaticInputSize; }
        constexpr const Derived& derived() const { return xpr_; }
       private:
        typename internals::ref_select<const Derived>::type xpr_;
    };
   public:
    // expose derivative types
    using FirstDerivative  = FirstDerivative_ <TrialFunction<SpSpace_, bspline>>;
    using SecondDerivative = SecondDerivative_<TrialFunction<SpSpace_, bspline>>;
  
    constexpr sp_scalar_trial_function_impl() noexcept = default;
    constexpr sp_scalar_trial_function_impl(SpSpace_& sp_space) noexcept : sp_space_(std::addressof(sp_space)) { }
    // assembly evaluation
    constexpr Scalar operator()(const InputType& sp_packet) const { return sp_packet.trial_value; }
    constexpr TrialSpace& function_space() { return *sp_space_; }
    constexpr const TrialSpace& function_space() const { return *sp_space_; }
    constexpr int input_size() const { return StaticInputSize; }
   private:
    TrialSpace* sp_space_;
};

}   // namespace internals

template <typename SpSpace_>
    requires(std::is_same_v<typename std::decay_t<SpSpace_>::space_category, bspline>)
struct TestFunction<SpSpace_, bspline> : public internals::sp_scalar_test_function_impl<SpSpace_> {
    using Base = internals::sp_scalar_test_function_impl<SpSpace_>;
    constexpr TestFunction() = default;
    constexpr TestFunction(SpSpace_& sp_space) : Base(sp_space) { }
};
// partial derivatives of scalar test function
template <typename SpSpace_>
struct PartialDerivative<TestFunction<SpSpace_, bspline>, 1> : public TestFunction<SpSpace_, bspline>::FirstDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<SpSpace_, bspline>& f, [[maybe_unused]] int i) :
        TestFunction<SpSpace_, bspline>::FirstDerivative(f) { }
};
template <typename SpSpace_>
struct PartialDerivative<TestFunction<SpSpace_, bspline>, 2> :
    public TestFunction<SpSpace_, bspline>::SecondDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TestFunction<SpSpace_, bspline>& f, [[maybe_unused]] int i, [[maybe_unused]] int j) :
        TestFunction<SpSpace_, bspline>::SecondDerivative(f) { }
};
template <typename SpSpace_> constexpr auto dx (const TestFunction<SpSpace_, bspline>& xpr) {
    return typename TestFunction<SpSpace_, bspline>::FirstDerivative(xpr);
}
template <typename SpSpace_> constexpr auto ddx(const TestFunction<SpSpace_, bspline>& xpr) {
    return typename TestFunction<SpSpace_, bspline>::SecondDerivative(xpr);
}
  
template <typename SpSpace_>
    requires(std::is_same_v<typename std::decay_t<SpSpace_>::space_category, bspline>)
struct TrialFunction<SpSpace_, bspline> : public internals::sp_scalar_trial_function_impl<SpSpace_> {
    using Base = internals::sp_scalar_trial_function_impl<SpSpace_>;
    using TrialSpace = typename Base::TrialSpace;
    static constexpr int local_dim = SpSpace_::local_dim;
    static constexpr int embed_dim = SpSpace_::embed_dim;
    constexpr TrialFunction() = default;
    constexpr TrialFunction(SpSpace_& sp_space) : Base(sp_space) { }
    // norm evaluation
    double l2_squared_norm() {
        TrialFunction u(*Base::sp_space_);
        TestFunction  v(*Base::sp_space_);
        auto assembler = integrate(*Base::sp_space_->triangulation())(u * v);
        return coeff_.dot(assembler.assemble() * coeff_);
    }
    double l2_norm() { return std::sqrt(l2_squared_norm()); }
    const Eigen::Matrix<double, Dynamic, 1>& coeff() const { return coeff_; }
    void set_coeff(const Eigen::Matrix<double, Dynamic, 1>& coeff) { coeff_ = coeff; }
   private:
    Eigen::Matrix<double, Dynamic, 1> coeff_;
};
// partial derivative of scalar trial function
template <typename SpSpace_>
struct PartialDerivative<TrialFunction<SpSpace_, bspline>, 1> :
    public TrialFunction<SpSpace_, bspline>::FirstDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TrialFunction<SpSpace_, bspline>& f, [[maybe_unused]] int i) :
        TrialFunction<SpSpace_, bspline>::FirstDerivative(f, i) { }
};
template <typename SpSpace_>
struct PartialDerivative<TrialFunction<SpSpace_, bspline>, 2> :
    public TrialFunction<SpSpace_, bspline>::SecondDerivative {
    PartialDerivative() = default;
    PartialDerivative(const TrialFunction<SpSpace_, bspline>& f, [[maybe_unused]] int i, [[maybe_unused]] int j) :
        TrialFunction<SpSpace_, bspline>::SecondDerivative(f, i, j) { }
};
template <typename SpSpace_> constexpr auto dx (const TrialFunction<SpSpace_, bspline>& xpr) {
    return typename TrialFunction<SpSpace_, bspline>::FirstDerivative(xpr);
}
template <typename SpSpace_> constexpr auto ddx(const TrialFunction<SpSpace_, bspline>& xpr) {
    return typename TrialFunction<SpSpace_, bspline>::SecondDerivative(xpr);
}

// representation of u(x) = \sum_{i=1}^{n_dofs} u_i \psi_i(x) with \{ \psi_i \}_i a B-Spline basis system
template <typename SpSpace_> class BsFunction : public ScalarFieldBase<SpSpace_::local_dim, BsFunction<SpSpace_>> {
    using Triangulation = typename SpSpace_::Triangulation;
    fdapde_static_assert(SpSpace_::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);  
   public:
    using BsSpace = std::decay_t<SpSpace_>;
    using Base = ScalarFieldBase<SpSpace_::local_dim, BsFunction<SpSpace_>>;
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
    explicit BsFunction(SpSpace_& sp_space) : sp_space_(std::addressof(sp_space)) {
        coeff_ = Eigen::Matrix<double, Dynamic, 1>::Zero(sp_space_->n_dofs());
    }
    BsFunction(SpSpace_& sp_space, const Eigen::Matrix<double, Dynamic, 1>& coeff) :
        sp_space_(std::addressof(sp_space)), coeff_(coeff) {
        fdapde_assert(coeff.size() > 0 && coeff.size() == sp_space_->n_dofs());
    }
    Scalar operator()(const InputType& p) const { 
        int e_id = sp_space_->triangulation().locate(p);
        if (e_id == -1) return std::numeric_limits<Scalar>::quiet_NaN();   // return NaN if point lies outside domain
        // map p to reference interval [-1, 1]
        double a = sp_space_->triangulation().range()[0], b = sp_space_->triangulation().range()[1];
        double ref_p;
        if constexpr (fdapde::is_subscriptable<InputType, int>) { ref_p = p[0]; }
	else { ref_p = p; }
        ref_p = (ref_p - (b - a) / 2) * 2 / (b + a);
        // get active dofs
        typename DofHandlerType::CellType cell = sp_space_->dof_handler().cell(e_id);
	std::vector<int> active_dofs = cell.dofs();
        Scalar value = 0;
        for (int i = 0, n = active_dofs.size(); i < n; ++i) {
            value += coeff_[active_dofs[i]] *
                     sp_space_->eval_shape_value(active_dofs[i], fdapde::Matrix<double, embed_dim, 1>(ref_p));
        }
        return value;
    }
    // norm evaluation
    double l2_squared_norm() {
        TrialFunction u(*sp_space_);
        TestFunction  v(*sp_space_);
        auto assembler = integrate(*(sp_space_->triangulation()))(u * v);
        return coeff_.dot(assembler.assemble() * coeff_);
    }
    double l2_norm() { return std::sqrt(l2_squared_norm()); }
    // getters
    const Eigen::Matrix<double, Dynamic, 1>& coeff() const { return coeff_; }
    constexpr BsSpace& function_space() { return *sp_space_; }
    constexpr const BsSpace& function_space() const { return *sp_space_; }
    constexpr int rows() const { return Rows; }
    constexpr int cols() const { return Cols; }
    constexpr int input_size() const { return StaticInputSize; }
    void set_coeff(const Eigen::Matrix<double, Dynamic, 1>& coeff) { coeff_ = coeff; }
    // linear algebra between fe functions
    friend constexpr BsFunction<SpSpace_> operator+(BsFunction<SpSpace_>& lhs, BsFunction<SpSpace_>& rhs) {
        return BsFunction<SpSpace_>(lhs.function_space(), lhs.coeff() + rhs.coeff());
    }
    friend constexpr BsFunction<SpSpace_> operator-(BsFunction<SpSpace_>& lhs, BsFunction<SpSpace_>& rhs) {
        return BsFunction<SpSpace_>(lhs.function_space(), lhs.coeff() - rhs.coeff());
    }
    // assignment from expansion coefficient vector
    BsFunction& operator=(const Eigen::Matrix<double, Dynamic, 1>& coeff) {
        fdapde_assert(coeff.size() > 0 && coeff.size() == sp_space_->n_dofs());
        coeff_ = coeff;
        return *this;
    }
   private:
    Eigen::Matrix<double, Dynamic, 1> coeff_;
    BsSpace* sp_space_;
};

// given a not sp_assembler_packet callable type Derived_, builds a map from a discrete set of points (e.g., quadrature
// nodes) to the evaluation of Derived_ at that points, so that the results is sp_assembler_packet evaluable
template <typename Derived_> struct SpMap : public ScalarBase<Derived_::StaticInputSize, SpMap<Derived_>> {
   private:
    using OutputType = decltype(std::declval<Derived_>().operator()(std::declval<typename Derived_::InputType>()));
    using Derived = std::decay_t<Derived_>;
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
   public:
    using InputType = internals::sp_assembler_packet<Derived::StaticInputSize>;
    using Scalar = double;
    static constexpr int StaticInputSize = Derived::StaticInputSize;
    using Base = ScalarBase<StaticInputSize, SpMap<Derived>>;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = Derived::XprBits | int(sp_assembler_flags::compute_physical_quad_nodes);
    static constexpr int ReadOnly = 1;
    static constexpr int Rows = 1;
    static constexpr int Cols = 1;

    constexpr SpMap() = default;
    constexpr SpMap(const Derived_& xpr) : xpr_(&xpr) { }
    template <typename CellIterator>
    void init(
      std::unordered_map<const void*, MatrixType>& buff, const MatrixType& nodes, [[maybe_unused]] CellIterator begin,
      [[maybe_unused]] CellIterator end) const {
        const void* ptr = reinterpret_cast<const void*>(xpr_);
        if (buff.find(ptr) == buff.end()) {
            Eigen::Matrix<double, Dynamic, Dynamic> mapped(nodes.rows(), Rows * Cols);
            for (int i = 0, n = nodes.rows(); i < n; ++i) { mapped(i, 0) = xpr_->operator()(nodes.row(i)); }
            buff[ptr] = mapped;
            map_ = &buff[ptr];
        } else {
            map_ = &buff[ptr];
        }
    }
    // bs assembler evaluation
    constexpr OutputType operator()(const InputType& sp_packet) const {
        return map_->operator()(sp_packet.quad_node_id, 0);
    }
    constexpr const Derived& derived() const { return xpr_; }
    constexpr int input_size() const { return StaticInputSize; }
    constexpr int rows() const { return Rows; }
    constexpr int cols() const { return Cols; }
   private:
    const Derived* xpr_;
    mutable const MatrixType* map_;
};
 
}   // namespace fdapde

#endif   // __SP_OBJECTS_H__
