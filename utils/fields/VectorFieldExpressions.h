#ifndef __VECTOR_FIELD_EXPRESSIONS_H__
#define __VECTOR_FIELD_EXPRESSIONS_H__

#include "../Symbols.h"
#include <functional>

namespace fdaPDE{
namespace core{

  // forward declarations
  template <typename T1, typename T2> class DotProduct;
  template <unsigned int M, unsigned int N> class VectConst;
  
  // Base class for any VectorField type
  struct VectBase {};
  
#define DEF_VECT_EXPR_OPERATOR(OPERATOR, FUNCTOR)			\
  template <int M, int N, typename E1, typename E2>			\
  VectBinOp<M,N, E1, E2, FUNCTOR> OPERATOR(				\
      const VectExpr<M,N, E1> &op1, const VectExpr<M,N, E2> &op2) {     \
  return VectBinOp<M,N, E1, E2, FUNCTOR>(		                \
	op1.get(), op2.get(), FUNCTOR());                               \
  }									\
  									\
  template <int M, int N, typename E>					\
  VectBinOp<M,N, VectConst<M,N>, E, FUNCTOR> OPERATOR(			\
      SVector<N> op1, const VectExpr<M, N, E> &op2) {                   \
  return VectBinOp<M,N, VectConst<M,N>, E, FUNCTOR>(                    \
	VectConst<M,N>(op1), op2.get(), FUNCTOR());                     \
  }									\
  									\
  template <int M, int N, typename E>					\
  VectBinOp<M,N, E, VectConst<M,N>, FUNCTOR> OPERATOR(			\
      const VectExpr<M,N, E> &op1, SVector<N> op2) {                    \
  return VectBinOp<M,N, E, VectConst<M,N>, FUNCTOR>(                    \
        op1.get(), VectConst<M,N>(op2), FUNCTOR());                     \
  }                                                                     \
  
  // Base class for vectorial expressions
  // M dimension of the space where the field is defined, N dimension of the arriving space
  template <int M, int N, typename E> struct VectExpr : public VectBase {
    // call operator[] on the base type E
    auto operator[](std::size_t i) const {
      return static_cast<const E&>(*this)[i];
    }
    // get underyling type composing the expression node
    const E& get() const { return static_cast<const E&>(*this); }
    // evaluate the expression at point p
    SVector<N> operator()(const SVector<M>& p) const {
      SVector<N> result;
      for(size_t i = 0; i < N; ++i){
	// trigger evaluation, call subscript of the underyling type. This will produce along the dimension i
	// a callable object, evaluate this passing the point p to get a double
	result[i] = operator[](i)(p);
      }
      return result;
    }
    // dot product between VectExpr and SVector
    virtual DotProduct<E, VectConst<M,N>> dot(const SVector<N>& op) const;
    // VectExpr - VectExpr dot product
    template <typename F>
    DotProduct<E,F> dot(const VectExpr<M,N,F>& op) const;
    // evaluate parametric nodes in the expression, does nothing if not redefined in derived classes
    template <typename T> void eval_parameters(T i) const { return; }
    // expose compile time informations
    static constexpr int rows = N;
    static constexpr int cols = 1;
    static constexpr int base = M; // dimensionality of base space
  };

  // an expression node representing a constant vector
  template <unsigned int M, unsigned int N>
  class VectConst : public VectExpr<M,N, VectConst<M,N>> {
  private:
    SVector<N> value_;
  public:
    VectConst(SVector<N> value) : value_(value) { }
    // return the stored value along direction i
    double operator[](std::size_t i) const { return value_[i]; }
  };

  // a parameter node
  template <unsigned int M, unsigned int N, typename F, typename T>
  class VectParam : public VectExpr<M,N, VectParam<M,N,F,T>> {
    // check F is callable with type T and returns an SVector<N>
    static_assert(std::is_same<decltype(std::declval<F>().operator()(T())), SVector<N>>::value);
  private:
    const F& f_;
    SVector<N> value_;
  public:
    // default constructor
    VectParam() = default;
    VectParam(const F& f) : f_(f) {};
    double operator[](std::size_t i) const { return value_[i]; }
    // evaluating the parameter makes a parametric node to act as a VectConst
    void eval_parameters(T i) { value_ = f_(i); }
  };
  
  // a generic binary operation node
  template <int M, int N, typename OP1, typename OP2, typename BinaryOperation>
  class VectBinOp : public VectExpr<M,N, VectBinOp<M,N, OP1, OP2, BinaryOperation>> {
  private:
    typename std::remove_reference<OP1>::type op1_;   // first  operand
    typename std::remove_reference<OP2>::type op2_;   // second operand
    BinaryOperation f_;                               // operation to apply
  public:
    // constructor
    VectBinOp(const OP1& op1, const OP2& op2, BinaryOperation f) : op1_(op1), op2_(op2), f_(f) { };
    // subscript operator. Let compiler to infer the return type (generally a FieldExpr)
    auto operator[](std::size_t i) const{
      return f_(op1_[i], op2_[i]);
    }
    // call parameter evaluation on operands
    template <typename T> const VectBinOp<M, N, OP1, OP2, BinaryOperation>& eval_parameters(T i) {
      op1_.eval_parameters(i); op2_.eval_parameters(i);
      return *this;
    }
  };  
  DEF_VECT_EXPR_OPERATOR(operator+, std::plus<> )
  DEF_VECT_EXPR_OPERATOR(operator-, std::minus<>)

  // support for double*VectExpr: multiplies each element of VectExpr by the scalar
  template <unsigned int M, unsigned int N> // represent a single scalar node in a vectorial expression.
  class VectScalar : public VectExpr<M,N, VectScalar<M,N>>{
  private:
    double value_;
  public:
    VectScalar(double value) : value_(value) { }
    double operator[](size_t i) const { return value_; }
  };
  template <int M, int N, typename E>
  VectBinOp<M,N, VectScalar<M,N>, E, std::multiplies<> >
  operator*(double op1, const VectExpr<M,N, E> &op2) {
    return VectBinOp<M,N, VectScalar<M,N>, E, std::multiplies<> >
      (VectScalar<M,N>(op1), op2.get(), std::multiplies<>());
  }  
  template <int M, int N, typename E>
  VectBinOp<M,N, E, VectScalar<M,N>, std::multiplies<> >
  operator*(const VectExpr<M,N, E> &op1, double op2) {
    return VectBinOp<M,N, E, VectScalar<M,N>, std::multiplies<> >
      (op1.get(), VectScalar<M,N>(op2), std::multiplies<>());
  }

  // dot product between a VectExpr and an (eigen) SVector.
  template <int M, int N, typename E>
  DotProduct<E, VectConst<M,N>> VectExpr<M,N, E>::dot(const SVector<N>& op) const {
    return DotProduct<E, VectConst<M,N>>(this->get(), VectConst<M,N>(op));
  }
  // dot product between a VectExpr and a VectExpr
  template <int M, int N, typename E>
  template <typename F>
  DotProduct<E,F> VectExpr<M,N, E>::dot(const VectExpr<M,N,F>& op) const {
    return DotProduct<E,F>(this->get(), op.get());
  }

  // allows for changes in the operands of an expression while keeping the same expression tree structure
  template <typename E> class VectPtr : public VectExpr<E::base, E::rows, VectPtr<E>> {
    static_assert(std::is_base_of<VectBase, E>::value);
  private:
    typename std::remove_reference<E>::type* ptr_;
  public:
    VectPtr(E* ptr) : ptr_(ptr) {};
    // delegate to pointed memory location
    auto operator[](std::size_t i) const{
      return ptr_->operator[](i);
    }
    // delegate to pointed memory location
    template <typename T> void eval_parameters(T i) {
      ptr_->eval_parameters(i);
      return;
    }
  };

}}

#endif // __VECTOR_FIELD_EXPRESSIONS_H__
