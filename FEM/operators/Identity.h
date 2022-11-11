#ifndef __IDENTITY_H__
#define __IDENTITY_H__

#include <type_traits>
#include "../../utils/fields/ScalarField.h"
using fdaPDE::core::ScalarField;
#include "../../MESH/Element.h"
using fdaPDE::core::MESH::Element;
#include "BilinearFormExpressions.h"
using fdaPDE::core::FEM::BilinearFormExpr;

namespace fdaPDE{
namespace core{
namespace FEM{

  // class representing the identity operator (reaction term)
  template <typename T = DefaultOperator>
  class Identity : public BilinearFormExpr<Identity<T>>{
    // perform compile-time sanity checks
    static_assert(std::is_same<DefaultOperator, T>::value || // implicitly set c_ = 1
		  std::is_base_of<ScalarBase, T>::value ||   // space-varying case
		  std::is_floating_point<T>::value);         // constant coefficient case
  private:
    T c_; // reaction term
  public:
    // constructors
    Identity() = default;
    Identity(const T& c) : c_(c) {};

    std::tuple<Identity<T>> getTypeList() const { return std::make_tuple(*this); }
    static constexpr bool is_space_varying = std::is_base_of<ScalarBase, T>::value;

    // approximates the contribution to the (i,j)-th element of the discretization matrix given by the reaction term:

    // NOTE: is important to use auto return type to let the compiler return the whole expression template produced by this
    // operator avoiding both type erause (e.g. by casting to some ScalarField object) as well as the creation of temporaries
    template <typename... Args>
    auto integrate(const std::tuple<Args...>& mem_buffer) const {
      IMPORT_MEM_BUFFER_SYMBOLS(mem_buffer);
      if constexpr(std::is_same<DefaultOperator, T>::value)
	return psi_i*psi_j;
      else
	return c_*psi_i*psi_j;
    }
  };
  // template argument deduction guide
  template <typename T> Identity(const T&) -> Identity<T>;
  
}}}
#endif // __IDENTITY_H__
