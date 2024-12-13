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

#ifndef __BS_MAP_H__
#define __BS_MAP_H__

#include "../fields/scalar_field.h"
#include "../fields/meta.h"
#include "bs_function.h"

namespace fdapde {
  
// anytime you compose a trial or test function with a functor which is not callable at fe_assembler_packet, we wrap it
// into a BsMap, a fe_assembler_packet callable type encoding the functor evaluated at a fixed set of (quadrature) nodes
template <typename Derived_>
struct BsMap :
    public std::conditional_t<
      meta::is_scalar_field_v<Derived_>, ScalarBase<Derived_::StaticInputSize, BsMap<Derived_>>,
      MatrixBase<Derived_::StaticInputSize, BsMap<Derived_>>> {
   private:
    using OutputType = decltype(std::declval<Derived_>().operator()(std::declval<typename Derived_::InputType>()));
    static constexpr bool is_scalar = meta::is_scalar_field_v<Derived_>;
    using Derived = std::decay_t<Derived_>;
   public:
    using InputType = internals::bs_assembler_packet<Derived::StaticInputSize>;
    using Scalar = double;
    static constexpr int StaticInputSize = Derived::StaticInputSize;
    using Base = std::conditional_t<
      is_scalar, ScalarBase<StaticInputSize, BsMap<Derived>>, MatrixBase<StaticInputSize, BsMap<Derived>>>;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = Derived::XprBits | int(bs_assembler_flags::compute_physical_quad_nodes);
    static constexpr int ReadOnly = 1;
    static constexpr int Rows = []() { if constexpr(is_scalar) return 1; else return Derived::Rows; }();
    static constexpr int Cols = []() { if constexpr(is_scalar) return 1; else return Derived::Cols; }();

    constexpr BsMap() = default;
    constexpr BsMap(const Derived_& xpr) : xpr_(&xpr) { }
    template <typename CellIterator>
    void init(
      std::unordered_map<const void*, DMatrix<double>>& buff, const DMatrix<double>& nodes,
      [[maybe_unused]] CellIterator begin, [[maybe_unused]] CellIterator end) const {
        const void* ptr = reinterpret_cast<const void*>(xpr_);
        if (buff.find(ptr) == buff.end()) {
            DMatrix<double> mapped(nodes.rows(), Rows * Cols);
            if constexpr (is_scalar) {
                for (int i = 0, n = nodes.rows(); i < n; ++i) { mapped(i, 0)  = xpr_->operator()(nodes.row(i)); }
            } else {
                for (int i = 0, n = nodes.rows(); i < n; ++i) {
                    auto tmp = xpr_->operator()(nodes.row(i));
                    for (int j = 0; j < tmp.size(); ++j) { mapped(i, j) = tmp[j]; }
                }
            }
            buff[ptr] = mapped;
            map_ = &buff[ptr];
        } else {
            map_ = &buff[ptr];
        }
    }
    // fe assembler evaluation

  // here we need to better state what is the scalar field interface and what the matrix field one
  
    constexpr OutputType operator()(const InputType& fe_packet) const {
        if constexpr (is_scalar) {
            return map_->operator()(fe_packet.quad_node_id, 0);
        } else {
            return map_->row(fe_packet.quad_node_id);
        }
    }
    constexpr auto eval(int i, const InputType& fe_packet) const { return map_->operator()(fe_packet.quad_node_id, i); }
    constexpr const Derived& derived() const { return xpr_; }
    constexpr int input_size() const { return StaticInputSize; }
    constexpr int rows() const { return Rows; }
    constexpr int cols() const { return Cols; }
   private:
    const Derived* xpr_;
    mutable const DMatrix<Scalar>* map_;
};

}   // namespace fdapde

#endif   // __BS_MAP_H__
