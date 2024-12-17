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
  
// given a not fe_assembler_packet callable type Derived_, builds a map from a discrete set of points (e.g., quadrature
// nodes) to the evaluation of Derived_ at that points
template <typename Derived_> struct BsMap : public ScalarBase<Derived_::StaticInputSize, BsMap<Derived_>> {
   private:
    using OutputType = decltype(std::declval<Derived_>().operator()(std::declval<typename Derived_::InputType>()));
    using Derived = std::decay_t<Derived_>;
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
   public:
    using InputType = internals::bs_assembler_packet<Derived::StaticInputSize>;
    using Scalar = double;
    static constexpr int StaticInputSize = Derived::StaticInputSize;
    using Base = ScalarBase<StaticInputSize, BsMap<Derived>>;
    static constexpr int NestAsRef = 0;
    static constexpr int XprBits = Derived::XprBits | int(bs_assembler_flags::compute_physical_quad_nodes);
    static constexpr int ReadOnly = 1;
    static constexpr int Rows = 1;
    static constexpr int Cols = 1;

    constexpr BsMap() = default;
    constexpr BsMap(const Derived_& xpr) : xpr_(&xpr) { }
    template <typename CellIterator>
    void init(
      std::unordered_map<const void*, MatrixType>& buff, const MatrixType& nodes, [[maybe_unused]] CellIterator begin,
      [[maybe_unused]] CellIterator end) const {
        const void* ptr = reinterpret_cast<const void*>(xpr_);
        if (buff.find(ptr) == buff.end()) {
            DMatrix<double> mapped(nodes.rows(), Rows * Cols);
            for (int i = 0, n = nodes.rows(); i < n; ++i) { mapped(i, 0) = xpr_->operator()(nodes.row(i)); }
            buff[ptr] = mapped;
            map_ = &buff[ptr];
        } else {
            map_ = &buff[ptr];
        }
    }
    // bs assembler evaluation
    constexpr OutputType operator()(const InputType& bs_packet) const {
        return map_->operator()(bs_packet.quad_node_id, 0);
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

#endif   // __BS_MAP_H__
