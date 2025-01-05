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

#ifndef __FDAPDE_SYMBOLS_H__
#define __FDAPDE_SYMBOLS_H__

namespace fdapde {

  // could be removed, constants placed where they have sense...

constexpr int Dynamic = -1;   // used when the size of a vector or matrix is not known at compile time
constexpr int random_seed = -1;

// algorithm computation policies
[[maybe_unused]] static struct tag_exact { } Exact;
[[maybe_unused]] static struct tag_not_exact { } NotExact;

// matrix related symbols
[[maybe_unused]] constexpr int Upper = 0;       // lower triangular view of matrix
[[maybe_unused]] constexpr int Lower = 1;       // upper triangular view of matrix
[[maybe_unused]] constexpr int UnitUpper = 2;   // lower triangular view of matrix with ones on the diagonal
[[maybe_unused]] constexpr int UnitLower = 3;   // upper triangular view of matrix with ones on the diagonal
  
[[maybe_unused]] constexpr int RowMajor = 0;
[[maybe_unused]] constexpr int ColMajor = 1;


template <typename T> class Duplet {
   private:
    int row_;
    T value_;
   public:
    Duplet() = default;
    Duplet(int row, const T& value) : row_(row), value_(value) { }

    int row() const { return row_; }
    const T& value() const { return value_; }
    T& value() { return value_; }
};

// hash function for std::pair (allow pairs as key of unordered_map). inspired from boost::hash
struct pair_hash {
    template <typename T1, typename T2> std::size_t operator()(const std::pair<T1, T2>& pair) const {
        std::size_t hash = 0;
        hash ^= std::hash<T1>()(pair.first ) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        hash ^= std::hash<T2>()(pair.second) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
};
  
// hash function for an standard container's iterator range
template <typename T>
    requires(requires(T) {
        typename T::value_type;
        typename T::const_iterator;
    })
struct std_container_hash {
    std::size_t operator()(const typename T::const_iterator& begin, const typename T::const_iterator& end) const {
        std::size_t hash = 0;
        for (auto it = begin; it != end; ++it) {
            hash ^= std::hash<typename T::value_type>()(*it) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};
// hash function for std::array<T, N>
template <typename T, int N> struct std_array_hash {
    std::size_t operator()(const std::array<T, N>& array) const {
      return std_container_hash<std::array<T, N>>()(array.begin(), array.end());
    };
};

// reverse a std::unordered_map (if the map represents a bijection)
template <typename V, typename W> std::unordered_map<V, W> reverse(const std::unordered_map<W, V>& in) {
    std::unordered_map<V, W> out;
    for (const auto& [key, value] : in) {
        fdapde_assert(
          out.find(value) == out.end());   // if this is not passed, in is not a bijection, cannot be inverted
        out[value] = key;
    }
    return out;
}

}   // namespace fdapde

#endif   // __FDAPDE_SYMBOLS_H__
