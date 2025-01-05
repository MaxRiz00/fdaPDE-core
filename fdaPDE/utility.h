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

#ifndef __FDAPDE_UTILITY_MODULE_H__
#define __FDAPDE_UTILITY_MODULE_H__

// clang-format off

// STL includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <limits>
#include <memory>   // for std::shared_ptr
#include <type_traits>
#include <optional>
#include <random>
// common STL containers
#include <array>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// utils include
#include "utility/assert.h"
#include "utility/xpr_helper.h" // move this in fields, and place ref here
#include "utility/symbols.h"
#include "utility/traits.h"

#include "utility/matrix.h"
#include "utility/binary.h"
#include "utility/mdarray.h"
#include "utility/binary_tree.h"
#include "utility/type_erasure.h"
#include "utility/numeric.h"

// clang-format on

#endif   // __FDAPDE_UTILITY_MODULE_H__
