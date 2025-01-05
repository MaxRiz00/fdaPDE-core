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

#ifndef __FDAPDE_OPTIMIZATION_MODULE_H__
#define __FDAPDE_OPTIMIZATION_MODULE_H__

// clang-format off

// include required modules
#include "linear_algebra.h"    // pull-in Eigen first
#include "utility.h"
#include "fields.h"

// callbacks
#include "optimization/callbacks.h"
#include "optimization/backtracking_line_search.h"
#include "optimization/wolfe_line_search.h"
// algorithms
#include "optimization/grid.h"
#include "optimization/newton.h"
#include "optimization/gradient_descent.h"
#include "optimization/bfgs.h"

// clang-format on

#endif   // __FDAPDE_OPTIMIZATION_MODULE_H__
