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

#ifndef __FDAPDE_FIELDS_MODULE_H__
#define __FDAPDE_FIELDS_MODULE_H__

// clang-format off

// include required modules
#include "utility.h"

// scalar fields logic, as the matrix field one will depend on it
#include "fields/scalar_field.h"
#include "fields/divergence.h"
#include "fields/dot.h"
#include "fields/laplacian.h"
#include "fields/norm.h"
#include "fields/space_time_field.h"
// matrix field logic
#include "fields/jacobian.h"
#include "fields/matrix_field.h"
#include "fields/gradient.h"
#include "fields/hessian.h"

#include "fields/polynomial.h"
#include "fields/spline.h"

// clang-format on

#endif   // __FDAPDE_FIELDS_MODULE_H__
