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

#ifndef __FDAPDE_SPLINES_MODULE_H__
#define __FDAPDE_SPLINES_MODULE_H__

// include required modules
#ifndef __FDAPDE_FIELDS_MODULE_H__

#    include "fields.h"

#endif   // __FDAPDE_FIELDS_MODULE_H__

#ifndef __FDAPDE_GEOMETRY_MODULE_H__

#    include "geometry.h"

#endif   // __FDAPDE_GEOMETRY_MODULE_H__

#include "splines/bs_assembler_base.h"
#include "splines/bs_bilinear_form_assembler.h"
#include "splines/bs_linear_form_assembler.h"
#include "splines/bs_integration.h"
#include "splines/bs_objects.h"
#include "splines/bs_space.h"
#include "splines/dof_handler.h"
#include "splines/spline_basis.h"

#endif   // __FDAPDE_SPLINES_MODULE_H__
