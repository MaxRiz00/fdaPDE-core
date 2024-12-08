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

#ifndef __BS_BILINEAR_FORM_ASSEMBLER_H__
#define __BS_BILINEAR_FORM_ASSEMBLER_H__

#include <unordered_map>

#include "bs_assembler_base.h"

namespace fdapde {
namespace internals {
  
// galerkin and petrov-galerkin spline assembly loop
template <typename Triangulation_, typename Form_, int Options_, typename... Quadrature_>
class bs_bilinear_form_assembly_loop :
    public bs_assembler_base<Triangulation_, Form_, Options_, Quadrature_...>,
    public assembly_xpr_base<bs_bilinear_form_assembly_loop<Triangulation_, Form_, Options_, Quadrature_...>> { };

}   // namespace internals
}   // namespace fdapde

#endif   // __BS_BILINEAR_FORM_ASSEMBLER_H__
