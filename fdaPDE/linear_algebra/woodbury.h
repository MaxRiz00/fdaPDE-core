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

#ifndef __WOODBURY_H__
#define __WOODBURY_H__

#include "../utils/symbols.h"

namespace fdapde {

// linear system solver based on the Sherman–Morrison–Woodbury formula
template <typename MatrixType>
    requires(is_eigen_dense_v<MatrixType>)
struct Woodbury {
    // solves linear system (A + U*C^{-1}*V)x = b, given already computed inversion of dense matrix C
    template <typename SparseSolver>
    MatrixType solve(
      const SparseSolver& invA, const MatrixType& U, const MatrixType& invC, const MatrixType& V, const MatrixType& b) {
        MatrixType y = invA.solve(b);   // y = A^{-1}b
        MatrixType Y = invA.solve(U);   // Y = A^{-1}U. Heavy step of the method. solver is more efficient for small q
        MatrixType G = invC + V * Y;    // compute dense matrix G = C^{-1} + V*A^{-1}*U = C^{-1} + V*y
        Eigen::PartialPivLU<MatrixType> invG;
        invG.compute(G);   // factorize qxq dense matrix G
        MatrixType t = invG.solve(V * y);
        // v = A^{-1}*U*t = A^{-1}*U*(C^{-1} + V*A^{-1}*U)^{-1}*V*A^{-1}*b by solving linear system A*v = U*t
        MatrixType v = invA.solve(U * t);
        return y - v;   // return system solution
    }
};

}   // namespace fdapde

#endif   // __WOODBURY_H__
