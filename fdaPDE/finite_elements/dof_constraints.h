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

#ifndef __FDAPDE_DOF_CONSTRAINTS_H__
#define __FDAPDE_DOF_CONSTRAINTS_H__

#include "header_check.h"

namespace fdapde {

// managment of affine constraints on degrees of freedom of type \sum_{j} c_ij * dof_ij = b_j
template <typename DofHandler> class DofConstraints {
   public:
    using DofHandlerType = std::decay_t<DofHandler>;
    static constexpr int local_dim = DofHandlerType::local_dim;
    static constexpr int embed_dim = DofHandlerType::embed_dim;
    static constexpr double eps = 1e30;

    DofConstraints() = default;
    DofConstraints(const DofHandlerType& dof_handler) : dof_handler_(&dof_handler) { }

    // guarantees that the linear system Ax = b is such that all (affine) constraints are respected
    template <typename T> void enforce_constraints(T&& t) const {
        if constexpr (fdapde::is_subscriptable<std::decay_t<T>, int>) {   // linear system dense rhs
            for (const fdapde::Duplet<double>& duplet : constraint_values_) { t[duplet.row()] = duplet.value() * eps; }
        } else {   // linear system sparse matrix
            for (const fdapde::Triplet<double>& triplet : constraint_pattern_) {
                t.coeffRef(triplet.row(), triplet.col()) = triplet.value() * eps;
            }
        }
        return;
    }
    template <typename SystemMatrix, typename SystemRhs>
    void enforce_constraints(SystemMatrix&& A, SystemRhs&& b) const {
        fdapde_assert(A.rows() == b.rows());
	enforce_constraints(A);
	enforce_constraints(b);
        return;
    }
    // set dirichlet constraint type on boundary nodes with marker_id = marker
    template <typename... Callable> void set_dirichlet_constraint(int marker, Callable&&... g) {
        int n_boundary_dofs = dof_handler_->n_boundary_dofs(marker);
        fdapde_assert(
          sizeof...(Callable) == dof_handler_->dof_multiplicity() && (marker == BoundaryAll || n_boundary_dofs > 0));

        for (typename DofHandlerType::boundary_dofs_iterator it = dof_handler_->boundary_dofs_begin(marker);
             it != dof_handler_->boundary_dofs_end(marker); ++it) {
            int dof_id = it->id();
            if (dof_id < dof_handler_->n_unique_dofs()) {
                int component_id = 0;   // for vector elements, the component to which this dof refers to
                (
                  [&]() {
                      int dof_id_ = dof_id + component_id * dof_handler_->n_unique_dofs();
                      constraint_pattern_.emplace_back(dof_id_, dof_id_, 1.0);
                      constraint_values_ .emplace_back(dof_id_, g(it->coord()));
                      component_id++;
                  }(),
                  ...);
            }
        }
        return;
    }
   private:
    const DofHandlerType* dof_handler_;
    std::vector<fdapde::Triplet<double>> constraint_pattern_;   // triplets (i, j, c_ij)
    std::vector<fdapde::Duplet<double>>  constraint_values_;    // pairs (i, b_i)
};

}   // namespace fdapde

#endif   // __FDAPDE_DOF_CONSTRAINTS_H__
