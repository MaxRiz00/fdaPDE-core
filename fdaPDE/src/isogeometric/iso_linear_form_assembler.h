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

#ifndef __FDAPDE_ISO_LINEAR_FORM_ASSEMBLER_H__
#define __FDAPDE_ISO_LINEAR_FORM_ASSEMBLER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {
  
// assembly loop for the discretization of integrals \int_D \langle f, \psi_i \rangle, with \psi_i \in test space
template <typename Mesh_, typename Form_, int Options_, typename... Quadrature_>
class iso_linear_form_assembly_loop :
    public iso_assembler_base<Mesh_, Form_, Options_, Quadrature_...>,
    public assembly_xpr_base<iso_linear_form_assembly_loop<Mesh_, Form_, Options_, Quadrature_...>> {
   public:
    using Base = iso_assembler_base<Mesh_, Form_, Options_, Quadrature_...>;
    using Form = typename Base::Form;
    using discretization_category = typename Base::discretization_category;
    static constexpr int local_dim = Base::local_dim;
    static constexpr int embed_dim = Base::embed_dim;
    using Base::dof_handler_;
    using Base::form_;
  
    iso_linear_form_assembly_loop() = default;
    iso_linear_form_assembly_loop(
      const Form_& form, typename Base::geo_iterator begin, typename Base::geo_iterator end,
      const Quadrature_&... quadrature) requires(sizeof...(quadrature) <= 1)
        : Base(form, begin, end, quadrature...) { }

    Eigen::Matrix<double, Dynamic, 1> assemble() const {
        Eigen::Matrix<double, Dynamic, 1> assembled_vec(dof_handler_->n_dofs());
        assembled_vec.setZero();
        assemble(assembled_vec);
        return assembled_vec;
    }
    void assemble(Eigen::Matrix<double, Dynamic, 1>& assembled_vec) const {
        //std::cout << "Assembling linear form..." << std::endl;
        using iterator = typename Base::dof_iterator;
        iterator begin(Base::begin_.index(), dof_handler_, Base::begin_.marker());
        iterator end  (Base::end_.index(),   dof_handler_, Base::end_.marker()  );

        // prepare assembly loop
        std::vector<int> active_dofs;
        int n = 1;

        for(int i = 0; i<local_dim; i++){
            n*= (Base::test_space_->degree())[i] + 1;
        }

        int q = Base::n_quadrature_nodes_;
        //std::cout << "n: " << n << ", q: " << q << std::endl;
        MdArray<double, MdExtents<Dynamic, Dynamic>> shape_values(n, q);

        MdArray<double, MdExtents<Dynamic>> metric_dets(q);

        
        //std::unordered_map<const void*, Eigen::Matrix<double, Dynamic, Dynamic>> iso_map_buff;
        if constexpr (Form::XprBits & int(iso_assembler_flags::compute_physical_quad_nodes)) {
            Base::distribute_quadrature_nodes(
              begin, end);   // distribute quadrature nodes on physical mesh (if required)
        }
        
        // start assembly loop
        internals::iso_assembler_packet<embed_dim> iso_packet {};  
	    int local_cell_id = 0;
        //std::cout<<"Start linear assembly loop..."<<std::endl;
        for (iterator it = begin; it != end; ++it) {
            iso_packet.cell_measure = it->parametric_measure();
            active_dofs = it->dofs();
            //std::cout << std::endl;
            Base::eval_param_shape_values(Base::test_space_->basis(), active_dofs, it, shape_values);
            //if constexpr(embed_dim == 2){
            //    for(int q_k = 0; q_k < Base::n_quadrature_nodes_; ++q_k) {
            //        metric_dets(q_k) = 1;
            //    }
            //}
            //else
            Base::eval_metric_determinant(it, metric_dets); // metric det(F^T F)

            
                
            // perform integration of weak form for (i, j)-th basis pair
            for (int i = 0; i < n; ++i) {   // test function loop
                double value = 0;
                for (int q_k = 0; q_k < Base::n_quadrature_nodes_; ++q_k) {
                    iso_packet.test_value = shape_values(i, q_k);
                    if constexpr (Form::XprBits & int(iso_assembler_flags::compute_physical_quad_nodes)) {
                        iso_packet.quad_node_id = local_cell_id * Base::n_quadrature_nodes_ + q_k;
                    }
                    //std::cout<<"Quad node: "<<q_k<<std::endl;
                    value += Base::quad_weights_(q_k, 0) * form_(iso_packet) * metric_dets(q_k) ;
                }
                assembled_vec[active_dofs[i]] += value * iso_packet.cell_measure; // ????? forse il determinante metrico
            }
	    local_cell_id++;
        }
        return;
    }
    constexpr int n_dofs() const { return dof_handler_->n_dofs(); }
    constexpr int rows() const { return n_dofs(); }
    constexpr int cols() const { return 1; }
};
  
}   // namespace internals
}   // namespace fdapde

#endif   // __FDAPDE_FE_LINEAR_FORM_ASSEMBLER_H__
