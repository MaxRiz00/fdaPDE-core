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

#ifndef __FDAPDE_ISO_BILINEAR_FORM_ASSEMBLER_H__
#define __FDAPDE_ISO_BILINEAR_FORM_ASSEMBLER_H__

#include "header_check.h"

namespace fdapde {
namespace internals {

template<typename IsoMesh_, typename Form_, int Options_, typename... Quadrature_>
class iso_bilinear_form_assembly_loop :
    public iso_assembler_base<IsoMesh_, Form_, Options_, Quadrature_...>,
    public assembly_xpr_base<iso_bilinear_form_assembly_loop<IsoMesh_, Form_, Options_, Quadrature_...>> {
    // detect trial and test spaces from bilinear form
    using TestSpace = test_space_t<Form_>;
    using TrialSpace = trial_space_t<Form_>;
    static_assert(TestSpace::local_dim == TrialSpace::local_dim && TestSpace::embed_dim == TrialSpace::embed_dim);
    static constexpr bool is_galerkin = std::is_same_v<TestSpace, TrialSpace>;
    static constexpr bool is_petrov_galerkin = !is_galerkin;
    using Base = iso_assembler_base<IsoMesh_, Form_, Options_, Quadrature_...>;
    using Form = typename Base::Form;
    
    using DofHandlerType = typename Base::DofHandlerType;
    using discretization_category = typename TestSpace::discretization_category;

    fdapde_static_assert(
        std::is_same_v<discretization_category FDAPDE_COMMA iso_tag>, TEST_AND_TRIAL_SPACE_MUST_HAVE_THE_SAME_DISCRETIZATION_CATEGORY);
    static constexpr int local_dim = Base::local_dim;
    static constexpr int embed_dim = Base::embed_dim;
    using Base::form_;
    using Base::test_space_;
    // private data members
    const DofHandlerType* trial_dof_handler_;
    constexpr const DofHandlerType* test_dof_handler() const { return Base::dof_handler_; }
    constexpr const DofHandlerType* trial_dof_handler() const {
        return is_galerkin ? Base::dof_handler_ : trial_dof_handler_;
    }
    const TrialSpace* trial_space_;

    public:

    inline int map_dof(int i, const std::vector<int>& map) {
        return map.empty() ? i : map[i];
    }
    
    iso_bilinear_form_assembly_loop() = default;
    
    iso_bilinear_form_assembly_loop(
        const Form_& form, typename Base::geo_iterator begin, typename Base::geo_iterator end, const Quadrature_&... quadrature)
        requires(sizeof...(quadrature)<=1) 
        : Base(form, begin, end, quadrature...), trial_space_(std::addressof(internals::trial_space(form_))){
            // print begin cell id
            //std::cout<<"ciaooo Begin cell id: "<<begin->id()<<std::endl;
            //std::cout<<"ciaoo End cell id: "<<end->id()<<std::endl;
            if constexpr(is_petrov_galerkin){
                trial_dof_handler_ = std::addressof(internals::trial_space(form_).dof_handler());
            }
            fdapde_assert(test_dof_handler()->n_dofs() != 0 && trial_dof_handler()->n_dofs() != 0);

            if constexpr (sizeof...(Quadrature_) == 0) {
                // default to higher-order quadrature
                /*
                if (test_space_->order() != trial_space_->order()) {
                    internals::get_iso_quadrature(
                    test_space_->order() > trial_space_->order() ? test_space_->order() : trial_space_->order(),
                    Base::quad_nodes_, Base::quad_weights_);
                }
                */
            }
    }


    Eigen::SparseMatrix<double> assemble() const {
        //std::cout << "MATRIX prima..." << std::endl;
        Eigen::SparseMatrix<double> assembled_mat(test_dof_handler()->n_dofs(), trial_dof_handler()->n_dofs());
        //std::cout << "MATRIX dopo..." << std::endl;
        std::vector<Eigen::Triplet<double>> triplet_list;
        //std::cout << "Assembling bilinear form..." << std::endl;
        assemble(triplet_list);
        //std::cout << "Assembling bilinear form done." << std::endl;
        // linearity of the integral is implicitly used here, as duplicated triplets are summed up (see Eigen docs)
        assembled_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
        assembled_mat.makeCompressed();
        return assembled_mat;
    }


    void assemble(std::vector<Eigen::Triplet<double>>& triplet_list) const{
        using iterator = typename Base::dof_iterator;
        iterator begin(Base::begin_.index(), test_dof_handler(), Base::begin_.marker());
        iterator end  (Base::end_.index()  , test_dof_handler(), Base::end_.marker());

        // prepare assembly loop
        std::vector<int> test_active_dofs, trial_active_dofs;
        int q = Base::n_quadrature_nodes_;
        int n1 = 1, n2 = 1;

        for(int i = 0; i<local_dim; i++){
            n1*= (test_space_->order())[i] + 1;
            n2*= (is_galerkin ? (test_space_->order())[i] : (trial_space_->order())[i]) + 1;
        }

        MdArray<double, MdExtents<Dynamic, Dynamic>> test_param_shape_values(n1,q), trial_param_shape_values(n2, q);
        MdArray<double, MdExtents<Dynamic, Dynamic, Dynamic>> 
            test_param_shape_grads(n1, q, local_dim), trial_param_shape_grads(n2, q, local_dim);
        MdArray<double, MdExtents<Dynamic, Dynamic, Dynamic, Dynamic>> 
            test_param_shape_hess(n1, q, local_dim, local_dim), trial_param_shape_hess(n2, q, local_dim, local_dim);

        MdArray<Eigen::Matrix<double, embed_dim, local_dim> , MdExtents< Dynamic>> 
            param_grad(q);

        MdArray<MdArray<double, MdExtents<embed_dim, local_dim, local_dim>>, MdExtents<Dynamic>> 
            param_hess(q);
        
        MdArray<double, MdExtents<Dynamic>> metric_dets(q);
        

        // distribute quadrature nodes on physical mesh (if required) ..... da capire
        //std::cout << "Distributing quadrature nodes..." << std::endl;

        // start assembly loop
        internals::iso_assembler_packet<embed_dim> iso_packet {};
        //std::cout << "Assembling cells..." << std::endl;
        int local_cell_id = 0;

        auto reduced_dof_map = test_dof_handler()->reduced_dof_map();
        for(iterator it = begin; it!= end; ++it) {
            //std::cout << "Assembling cell " << it->id() << std::endl;
            test_active_dofs = it->dofs();

            //std::cout << "Test active dofs: ";
            //for(auto dof : test_active_dofs) std::cout<<dof<<", ";
            //std::cout << std::endl;
            if constexpr (is_petrov_galerkin) { trial_active_dofs = trial_dof_handler()->active_dofs(it->id()); }
            // update the iso_packet
            iso_packet.cell_measure = it->parametric_measure();
            //std::cout << "Cell measure: " << iso_packet.cell_measure << std::endl;

            if constexpr (Form::XprBits & int(iso_assembler_flags::compute_shape_values)) {
                //std::cout<<"Computing values..."<<std::endl;
                Base::eval_param_shape_values(test_space_->basis(), test_active_dofs, it, test_param_shape_values);
                Base::eval_param_shape_values(
                    trial_space_->basis(), is_petrov_galerkin ? trial_active_dofs : test_active_dofs, it,
                    trial_param_shape_values);
            }
            if constexpr ( int(iso_assembler_flags::compute_shape_grad)) { //Form::XprBits &
                //std::cout<<"Computing grads..."<<std::endl;
                Base::eval_param_shape_grads(test_space_->basis(), test_active_dofs, it, test_param_shape_grads);
                Base::eval_param_shape_grads(
                  trial_space_->basis(), is_petrov_galerkin ? trial_active_dofs : test_active_dofs, it, trial_param_shape_grads);
                
            }
            if constexpr (Form::XprBits & int(iso_assembler_flags::compute_shape_hess)) { //Form::XprBits &
                //std::cout<<"Computing hessians..."<<std::endl;
                Base::eval_param_shape_hess(test_space_->basis(), test_active_dofs, it, test_param_shape_hess);
                Base::eval_param_shape_hess(
                  trial_space_->basis(), is_petrov_galerkin ? trial_active_dofs : test_active_dofs, it,
                  trial_param_shape_hess);
            }

            Base::eval_param_grad(it, param_grad); // F
            Base::eval_metric_determinant(it, metric_dets); // metric sqrt det(F^T F)
            Base::eval_param_hess(it, param_hess); 
            
            // precompute F(F^T F)^-1 for each q_k
            MdArray<Eigen::Matrix<double, embed_dim, local_dim> , MdExtents< Dynamic>> grad_transf(q);
            for (int q_k = 0; q_k < Base::n_quadrature_nodes_; ++q_k) {
                grad_transf(q_k) = param_grad(q_k) * (param_grad(q_k).transpose() * param_grad(q_k)).inverse();
            }
            //std::cout << "Grad transf: " << grad_transf(0).rows() << " x " << grad_transf(0).cols() << std::endl;
            // print Form::XprBits
            //std::cout << "Form::XprBits = 0x" << std::hex << Form::XprBits << std::endl;
            //std::cout << "compute_shape_hess = 0x" << std::hex << int(fdapde::iso_assembler_flags::compute_shape_hess) << std::endl;

            std::map<std::pair<int, int>, double> mat_entries;

            // perform integration of weak form for (i,j)-th basis pair
            for(int i = 0; i<n2; ++i){
                for(int j = 0; j<n1; ++j){
                    double value = 0;
                    //std::cout<<"i = "<<i<<std::endl;
                    //std::cout<<"j = "<<j<<std::endl;
                    for (int q_k = 0; q_k < Base::n_quadrature_nodes_; ++q_k) {
                        if constexpr (Form::XprBits & int(iso_assembler_flags::compute_shape_values)) {
                            iso_packet.trial_value = trial_param_shape_values(i, q_k) ;
                            iso_packet.test_value  = test_param_shape_values (j, q_k) ;
                        }
                        if constexpr (Form::XprBits & int(iso_assembler_flags::compute_shape_grad)) {
                            auto temp_trial_grad = trial_param_shape_grads.template slice<0,1>(i, q_k); 
                            auto temp_test_grad  = test_param_shape_grads.template slice<0,1>(j, q_k);

                            // Project onto the tangent plane
                            Eigen::Matrix<double, embed_dim, 1> n = ((param_grad(q_k).col(0)).cross(param_grad(q_k).col(1))).normalized();
                            Eigen::Matrix<double, embed_dim, embed_dim> P = Eigen::Matrix<double, embed_dim, embed_dim>::Identity() - n * n.transpose();

                            // dxi/dx
                            Eigen::Matrix<double, local_dim, embed_dim> dxi_dx = (param_grad(q_k).transpose() * param_grad(q_k)).inverse() * param_grad(q_k).transpose();

                            // Compute pullback gradient

                            Eigen::Matrix<double, embed_dim, 1> trial_grad_ , test_grad_ ;
                            trial_grad_.setZero();
                            test_grad_.setZero();

                            for(int k = 0; k < embed_dim; ++k) {
                                for(int l = 0; l < local_dim; ++l) {
                                    trial_grad_(k) += dxi_dx(l,k) * temp_trial_grad(l);
                                    test_grad_(k)  += dxi_dx(l,k) * temp_test_grad(l);
                                }
                            }
                            // project onto the tangent plane
                            auto phys_trial_grad = P * trial_grad_;
                            auto phys_test_grad  = P * test_grad_;

                            // assign to 



                            /*

                            //iso_packet.param_grad = param_grad(q_k);
                            // convert to eigen matrix
                            Eigen::Matrix<double, local_dim, 1> trial_grad, test_grad;

                            for(int k = 0; k < local_dim; ++k) {
                                trial_grad(k) =  temp_trial_grad(k);
                                test_grad(k)  = temp_test_grad(k);
                            }
                            // convert to physical gradient
                            // print the dimensions of grad_transf(q_k)
                            //std::cout << "Grad transf: " << grad_transf(q_k).rows() << " x " << grad_transf(q_k).cols() << std::endl;
                            Eigen::Matrix<double, 3, 1> phys_trial_grad = grad_transf(q_k) *  trial_grad;
                            Eigen::Matrix<double, 3, 1> phys_test_grad  = grad_transf(q_k) * test_grad;
                            //std::cout << "Trial grad: " << i <<": "<< phys_trial_grad.transpose() << std::endl;
                            */

                            // assign to iso_packet
                            //iso_packet.trial_grad.resize(embed_dim);
                            //iso_packet.test_grad.resize(embed_dim);
                            for(int k = 0; k < embed_dim; ++k) {
                                iso_packet.trial_grad(k) = phys_trial_grad(k);
                                iso_packet.test_grad(k)  = phys_test_grad(k);
                            }
                            //std::cout << "Test grad: "<< j <<": " << phys_test_grad.transpose() << std::endl;

                        }
                        // print Form::XprBits bitmask
                        //std::cout << "Form::XprBits = 0x" << std::hex << Form::XprBits << std::endl;
                        if constexpr (Form::XprBits & int(iso_assembler_flags::compute_shape_hess)) {
                            auto temp_trial_grad = trial_param_shape_grads.template slice<0,1>(i, q_k); 
                            auto temp_test_grad  = test_param_shape_grads.template slice<0,1>(j, q_k);
                            Eigen::Matrix<double, local_dim, 1> trial_grad, test_grad;
                            for(int k = 0; k < local_dim; ++k) {
                                trial_grad(k) =  temp_trial_grad(k);
                                test_grad(k)  = temp_test_grad(k);
                            }
                            auto param_trial_hess = trial_param_shape_hess.template slice<0,1>(i, q_k);
                            auto param_test_hess  = test_param_shape_hess.template slice<0,1>(j, q_k);
                            Eigen::Matrix<double, local_dim, local_dim> trial_hess, test_hess;
                            for(int k = 0; k < local_dim; ++k) {
                                for(int l = 0; l < local_dim; ++l) {
                                    trial_hess(k,l) = param_trial_hess(k,l);
                                    test_hess(k,l)  = param_test_hess(k,l);
                                    //std::cout << "Trial hess: " << k << ": " << trial_hess(k,l) << std::endl;
                                }
                            }
                            
                            // --- START curvature-corrected Hessian logic ---
                            // Compute dξ/dx (inverse Jacobian)
                            Eigen::Matrix<double, local_dim, embed_dim> dxi_dx = (param_grad(q_k).transpose() * param_grad(q_k)).inverse() * param_grad(q_k).transpose();
                            // Compute pullback Hessian
                            Eigen::Matrix<double, embed_dim, embed_dim> phys_trial_hess;
                            Eigen::Matrix<double, embed_dim, embed_dim> phys_test_hess;
                            //phys_trial_hess = param_grad(q_k) * trial_hess * param_grad(q_k).transpose();
                            //phys_test_hess  = param_grad(q_k) * test_hess * param_grad(q_k).transpose();

                            // Project onto the tangent plane
                            Eigen::Matrix<double, embed_dim, 1> n = ((param_grad(q_k).col(0)).cross(param_grad(q_k).col(1))).normalized();
                            Eigen::Matrix<double, embed_dim, embed_dim> P = Eigen::Matrix<double, embed_dim, embed_dim>::Identity() - n * n.transpose();
                            
                            phys_trial_hess.setZero();
                            phys_test_hess.setZero();
                            // Add curvature correction term
                            for (int ii = 0; ii < embed_dim; ++ii) {
                                for (int jj = 0; jj < embed_dim; ++jj) {
                                    for(int k =0 ; k<local_dim; ++k){
                                        for(int n =0 ; n<local_dim; ++n){
                                            phys_trial_hess(ii,jj) += dxi_dx(n,ii) * dxi_dx(k,jj) * trial_hess(k,n);
                                            phys_test_hess(ii,jj)  += dxi_dx(n,ii) * dxi_dx(k,jj) * test_hess(k,n);
                                        }

                                    }

                                    for (int alpha = 0; alpha < local_dim; ++alpha) {
                                        double d2xi = 0.0;
                                        for (int kk = 0; kk < embed_dim; ++kk) {
                                            for (int beta = 0; beta < local_dim; ++beta) {
                                                for (int gamma = 0; gamma < local_dim; ++gamma) {
                                                    //std::cout<<"alpha: "<<alpha<<", beta: "<<beta<<", gamma: "<<gamma<<std::endl;
                                                    //std::cout<<param_hess(q_k)(kk, beta, gamma)<<std::endl;
                                                    d2xi -= dxi_dx(alpha, kk) * param_hess(q_k)(kk, beta, gamma) * dxi_dx(beta, ii) * dxi_dx(gamma, jj);
                                                }
                                            }
                                        }
                                        phys_trial_hess(ii, jj) += trial_grad(alpha) * d2xi;
                                        phys_test_hess(ii, jj)  += test_grad(alpha)  * d2xi;
                                    }
                                }
                            }
                                

                            auto phys_trial_hess_ = P * phys_trial_hess * P; //P * phys_trial_hess * P;
                            auto phys_test_hess_  = P * phys_test_hess * P; ; //P * phys_test_hess * P;
                            
                            /*
                            std::cout<<" ----------- "<<std::endl;

                            std::cout<<"Trial trace 1: "<<phys_trial_hess_.trace()<<std::endl;
                            std::cout<<"Test trace 1: "<<phys_test_hess_.trace()<<std::endl;

                            //std::cout<<"P: "<<P<<std::endl;
                            */

                            

                            /*
                            // other method compute directly the laplacian
                            auto trial_transf = metric_dets(q_k) *  ((param_grad(q_k).transpose() * param_grad(q_k)).inverse())* trial_grad;
                            auto test_transf  = metric_dets(q_k) *  ((param_grad(q_k).transpose() * param_grad(q_k)).inverse() * test_grad);
                            //std::cout << "Trial transf: " << trial_transf.transpose() << std::endl;

                            auto basis = test_space_->basis();


                            auto der_trial0 = basis[test_active_dofs[i]].derive(0);
                            auto der_test0  = basis[test_active_dofs[j]].derive(0);
                            auto der_trial1 = basis[test_active_dofs[i]].derive(1);
                            auto der_test1  = basis[test_active_dofs[j]].derive(1);
                            
                            double trial_div = der_trial0(trial_transf) + der_trial1(trial_transf);
                            double test_div  = der_test0(test_transf) + der_test1(test_transf);
                            //std::cout << "Trial div: " << trial_div << std::endl;
                            //std::cout << "Test div: " << test_div << std::endl;

                            for(int k = 0; k < embed_dim; ++k) {
                                for(int l = 0; l < embed_dim; ++l) {
                                    iso_packet.trial_hess(k,l) = 0;
                                    iso_packet.test_hess(k,l)  = 0;
                                }
                            }

                            // per farlo funzionare
                            iso_packet.trial_hess(0,0) = trial_div/metric_dets(q_k);
                            iso_packet.test_hess(0,0)  = test_div/metric_dets(q_k);

                            //std::cout<<"Trial trace 2: "<<trial_div/metric_dets(q_k)<<std::endl;
                            //std::cout<<"Test trace 2: "<<test_div/metric_dets(q_k)<<std::endl;
                            */
                            
                            

                            // assign to iso_packet
                            
                            for(int k = 0; k < embed_dim; ++k) {
                                for(int l = 0; l < embed_dim; ++l) {
                                    iso_packet.trial_hess(k,l) = phys_trial_hess_(k,l) ;
                                    iso_packet.test_hess(k,l)  = phys_test_hess_(k,l) ;
                                }
                            }
                                
                                
                                
                            // --- END curvature-corrected Hessian logic ---
                            //std::cout << "Trial hess: " << i << ": " << phys_trial_hess << std::endl;
                        }


                        
                        if constexpr (Form::XprBits & int(iso_assembler_flags::compute_physical_quad_nodes)) {
                            iso_packet.quad_node_id = local_cell_id * Base::n_quadrature_nodes_ + q_k;
                        }
                        value += Base::quad_weights_(q_k, 0) * form_(iso_packet) * metric_dets(q_k);
                    }

                    
                    //std::cout<<"Couple: "<<reduced_dof_map[test_active_dofs[j]]<<", "<<reduced_dof_map[is_galerkin ? test_active_dofs[i] : trial_active_dofs[i]]<<std::endl;
                    triplet_list.emplace_back(
                        test_active_dofs[j],
                        is_galerkin ? test_active_dofs[i] : trial_active_dofs[i],
                        value * iso_packet.cell_measure);
                    

                }   

            }
            local_cell_id++;
        }
        return;

    }

    constexpr int n_dofs() const { return trial_dof_handler()->n_dofs(); }
    constexpr int rows() const { return test_dof_handler()->n_dofs(); }
    constexpr int cols() const { return trial_dof_handler()->n_dofs(); }
    constexpr const TrialSpace& trial_space() const { return *trial_space_; } 
    
};



} // namespace internals
} // namespace fdapde

#endif // __FDAPDE_ISO_BILINEAR_FORM_ASSEMBLER_H__