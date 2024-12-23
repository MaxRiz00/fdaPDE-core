
#ifndef __NURBS_H__
#define __NURBS_H__

#include "../linear_algebra/constexpr_matrix.h"
#include "../linear_algebra/mdarray.h"
#include "../fields/matrix_field.h"
#include "scalar_field.h"
#include "../splines/spline_basis.h"

namespace fdapde{

    // Fai un check computazionale confrontando le nurbs evaluation con quelle di  DeGaspari

// multi-contract function (iterative version)    
template<int M>
inline double multicontract(const MdArray<double, full_dynamic_extent_t<M>>& weights,
                     const std::array<std::vector<double>, M>& parts) {

    double contracted_value = 0.0;
    std::array<int, M> sizes = {0};
    std::size_t total_size = 1;

    for (int i = 0; i < M; ++i) {
        sizes[i] = parts[i].size();
        total_size *= sizes[i];
    }

    for(int flat_idx = 0; flat_idx < total_size; ++flat_idx) {
        std::array<int, M> multi_idx = {0};
        int temp = flat_idx;
        for (int i = M - 1; i >= 0; --i) {
            multi_idx[i] = temp % sizes[i];
            temp /= sizes[i];
        }
        double sub_value = weights(multi_idx);
        for (int i = 0; i < M; ++i) {
            sub_value *= parts[i][multi_idx[i]];
        }
        contracted_value += sub_value;
    }

    return contracted_value;
}



template<int M>
class Nurbs: public ScalarBase<M,Nurbs<M>> {
        public:

            using Base = ScalarBase<M,Nurbs<M>>;
            static constexpr int StaticInputSize = M;
            static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporaries
            static constexpr int XprBits = 0;
            static constexpr int Order = Dynamic;
            using Scalar = double;
            using InputType = cexpr::Vector<Scalar, StaticInputSize>;

        private:
            std::array<std::shared_ptr<SplineBasis>, M> spline_basis_;
            MdArray<double,full_dynamic_extent_t<M>> weights_;
            std::array<int,M> index_ ;
            int order_ = 0;

            double num0_ = 0.0;
            std::array<std::size_t, M> minIdx_;
            std::array<int, M> extents_;


        public:
            Nurbs() = default;

            Nurbs(std::array<std::vector<double>,M>& knots, MdArray<double,full_dynamic_extent_t<M>>& weights, std::array<int,M>& index, int order): 
                 index_(index), order_(order){

                // build a spline basis for each dimension

                std::array<std::size_t, M> maxIdx;


                for (std::size_t i = 0; i < M; ++i) {
                    // create a spline basis for each dimension
                    spline_basis_[i] = std::make_shared<SplineBasis>(knots[i], order_);

                    // compute the minIdx and extents for each dimension
                    minIdx_[i] = (index_[i] >= order_)? (index_[i]-order_) : 0;
                    extents_[i] = (index_[i] + order_ < weights.extent(i))? (index_[i]+order_+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]);
                    maxIdx[i] = (minIdx_[i] + extents_[i]-1);
                }

                // initialize the gradient
                for (std::size_t i = 0; i < M; ++i){
                    gradient_[i] = FirstDerivative(spline_basis_, weights, index, i);
                    }

                // todo: initialize the hessian

                // allocate space for the weights
                weights_.resize(extents_);
                // fill the block of weights with only the necessary values;
                weights_ = weights.block(minIdx_, maxIdx);
                // compute the starting numerator  
                num0_ = weights(index_); 

            };

            // constructor for the 1d knot passed as std::vector<double>, using the previous constructor
            Nurbs(std::vector<double>& knots, MdArray<double,full_dynamic_extent_t<M>>& weights, std::array<int,M>& index, int order): 
            Nurbs(std::array<std::vector<double>,M>{knots}, weights, index, order) {
                fdapde_static_assert(M == 1, THIS_METHOD_IS_ONLY_FOR_1D_NURBS);
            };

            Nurbs(std::array<std::shared_ptr<SplineBasis>, M> spline_basis, MdArray<double,full_dynamic_extent_t<M>>& weights,std::array<int,M>& index) : spline_basis_(spline_basis), index_(index) { 
                // questo viene usato da NurbsBasis

                //decalre maxIdx
                std::array<std::size_t, M> maxIdx;

                // vector of knots
                std::array<std::vector<double>, M> knots;

                order_ = spline_basis[0]->order();

                // like in the previous constructor
                for (std::size_t i = 0; i < M; ++i) {
                    // compute the minIdx and extents for each dimension
                    minIdx_[i] = (index_[i] >= order_)? (index_[i]-order_) : 0;
                    extents_[i] = (index_[i] + order_ < weights.extent(i))? (index_[i]+order_+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]);
                    maxIdx[i] = (minIdx_[i] + extents_[i]-1);
                    
                }

                // allocate for the gradient
                for (std::size_t i = 0; i < M; ++i){
                    gradient_[i] = FirstDerivative(spline_basis_, weights, index, i);
                    }

                // allocate space for the weights
                weights_.resize(extents_);
                // fill the block of weights with only the necessary values;
                weights_ = weights.block(minIdx_, maxIdx);

                // compute the starting numerator
                num0_ = weights(index_);

            // STAI SALTANDO L'INIZIALIZZAZIONE DELLA BASE SPLINE, QUINDI SE CREO 1000 Nurbs questo è più efficiente
            
            }

            //evaluates the NURBS at a given point, funziona
            template <typename InputType_>
            requires(fdapde::is_subscriptable<InputType_, int>)
            constexpr Scalar operator()(const InputType_& p_) const {
                
                double num = num0_;
                std::array<std::vector<double>,M> spline_evaluation {};
                double den;

                // Initialize each vector to the appropriate size
                //for (std::size_t i = 0; i < M; ++i) {
                //    spline_evaluation[i].resize(extents_[i]);
                //}

                for(std::size_t i=0;i<M;i++){

                    auto basis_eval = spline_basis_[i]->evaluate_basis(p_[i]);
                    spline_evaluation[i].resize(extents_[i]);
                    for(std::size_t j = 0; j<extents_[i]; j++ ){
                        // compute a spline basis function
                        spline_evaluation[i][j] = basis_eval[minIdx_[i]+j]; // rivedi 

                        // metti in spline un overload del call operator che accetta un double
                    }
                    //numerator update
                    num *= spline_evaluation[i][index_[i] - minIdx_[i]]; 
                    //spline evaluation for i-th dimension
                }
                // avoid division by 0
                if(num == 0)
                    return 0;
                // compute the sum that appears at the denominator of the formula
                den = multicontract<M>(weights_, spline_evaluation);

                return num/den;
            }; 

        private:
            class FirstDerivative: public MatrixBase<M,FirstDerivative> {
            public:
            
                using Base = MatrixBase<M,FirstDerivative>;
                static constexpr int StaticInputSize = M;
                static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporaries
                static constexpr int XprBits = 0;
                static constexpr int Order = Dynamic;
                using Scalar = double;
                using InputType = cexpr::Vector<Scalar, StaticInputSize>;

            private:
                std::array<std::shared_ptr<SplineBasis>, M> spline_basis_;
                MdArray<double,full_dynamic_extent_t<M>> weights_;
                std::array<int,M> index_ ;
                int order_;

                std::array<std::size_t, M> minIdx_;
                double num0_ = 0.0;
                std::array<int, M> extents_;

                size_t i_ = 0; // index of the derivative, questa i è necessaria ? se mi interessa il gradiente in una direzione si

            public:
                FirstDerivative() = default;

                FirstDerivative(std::array<std::shared_ptr<SplineBasis>, M> spline_basis, const MdArray<double,full_dynamic_extent_t<M>>& weights, const std::array<int,M>& index, std::size_t i): 
                     spline_basis_(spline_basis), index_(index), i_(i){

                    // build a spline basis for each dimension

                    order_ = spline_basis[0]->order();

                    std::array<std::size_t, M> maxIdx;
                    for (std::size_t i = 0; i < M; ++i) {

                        // compute the minIdx and extents for each dimension
                        minIdx_[i] = (index_[i] >= order_)? (index_[i]-order_) : 0;
                        extents_[i] = (index_[i] + order_ < weights.extent(i))? (index_[i]+order_+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]);
                        maxIdx[i] = (minIdx_[i] + extents_[i]-1);
                    }

                    // allocate space for the weights
                    weights_.resize(extents_);
                    // fill the block of weights with only the necessary values;
                    weights_ = weights.block(minIdx_, maxIdx);
                    // compute the starting numerator  
                    num0_ = weights(index_);  

                    //std::cout << "NURBS derivative correctly initialized! " << std::endl;

                };

                // evalutes the first order partial derivative of the NURBS at a given point, funziona
                constexpr Scalar operator()(const InputType& p) const {

                    //std::cout<<"Inizio a calcolare"<<std::endl;

                    double num = num0_;
                    std::array<std::vector<double>,M> spline_evaluation {};
                    double den;

                    double num_derived = 0.0;
                    double den_derived = 0.0;


                    for(std::size_t i=0;i<M;i++){
                        // spline evaluation for i-th dimension
                        auto basis_eval = spline_basis_[i]->evaluate_basis(p[i]);
                        spline_evaluation[i].resize(extents_[i]);
                        for(std::size_t j = 0; j<extents_[i]; j++ ){
                        // compute a spline basis function
                        spline_evaluation[i][j] = basis_eval[minIdx_[i]+j]; // rivedi 

                        // metti in spline un overload del call operator che accetta un double
                    }
                        if (i!=i_)
                            num*=spline_evaluation[i][index_[i] - minIdx_[i]];
                    
                    }

                    

                    //compute the derivative of the i_th spline
                    //num_derived = num * Spline(knots_[i_], index_[i_], order_).gradient(1)(p[i_]);
                    num_derived = num * (*spline_basis_[i_])[index_[i_]].gradient(1)(p[i_]);

                    // compute the non derived numerator
                    num*=spline_evaluation[i_][index_[i_] - minIdx_[i_]];

                    if (num== 0 && num_derived == 0)
                        return 0;


                    // compute the sum that appears at the denominator of the formula
                    den = multicontract<M>(weights_, spline_evaluation);

                    // by replacing the i-th evaluations with their derivatives we get the derivative of the NURBS denominator
                    auto der_eval = spline_basis_[i_]->evaluate_der_basis(p[i_],1);

                    for (std::size_t j = 0; j<extents_[i_]; j++ ){
                        // extract the knots
                        //auto knots = spline_basis_[i_]->knots_vector();
                        //auto spline = Spline(knots, minIdx_[i_]+j, order_);
                        spline_evaluation[i_][j] = der_eval[ minIdx_[i_]+j];
                    }

                    den_derived = multicontract<M>(weights_, spline_evaluation);

                    //  ( N )'      N'D - ND'
                    //  (---)   =  ----------
                    //  ( D )         D^2
                    // where f' = df/dx_i
                    return (num_derived*den - num*den_derived)/(den*den);
                };
            };

            private:

                // gradient and hessian$
                VectorField<M, M, FirstDerivative> gradient_; //gradient
                //MatrixField<M,M,M,NurbsSecondDerivative>> hessian_; //hessian accedi con hessian_(i,j)

            public:
                constexpr FirstDerivative first_derive(int i=0) const { return gradient_[i]; }
                //constexpr SecondDerivarive second_derive(int i, int j) const { return SecondDerivative(knots_, weights_, index_, order_, i, j); }

                // getters
                constexpr int order() const { return order_; }
                constexpr int size() const { return weights_.size(); }
                constexpr const MdArray<double,full_dynamic_extent_t<M>>& weights() const { return weights_; }
                constexpr const std::array<int,M>& index() const { return index_; }
                constexpr const std::array<std::shared_ptr<SplineBasis>, M>& spline_basis() const { return spline_basis_; }
                
    
    };

}// namespace fdapde

#endif


