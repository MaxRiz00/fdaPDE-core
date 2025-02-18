
#ifndef __NURBS_H__
#define __NURBS_H__

#include "header_check.h"

#include "linear_algebra.h"
#include "fields.h"
#include "splines.h"

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
class Nurbs: public ScalarFieldBase<M,Nurbs<M>> {
        public:

            using Base = ScalarFieldBase<M,Nurbs<M>>;
            static constexpr int StaticInputSize = M;
            static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporaries
            static constexpr int XprBits = 0;
            static constexpr int Order = Dynamic;
            using Scalar = double;
            using InputType = Vector<Scalar, StaticInputSize>;

        private:
            std::array<std::shared_ptr<BSplineBasis>, M> spline_basis_;
            MdArray<double,full_dynamic_extent_t<M>> weights_;
            std::array<int,M> index_ ;
            int order_ = 0;

            double num0_ = 0.0;
            std::array<std::size_t, M> minIdx_;
            std::array<int, M> extents_;


        public:
            Nurbs() = default;


            template <typename KnotsVectorType>
            requires(requires(KnotsVectorType knots, int i) {
                    { knots[i] } -> std::convertible_to<double>;
                    { knots.size() } -> std::convertible_to<std::size_t>;
                })
            Nurbs(std::array<KnotsVectorType,M>&& knots, MdArray<double,full_dynamic_extent_t<M>>& weights, std::array<int,M>&& index, int order): 
                 index_(std::move(index)), order_(order){
                

                // we suppose the knots are not padded
                // build a spline basis for each dimension

                std::array<std::size_t, M> maxIdx;
                
                for (std::size_t i = 0; i < M; ++i) {

                    std::vector<double> knots_ ;
                    // pad the knot vector to obtain a full basis for the whole knot span [u_0, u_n]
                    auto n=knots[i].size();
                    knots_.resize(n + 2 * order_);
                    for (std::size_t j = 0; j < n + 2 * order_; ++j) {
                        if (j < order_) {
                            knots_[j] = knots[i][0];
                        } else {
                            if (j < n + order_) {
                                knots_[j] = knots[i][j - order_];
                            } else {
                                knots_[j] = knots[i][n - 1];
                            }
                        }
                    }
                    
                    spline_basis_[i] = std::make_shared<BSplineBasis>(knots_, order_);

                    // compute the minIdx and extents for each dimension
                    minIdx_[i] = (index_[i] >= order_)? (index_[i]-order_) : 0;
                    extents_[i] = (index_[i] + order_ < weights.extent(i))? (index_[i]+order_+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]);
                    maxIdx[i] = (minIdx_[i] + extents_[i]-1);
                }

                // initialize the gradient
                for (std::size_t i = 0; i < M; ++i){
                    gradient_[i] = FirstDerivative(spline_basis_, weights, index, i);
                    }

                // initialize the hessian
                for (std::size_t i = 0; i < M; ++i){
                    for (std::size_t j = 0; j < M; ++j){
                        hessian_(i,j) = SecondDerivative(spline_basis_, weights, index, i, j);
                    }
                }

                // allocate space for the weights
                weights_.resize(extents_);
                // fill the block of weights with only the necessary values;
                weights_ = weights.block(minIdx_, maxIdx);
                // compute the starting numerator  
                num0_ = weights(index_); 

            };

            // constructor for the 1d knot passed as std::vector<double>, using the previous constructor, non so se va bene ?????
            template <typename KnotsVectorType>
            requires(requires(KnotsVectorType knots) {
                    { knots.size() } -> std::convertible_to<std::size_t>;
                })
            Nurbs(KnotsVectorType& knots, MdArray<double,full_dynamic_extent_t<M>>& weights, int index, int order): 
            Nurbs(std::array<std::vector<double>,M>{std::move(knots)}, weights, std::array<int,M>{index}, order) {
                fdapde_static_assert(M == 1, THIS_METHOD_IS_ONLY_FOR_1D_NURBS);
            };

            // constructor with shared pointers, used by the NurbsBasis
            Nurbs(std::array<std::shared_ptr<BSplineBasis>, M> spline_basis, MdArray<double,full_dynamic_extent_t<M>>& weights,std::array<int,M>& index) : spline_basis_(spline_basis), index_(index) { 
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
                /*
                std::cout<<"Knots padded ECCOCI"<<std::endl;

                // print the index
                for(int j=0;j<M;++j){
                    std::cout<<index_[j]<<" ";
                }
                std::cout<<std::endl;
                */

                // allocate for the gradient
                for (std::size_t i = 0; i < M; ++i){
                    gradient_[i] = FirstDerivative(spline_basis_, weights, index, i);
                    }
                //std::cout<<"Knots padded ECCOCI1"<<std::endl;

                // initialize the hessian
                for (std::size_t i = 0; i < M; ++i){
                    for (std::size_t j = 0; j < M; ++j){
                        hessian_(i,j) = SecondDerivative(spline_basis_, weights, index, i, j);
                    }
                }
                //std::cout<<"Knots padded ECCOCI2"<<std::endl;

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
            requires(internals::is_subscriptable<InputType_, int>)
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
            class FirstDerivative: public MatrixFieldBase<M,FirstDerivative> {
            public:
            
                using Base = MatrixFieldBase<M,FirstDerivative>;
                static constexpr int StaticInputSize = M;
                static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporaries
                static constexpr int XprBits = 0;
                static constexpr int Order = Dynamic;
                using Scalar = double;
                using InputType = Vector<Scalar, StaticInputSize>;

            private:
                std::array<std::shared_ptr<BSplineBasis>, M> spline_basis_;
                MdArray<double,full_dynamic_extent_t<M>> weights_;
                std::array<int,M> index_ ;
                int order_;

                std::array<std::size_t, M> minIdx_;
                double num0_ = 0.0;
                std::array<int, M> extents_;

                size_t i_ = 0; // index of the derivative, questa i è necessaria ? se mi interessa il gradiente in una direzione si

            public:
                FirstDerivative() = default;

                FirstDerivative(std::array<std::shared_ptr<BSplineBasis>, M> spline_basis, const MdArray<double,full_dynamic_extent_t<M>>& weights, const std::array<int,M>& index, std::size_t i): 
                     spline_basis_(spline_basis), index_(index), i_(i){

                    // build a spline basis for each dimension

                    order_ = spline_basis[0]->order();

                    std::array<std::size_t, M> maxIdx;
                    //std::cout<<"NURBS derivative initialized"<<std::endl;
                    for (std::size_t i = 0; i < M; ++i) {

                        // compute the minIdx and extents for each dimension
                        minIdx_[i] = (index_[i] >= order_)? (index_[i]-order_) : 0;
                        extents_[i] = (index_[i] + order_ < weights.extent(i))? (index_[i]+order_+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]);
                        maxIdx[i] = (minIdx_[i] + extents_[i]-1);
                    }

                    //std::cout<<"NURBS derivative initialized2"<<std::endl;

                    // allocate space for the weights
                    weights_.resize(extents_);
                    // fill the block of weights with only the necessary values;
                    weights_ = weights.block(minIdx_, maxIdx);

                    //std::cout<<"NURBS derivative initialized3"<<std::endl;
                    // compute the starting numerator  
                    
                    num0_ = weights(index_);  
                    //num0_ = weights_(index_ - minIdx_); // rivedi

                    //std::cout << "NURBS derivative correctly initialized! " << std::endl;

                };

                // evalutes the first order partial derivative of the NURBS at a given point, funziona
                constexpr Scalar operator()(const std::array<Scalar,StaticInputSize> & p) const { // attaentzione, qua devi usare InputType

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

                // overload the call operator for the first derivative for a double
                //constexpr Scalar operator()(double p) const {
                //    return operator()(InputType{p});
                //}


            class SecondDerivative: public MatrixFieldBase<M,SecondDerivative> {
            public:
            
                using Base = MatrixFieldBase<M,SecondDerivative>;
                static constexpr int StaticInputSize = M;
                static constexpr int NestAsRef = 0;   // avoid nesting as reference, .derive() generates temporaries
                static constexpr int XprBits = 0;
                static constexpr int Order = Dynamic;
                using Scalar = double;
                using InputType = Vector<Scalar, StaticInputSize>;

            private:
                std::array<std::shared_ptr<BSplineBasis>, M> spline_basis_;
                MdArray<double,full_dynamic_extent_t<M>> weights_;
                std::array<int,M> index_ ;
                int order_ = 0;

                std::array<std::size_t, M> minIdx_;
                double num0_ = 0.0;
                std::array<int, M> extents_;

                size_t i_ = 0; // first index of the 2nd derivative
                size_t j_ = 0; // second index of the 2nd derivative

            public:
                SecondDerivative() = default;

                SecondDerivative(std::array<std::shared_ptr<BSplineBasis>, M> spline_basis, const MdArray<double,full_dynamic_extent_t<M>>& weights, const std::array<int,M>& index, std::size_t i, std::size_t j): 
                     spline_basis_(spline_basis), index_(index), i_(i), j_(j){

                    // build a spline basis for each dimension

                    order_ = spline_basis[0]->order();

                    //std::cout<<"Ecco l'ordine: "<<order_<<std::endl;

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

                    //std::cout << "NURBS  2nd derivative correctly initialized! " << std::endl;

                };


                // evalutes the hessian(i,j) of the NURBS at a given point
                constexpr Scalar operator()(const InputType& p) const {

                    double num = num0_; // numerator of the NURBS formula
                    double num_der_i; // partial derivative of num w.r.t. i-th coordinate
                    double num_der_j; // partial derivative of num w.r.t. j-th coordinate
                    double num_der_ij; // mixed partial derivative of num
                    std::array<std::vector<double>,M> spline_evaluation {}; // pointwise evaluation of all splines along each coordinate
                    double den; // denominator of the NURBS formula
                    double den_der_i; // partial derivative of den w.r.t. i-th coordinate
                    double den_der_j; // partial derivative of den w.r.t. j-th coordinate
                    double den_der_ij; // mixed partial derivative of den

                    for(std::size_t i=0;i<M;i++){
                        // spline evaluation for i-th dimension
                        auto basis_eval = spline_basis_[i]->evaluate_basis(p[i]);
                        spline_evaluation[i].resize(extents_[i]);
                        
                        for(std::size_t j = 0; j<extents_[i]; j++ )
                            spline_evaluation[i][j] = basis_eval[minIdx_[i]+j]; // rivedi 
                        
                        if (i!=i_ && i!=j_)
                            num*=spline_evaluation[i][index_[i] - minIdx_[i]];
                    }

                    if (i_!=j_){
                        auto der_i = (*spline_basis_[i_])[index_[i_]].gradient(1)(p[i_]);
                        auto der_j = (*spline_basis_[j_])[index_[j_]].gradient(1)(p[j_]);
                        num_der_i = num * der_i * spline_evaluation[j_][index_[j_] - minIdx_[j_]];
                        num_der_j = num * der_j * spline_evaluation[i_][index_[i_] - minIdx_[i_]];
                        num_der_ij = num * der_i * der_j;

                        num*=spline_evaluation[i_][index_[i_] - minIdx_[i_]]*spline_evaluation[j_][index_[j_] - minIdx_[j_]];

                        if (num== 0 && num_der_i == 0 && num_der_j == 0 && num_der_ij == 0)
                            return 0;
                        
                        // denominator evaluation
                        den = multicontract<M>(weights_, spline_evaluation);

                        // by replacing the i-th evaluations with their derivatives we get the derivative of the NURBS denominator
                        auto der_eval_i = spline_basis_[i_]->evaluate_der_basis(p[i_],1);
                        auto der_eval_j = spline_basis_[j_]->evaluate_der_basis(p[j_],1);

                        auto spline_eval_temp = spline_evaluation;

                        for (std::size_t j = 0; j<extents_[i_]; j++ ){
                            // extract the knots
                            spline_eval_temp[i_][j] = der_eval_i[ minIdx_[i_]+j];
                        }

                        // compute the derivative of the denominator w.r.t. i-th coordinate
                        den_der_i = multicontract<M>(weights_, spline_eval_temp);

                        spline_eval_temp = spline_evaluation;

                        for (std::size_t j = 0; j<extents_[j_]; j++ ){
                            // extract the knots
                            spline_eval_temp[j_][j] = der_eval_j[ minIdx_[j_]+j];
                        }
                        
                        // compute the derivative of the denominator w.r.t. j-th coordinate
                        den_der_j = multicontract<M>(weights_, spline_eval_temp);


                        for (std::size_t j = 0; j<extents_[i_]; j++ ){
                            // extract the knots
                            spline_eval_temp[i_][j] = der_eval_i[ minIdx_[i_]+j];
                        }

                        // compute the mixed partial derivative of the denominator
                        den_der_ij = multicontract<M>(weights_, spline_eval_temp);
                    }


                    else{
                        num_der_i = num_der_j = num * (*spline_basis_[i_])[index_[i_]].gradient(1)(p[i_]);

                        num_der_ij = num * (*spline_basis_[i_])[index_[i_]].gradient(2)(p[i_]);

                        num*=spline_evaluation[i_][index_[i_] - minIdx_[i_]];
                        if (num== 0 && num_der_i == 0  && num_der_ij == 0)
                            return 0;

                        // denominator evaluation    
                        den = multicontract<M>(weights_, spline_evaluation);

                        auto der_eval_i = spline_basis_[i_]->evaluate_der_basis(p[i_],1);

                        auto spline_eval_temp = spline_evaluation;

                        for (std::size_t j = 0; j<extents_[i_]; j++ ){
                            // extract the knots
                            spline_eval_temp[i_][j] = der_eval_i[ minIdx_[i_]+j];
                
                        }
                        // compute the derivative of the denominator w.r.t. i-th coordinate
                        den_der_i = den_der_j = multicontract<M>(weights_, spline_eval_temp);

                        auto der_eval_ij = spline_basis_[i_]->evaluate_der_basis(p[i_],2);

                        for (std::size_t j = 0; j<extents_[i_]; j++ ){
                            // extract the knots
                            spline_eval_temp[i_][j] = der_eval_ij[ minIdx_[i_]+j];
                        }

                        // compute the mixed (2nd derivative on i-th coordinate) partial derivative of the denominator
                        den_der_ij = multicontract<M>(weights_, spline_eval_temp);

                    }

                    //  ( N )'°     D(N'°D - N'D° - N°D' -ND'°) + 2D'D°N
                    //  (---)   =   ------------------------------------
                    //  ( D )                       D^3
                    // where f' = df/dx_i and f° = df/dx_j

                    return (den*(num_der_ij*den - num_der_i*den_der_j - num_der_j*den_der_i - num*den_der_ij) + 2*den_der_i*den_der_j*num)/(den*den*den);

                }

            };
            private:

                // gradient and hessian$
                VectorField<M, M, FirstDerivative> gradient_; //gradient
                MatrixField<M,M,M, SecondDerivative> hessian_; //hessian accedi con hessian_(i,j)

            public:
                constexpr FirstDerivative derive(int i=0) const { return gradient_[i]; }
                constexpr SecondDerivative deriveTwice(int i=0, int j=0) const { return hessian_(i,j); }

                // getters
                constexpr int order() const { return order_; }
                constexpr int size() const { return weights_.size(); }
                constexpr const MdArray<double,full_dynamic_extent_t<M>>& weights() const { return weights_; }
                constexpr const std::array<int,M>& index() const { return index_; }
                constexpr const std::array<std::shared_ptr<BSplineBasis>, M>& spline_basis() const { return spline_basis_; }

                constexpr Scalar operator()(double p) const { return operator()(std::vector<double>{p}); }
    
    };

}// namespace fdapde

#endif


