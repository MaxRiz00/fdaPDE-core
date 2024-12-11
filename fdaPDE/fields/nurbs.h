
#ifndef __NURBS_H__
#define __NURBS_H__

#include "../linear_algebra/constexpr_matrix.h"
#include "../linear_algebra/mdarray.h"
#include "scalar_field.h"
#include "spline.h"


namespace fdapde{

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
            std::array<std::vector<double>,M> knots_ ;
            //Eigen::Tensor<double,M> weights_;
            MdArray<double,full_dynamic_extent_t<M>> weights_;
            std::array<M, int> index_ ;
            int order_;

            std::array<std::size_t, M> minIdx_;
            Eigen::Array<Eigen::Index, M> extents_;
            double num0_;

        public:
            Nurbs() = default;

            Nurbs(const std::array<std::vector<double>,M>& knots, const MdArray<double,full_dynamic_extent_t<M>>& weights, const std::array<M, int>& index, int order): 
                knots_(knots), index_(index), order_(order){

                for (std::size_t i = 0; i < M; ++i) {
                    minIdx_[i] = (index_[i] >= order)? (index_[i]-order) : 0;
                    extents_[i] = (index_[i] + order < weights.extent(i))? (index_[i]+order+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]); // cambiare dimension

                }

                // allocate space for the weights
                weights_.resize(extents_);

                // Store the pairs
                std::vector<std::pair<size_t, size_t>> pairs;
                // Create pairs
                for (size_t i = 0; i < M; ++i) {
                    pairs.emplace_back(minIdx_[i], minIdx_[i] + extents_[i]);
                }

                //per ritornare un blocco, cioè un MdArray con lo stesso ordine dell'MdArray di partenza: auto blk = md.block({a1, a2}, {b1, b2});
                // copy the block of weights
                weights_ = weights.block(pairs);
                num0_ = weights(index_);

                //weights_ = weights.slice(minIdx_, extents_);
                //num0_ = weights(index_);

                // alloca per il gradiente 
                // alloca per l'hessiana

            };
            /*

            //evaluates the NURBS at a given point 
            inline double operator()(const SVector<M>& x) const {
                
                double num = num0_;
                std::array<std::vector<double>,M> spline_evaluation;
                double den;

                //builds the M-dim SVector containing in each position the set of spline basis evaluated
                // and compute the NURBS numerator
                for(std::size_t i=0;i<M;i++){
                    //resize i-th tensor according to i-th weights dimension

                    for(std::size_t j = 0; j<extents_[i]; j++ ){
                        spline_evaluation[i][j] = Spline(knots_[i],min_Idx[i]+j, order_) (x);
                    }

                    //numerator update
                    num *= spline_evaluation[i][index_[i] - minIdx_[i]];
                    //spline evaluation for i-th dimension
                }
                // avoid division by 0
                if(num == 0)
                    return 0;
                // compute the sum that appears at the denominator of the formula
                
                den = multicontract<M,0>(weights_, spline_evaluation);
                return num/den;
            }; */

            

    };
}


