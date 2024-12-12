
#ifndef __NURBS_H__
#define __NURBS_H__

#include "../linear_algebra/constexpr_matrix.h"
#include "../linear_algebra/mdarray.h"
#include "scalar_field.h"
#include "spline.h"

#include <tuple>

namespace fdapde{

// Funziona ma è migliorabile, magari una versione iterativa? Comunque la ricorsione non fa più di 3 passaggi
template <int N, int J = 0>
inline double multicontract(
    const MdArray<double, full_dynamic_extent_t<N - J>> &weights,
    const std::array<std::vector<double>, N> &parts
) {
    if constexpr (N == J + 1) {
        // Base case sum_{i=0}^{N-1} w_i * p_i
        double contracted_value = 0.0;
        for (int i = 0; i < weights.extent(0); ++i) {
            contracted_value += weights(i) * parts[J][i];
        }
        return contracted_value;
    } else {
        // Perform contraction along dimension J
        auto extent_J = weights.extent(J);
        //fdapde_assert(parts[J].extent(0) == extent_J); non funzione xk non ho l'extent

        double contracted_value = 0.0;

        for (int idx = 0; idx < extent_J; ++idx) {
            auto sub_block = weights.template slice<0>(idx);
            contracted_value += parts[J][idx] * multicontract<N, J+1 >(sub_block, parts);
        }
        return contracted_value;
    }
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
            std::array<std::vector<double>,M> knots_ ;
            MdArray<double,full_dynamic_extent_t<M>> weights_;
            std::array<int,M> index_ ;
            int order_;

            std::array<std::size_t, M> minIdx_;
            std::array<Eigen::Index, M> extents_; // fare il min_idx e max idx ??? al posto di minIdx e extents
            double num0_;

        public:
            Nurbs() = default;

            Nurbs(std::array<std::vector<double>,M>& knots, MdArray<double,full_dynamic_extent_t<M>>& weights, std::array<int,M>& index, int order): 
                knots_(knots), index_(index), order_(order){

                for (std::size_t i = 0; i < M; ++i) {
                    minIdx_[i] = (index_[i] >= order)? (index_[i]-order) : 0;
                    extents_[i] = (index_[i] + order < weights.extent(i))? (index_[i]+order+1-minIdx_[i]) : (weights.extent(i)-minIdx_[i]); 

                }

                // allocate space for the weights
                std::apply([&](auto... dims) { weights_.resize(dims...); }, extents_);

                // Create pairs of min and max indices 
                std::array<std::pair<size_t, size_t>, M> pairs; 
                for (size_t i = 0; i < M; ++i) {
                    pairs[i] = {minIdx_[i], minIdx_[i] + extents_[i]-1};
                    std::cout<<"Pair["<<i<<"] = "<<pairs[i].first<<","<<pairs[i].second<<std::endl;
                }

                
                // fill the block of weights with only the necessary values
                weights_ = std::apply([&](auto... slicers) { return weights.block(slicers...); }, pairs);

                // compute the starting numerator
                num0_ = std::apply([&](auto... idx) { return weights_(idx...); }, index_);      
                
                // alloca per il gradiente 
                // alloca per l'hessiana  

                std::cout << "NURBS correctly initialized! " << std::endl; 

            };
            

            //evaluates the NURBS at a given point, funziona
            template <typename InputType_>
            requires(fdapde::is_subscriptable<InputType_, int>)
            constexpr Scalar operator()(const InputType_& p_) const {
                
                double num = num0_;
                std::array<std::vector<double>,M> spline_evaluation {};
                double den;

                // Initialize each vector to the appropriate size
                for (std::size_t i = 0; i < M; ++i) {
                    spline_evaluation[i].resize(extents_[i]);
                }

                for(std::size_t i=0;i<M;i++){

                    // pad the spline evaluation with order + 1 repetitions of the first and last knot
                    int n = knots_[i].size();
                    std::vector<double> used_knots;
                    used_knots.resize(n + 2 * order_);
                    // pad the knot vector to obtain a full basis for the whole knot span [knots[0], knots[n-1]]
                    for (int k = 0; k < n + 2 * order_; ++k) {
                        if (k < order_) {
                            used_knots[k] = knots_[i][0];
                        } else {
                            if (k < n + order_) {
                                used_knots[k] = knots_[i][k - order_];
                            } else {
                                used_knots[k] = knots_[i][n - 1];
                            }
                        }
                    }

                    for(std::size_t j = 0; j<extents_[i]; j++ ){
                        auto spline = Spline(used_knots, minIdx_[i]+j, order_);
                        spline_evaluation[i][j] = spline(std::vector<double>{p_[i]}); 
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

    };
}

#endif


