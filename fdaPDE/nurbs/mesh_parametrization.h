#ifndef __MESH_PARAMETRIZATION_H__
#define __MESH_PARAMETRIZATION_H__

#include "nurbs_basis.h"


namespace fdapde{
template<int M, int N>
class MeshParametrization: public MatrixBase<M,MeshParametrization<M, N>> {
    
    private:
    NurbsBasis<M> nurbs_basis_;
    MdArray<double,full_dynamic_extent_t<M>> control_points_;

    public:
    MeshParametrization() = default;

    //constructor with i (as in curti de gaspari)

    MeshParametrization(std::array<std::vector<double>,M>& knots,MdArray<double,full_dynamic_extent_t<M>>& weights, 
        MdArray<double,full_dynamic_extent_t<M>>& control_points, int order): 
        nurbs_basis_(knots, weights, order), control_points_(control_points) {};

    inline double operator() (const std::array<double,M>& u) const {
        double x = 0.0;
        for(const auto &nurb : nurbs_basis_){
            x += nurb(x)*control_points_(nurb.index());
        }
        return x;
    };

};
} // namespace fdapde

#endif // __MESH_PARAMETRIZATION_H__