#ifndef __FDAPDE_ISO_CELL_H__
#define __FDAPDE_ISO_CELL_H__


namespace fdapde {

// Md Hypercube
template<int LocalDim_, int EmbedDim_> class IsoCell{
    static_assert(LocalDim_ >= 0 && LocalDim_ <= 3);
    public:
    static constexpr int local_dim = LocalDim_;
    static constexpr int embed_dim = EmbedDim_;
    static constexpr int n_nodes = 1 << LocalDim_;
    static constexpr int n_edges = LocalDim_ * (1 << (LocalDim_ - 1));
    static constexpr int n_faces = LocalDim_ * (LocalDim_ - 1) / 2 * (1 << (LocalDim_ - 2));
    static constexpr int n_nodes_per_face = 1 << (LocalDim_ - 1);
    using BoundaryCellType = std::conditional_t<LocalDim_ == 0, IsoCell<0, EmbedDim_>, IsoCell<LocalDim_ - 1, EmbedDim_>>;
    using NodeType = Eigen::Matrix<double, embed_dim, 1>;

    IsoCell() = default;

    IsoCell(const std::array<int, local_dim>& left_coords, const std::array<int, local_dim>& right_coords ): 
        left_coords_(left_coords), right_coords_(right_coords){} 


    std::array<double, LocalDim_> affine_map(const std::array<double, LocalDim_> & p) const {
            std::array<double, LocalDim_> x;
            for(std::size_t i = 0; i < LocalDim_; ++i){
                x[i] = 0.5*(right_coords_[i] + left_coords_[i] + (right_coords_[i] - left_coords_[i]) * p[i]);
            }
            return x;
        }

    double parametric_measure() const {
        double measure = 1.0;
        for(std::size_t i = 0; i < LocalDim_; i++){
            measure *= right_coords_[i] - left_coords_[i];
        }
        return measure/(1<<LocalDim_);
    }

    protected:

    std::array<int, LocalDim_> left_coords_ ; // coordinates of the left corner of the element
    std::array<int, LocalDim_> right_coords_ ; // coordinates of the right corner of the element
    // capire se salvare altre quantità

};


};

#endif // __FDAPDE_ISO_CELL_H__