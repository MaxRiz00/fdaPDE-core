#ifndef __FDAPDE_ISO_CUBE_H__
#define __FDAPDE_ISO_CUBE_H__


namespace fdapde {

template <typename MeshType> class IsoCube: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 2, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using Base = IsoCell<MeshType::local_dim, MeshType::embed_dim>;
    public:
    // constructor
    IsoCube() = default;
    IsoCube(int id, const MeshType* mesh) : id_(id), mesh_(mesh), boundary_(false) {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        this->left_coords_ = mesh_->compute_lr_vertices_(id_)[0];
        this->right_coords_ = mesh_->compute_lr_vertices_(id_)[1];
        // initialize = (){}; // da capire cosa inizializzare
    }

    /*

    class EdgeType : public IsoCell<MeshType::local_dim, MeshType::embed_dim>::BoundaryCellType {
        private:
        int edge_id_;
        const MeshType* mesh_;
        public:
        using CoordsType = Eigen::Matrix<double, MeshType::embed_dim, MeshType::local_dim>; 
        EdgeType() = default; //// da capire come adattare!!!
        EdgeType(int edge_id, const MeshType* mesh) : edge_id_(edge_id), mesh_(mesh) {
            for (int i = 0; i < this->n_nodes; ++i) { this->coords_.col(i) = mesh_->node(mesh_->edges()(edge_id_, i)); }
            //this->initialize();
        }
    };
    */

   // Affine map from reference domain [-1, 1]^M to parametric domain [left_coords, right_coords]^M
    // left_coords
    // map from refernce to parameric domain, map_to_parametric, left_coord e right_coord li prende dalla mesh
    Eigen::Matrix<double, MeshType::embed_dim, 1> parametrization(const std::array<double, MeshType::local_dim>& p) const {
        return mesh_->eval_param(affine_map(p));
    }

    Eigen::Matrix<double, MeshType::embed_dim, MeshType::local_dim, Eigen::RowMajor> parametrization_gradient(const std::array<double, MeshType::local_dim>& p) const {
        return mesh_->eval_param_derivative(affine_map(p));
    }

    // Metric tensor F^T * F
    Eigen::Matrix<double, MeshType::local_dim, MeshType::local_dim, Eigen::RowMajor> metric_tensor(const std::array<double, MeshType::local_dim>& p) const {
        auto F = parametrization_gradient(affine_map(p));
        return F.transpose() * F; 
    }

    // metric determinant sqrt(det(F^T * F)), array diventano matrici eigen
    double metric_determinant(const std::array<double, MeshType::local_dim>& p) const {
        return std::sqrt(metric_tensor(affine_map(p)).determinant()); 
    }


    //getters 
    int id() const { return id_; }
    Eigen::Matrix<int, 1, 2 * MeshType::local_dim> neighbors() const { return mesh_->neighbors().row(id_); }
    Eigen::Matrix<int, 1, MeshType::local_dim> node_ids() const { return mesh_->cells().row(id_); }
    bool on_boundary() const { return boundary_; }
    operator bool() const { return mesh_ != nullptr; }

    protected:
    int id_ = 0;   // segment ID in the physical mesh
    const MeshType* mesh_ = nullptr;
    bool boundary_ = false;   // true if the element has at least one vertex on the boundary
};
    
    
}; // namespace fdapde

#endif // __FDAPDE_ISO_SQUARE_H__