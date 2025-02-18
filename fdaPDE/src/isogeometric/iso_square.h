#ifndef __FDAPDE_ISO_SQUARE_H__
#define __FDAPDE_ISO_SQUARE_H__


namespace fdapde {

template <typename MeshType> class IsoSquare: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 2, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using Base = IsoCell<MeshType::local_dim, MeshType::embed_dim>;
    public:
    // constructor
    IsoSquare() = default;
    IsoSquare(int id, const MeshType* mesh) : id_(id), mesh_(mesh), boundary_(false) {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        //prova se sta in una riga
        auto [left_coords, right_coords] = mesh_->compute_lr_vertices_(id_);
        this->left_coords_ = left_coords;
        this->right_coords_ = right_coords;
        // initialize = (){}; // da capire cosa inizializzare
    }

    // view of an edge
    class EdgeType : public IsoCell<MeshType::local_dim, MeshType::embed_dim>::BoundaryCellType{
        private:
        int edge_id_;
        const MeshType* mesh_;
        public:
        EdgeType() = default;
        EdgeType(int edge_id, const MeshType* mesh): edge_id_(edge_id), mesh_(mesh){
           this->left_coords_[0] =  mesh_->parametric_nodes()(mesh_->edges()(edge_id,0), 0);
           this->right_coords_[0] = mesh_->parametric_nodes()(mesh_->edges()(edge_id,0), 1);
        }
        bool on_boundary() const { return mesh_->is_edge_on_boundary(edge_id_);}
        Eigen::Matrix<int, Dynamic, 1> node_ids() const { return mesh_->edges().row(edge_id_); }
        int id() const { return edge_id_; }
        Eigen::Matrix<int, Dynamic, 1> adjacent_cells() const { return mesh_->edge_to_cells().row(edge_id_); }
        int marker() const {   // mesh edge's marker
            return mesh_->edges_markers().size() > edge_id_ ? mesh_->edges_markers()[edge_id_] : Unmarked;
        }
    };

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
    EdgeType edge(int n){
        fdapde_assert(n<this->n_edges);
        return EdgeType(mesh_->cell_to_edes()(id_,n), mesh_);
    }

    // cell marker
    int marker() const {return mesh_->cell_markers().size() ? mesh_->cells_markers()[id_] : Unmarked; }

    class edge_iterator: public internals::index_iterator<edge_iterator, EdgeType>{
        using Base = internals::index_iterator<edge_iterator, EdgeType>;
        using Base::index_;
        friend Base;
        const IsoSquare* sq_;
        // access to i-th square edge
        edge_iterator& operator() (int i) {
            Base::val_ = sq_->edge(i);
            return *this;
        }
        public:
        edge_iterator(int index, const IsoSquare* sq){
            if(index_ < sq_->n_edges) operator()(index_);
        }

    };

    edge_iterator edges_begin() const {return edge_iterator(0, this);}
    edge_iterator edges_end() const {return edge_iterator(this->n_edges, this);}

    protected:
    int id_ = 0;   // segment ID in the physical mesh
    const MeshType* mesh_ = nullptr;
    bool boundary_ = false;   // true if the element has at least one vertex on the boundary
};
    
    
}; // namespace fdapde

#endif // __FDAPDE_ISO_SQUARE_H__