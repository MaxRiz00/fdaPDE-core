#ifndef __FDAPDE_ISO_SQUARE_H__
#define __FDAPDE_ISO_SQUARE_H__

#include "header_check.h"

namespace fdapde {

template <typename MeshType> class IsoSquare: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 2, THIS_CLASS_IS_FOR_2D_MESHES_ONLY);
    using Base = IsoCell<MeshType::local_dim, MeshType::embed_dim>;
    public:
    // constructor
    IsoSquare() = default;
    IsoSquare(int id, const MeshType* mesh) : IsoCell<MeshType::local_dim, MeshType::embed_dim>(
            mesh->compute_lr_vertices(id)[0],  // left_coords
            mesh->compute_lr_vertices(id)[1])   // right_coords
        , id_(id), mesh_(mesh), boundary_(false)  {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        //auto [left_coords, right_coords] = mesh_->compute_lr_vertices(id_);
        // print left_coords and right_coords
        //this->left_coords_ = left_coords;

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
           this->left_coords_(0) =  mesh_->parametric_nodes()(mesh_->edges()(edge_id,0), 0);
           this->right_coords_(0) = mesh_->parametric_nodes()(mesh_->edges()(edge_id,0), 1);
        }
        bool on_boundary() const { return mesh_->is_edge_on_boundary(edge_id_);}
        Eigen::Matrix<int, Dynamic, 1> node_ids() const { return mesh_->edges().row(edge_id_); }
        int id() const { return edge_id_; }
        Eigen::Matrix<int, Dynamic, 1> adjacent_cells() const { return mesh_->edge_to_cells().row(edge_id_); }
        int marker() const {   // mesh edge's marker
            return mesh_->edges_markers().size() > edge_id_ ? mesh_->edges_markers()[edge_id_] : Unmarked;
        }

        Eigen::Matrix<double, Eigen::Dynamic, MeshType::embed_dim> evaluation(int n) const {
            Eigen::Matrix<double, Eigen::Dynamic, MeshType::embed_dim> res(n, MeshType::embed_dim);
            auto nodes = node_ids(); // Expected to be Eigen::Matrix<int, 2, 1>
            Eigen::Matrix<double, Eigen::Dynamic, MeshType::local_dim> parametric_nodes = mesh_->parametric_nodes();
            Eigen::Matrix<double, 1, MeshType::local_dim> n1 = parametric_nodes.row(nodes(0));
            Eigen::Matrix<double, 1, MeshType::local_dim> n2 = parametric_nodes.row(nodes(1));

            Eigen::Matrix<double, Eigen::Dynamic, MeshType::local_dim> interpolated_points(n, MeshType::local_dim);
            for (int i = 0; i < n; ++i) {
                double t = static_cast<double>(i) / (n - 1);  // Normalized parameter (0 to 1)
                interpolated_points.row(i) = (1 - t) * n1 + t * n2;  // Linear interpolation
            }
            // Evaluate mapped coordinates in the embedded space
            for (int i = 0; i < n; ++i) {
                Eigen::Matrix<double,  MeshType::local_dim,1> p;
                for (int j = 0; j < MeshType::local_dim; ++j) {
                    p(j) = interpolated_points(i, j);
                }
                res.row(i) = mesh_->eval_param(p);
            }
            return res;
        }
        
    };

    // Affine map from reference domain [-1, 1]^M to parametric domain [left_coords, right_coords]^M
    // left_coords
    // map from refernce to parameric domain, map_to_parametric, left_coord e right_coord li prende dalla mesh
    Eigen::Matrix<double, MeshType::embed_dim, 1> parametrization(const Eigen::Matrix<double, MeshType::local_dim,1>& p) const {
        return mesh_->eval_param(this->affine_map(p));
    }

    Eigen::Matrix<double, MeshType::embed_dim, MeshType::local_dim> parametrization_gradient(const Eigen::Matrix<double, MeshType::local_dim,1>& p) const {
        return mesh_->eval_param_derivative(this->affine_map(p));
    }

    // Metric tensor F^T * F
    Eigen::Matrix<double, MeshType::local_dim, MeshType::local_dim> metric_tensor(const Eigen::Matrix<double, MeshType::local_dim,1>& p) const {
        auto F = parametrization_gradient(p);
        return F.transpose() * F; 
    }

    // metric determinant sqrt(det(F^T * F)), array diventano matrici eigen
    double metric_determinant(const Eigen::Matrix<double, MeshType::local_dim,1>& p) const {
        return std::sqrt(metric_tensor(p).determinant()); 
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