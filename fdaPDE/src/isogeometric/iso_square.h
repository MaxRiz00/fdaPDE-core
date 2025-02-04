#ifndef __FDAPDE_ISO_SEGMENT_H__
#define __FDAPDE_ISO_SEGMENT_H__

#include "header_check.h"

namespace fdapde {

template <typename MeshType> class IsoSegment: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 2, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    using Base = IsoCell<MeshType::local_dim, MeshType::embed_dim>;
    public:
    // constructor
    IsoSegment() = default;
    IsoSegment(int id, const MeshType* mesh) : id_(id), mesh_(mesh), boundary_(false) {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        left_coords_ = mesh_->compute_lr_vertices_(id_)[0];
        right_coords_ = mesh_->compute_lr_vertices_(id_)[1];
        // initialize = (){}; // da capire cosa inizializzare
    }

    class EdgeType : public IsoCell<MeshType::local_dim, MeshType::embed_dim>::BoundaryCellType {
        private:
        int edge_id_;
        const MeshType* mesh_;
        public:
        using CoordsType = Eigen::Matrix<double, MeshType::embed_dim, MeshType::local_dim>; 
        EdgeType() = default; //// da capire come adattare!!!
        EdgeType(int edge_id, const MeshType* mesh) : edge_id_(edge_id), mesh_(mesh) {
            for (int i = 0; i < this->n_nodes; ++i) { this->coords_.col(i) = mesh_->node(mesh_->edges()(edge_id_, i)); }
            this->initialize();
        }
    };


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
}
    
    
}; // namespace fdapde

#endif // __FDAPDE_ISO_SEGMENT_H__