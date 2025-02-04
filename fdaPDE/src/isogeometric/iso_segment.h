#ifndef __FDAPDE_ISO_SEGMENT_H__
#define __FDAPDE_ISO_SEGMENT_H__

#include "header_check.h"

namespace fdapde {

template <typename MeshType> class IsoSegment: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 1, THIS_CLASS_IS_FOR_INTERVAL_MESHES_ONLY);
    public:
    // constructor
    IsoSegment() = default;
    IsoSegment(int id, const MeshType* mesh) : id_(id), mesh_(mesh), boundary_(false) {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        left_coords_ = mesh_->compute_lr_vertices_(id_)[0];
        right_coords_ = mesh_->compute_lr_vertices_(id_)[1];
        // initialize = (){}; // da capire cosa inizializzare
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
}
    
    
}; // namespace fdapde

#endif // __FDAPDE_ISO_SEGMENT_H__