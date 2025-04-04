#ifndef __FDAPDE_NURBS_DOF_HANDLER_H__
#define __FDAPDE_NURBS_DOF_HANDLER_H__


#include "header_check.h"

namespace fdapde{

template <int LocalDim, int EmbedDim, typename DiscretizationCategory> class DofHandler;

template<int N> class DofHandler<2, N, iso_tag> {

    public:
    using MeshType = IsoMesh<2, N>;
    static constexpr int local_dim = MeshType::local_dim;
    static constexpr int embed_dim = MeshType::embed_dim;

    protected:
    int flatten(const std::array<int, local_dim>& multi_idx) const {
        const auto& dims = mesh_->n_control_points();
        int id = 0;
        int stride = 1;
        for (int d = local_dim - 1; d >= 0; --d) {
            id += multi_idx[d] * stride;
            stride *= dims[d];
        }
        return id;
    }
    std::array<int,local_dim> unflatten(int id) const {
        const auto& dims = mesh_->n_control_points();
        std::array<int,local_dim> multi_idx;
        for (int d = local_dim - 1; d >= 0; --d) {
            multi_idx[d] = id % dims[d];
            id /= dims[d];
        }
        return multi_idx;
    }
    

    public:
    // a geometrical segment with attached dofs
    struct CellType : public IsoSquare<MeshType> {
        using Base = IsoSquare<MeshType>;
        const DofHandler* dof_handler_;
        public:
        static constexpr int local_dim = 2;
        static constexpr int embed_dim = N;

        CellType() : dof_handler_(nullptr) { }
        CellType(int cell_id, const DofHandler* dof_handler):
            Base(cell_id, dof_handler->mesh()), dof_handler_(dof_handler) { } 
            std::vector<int> dofs() const {
                return dof_handler_->active_dofs(Base::id());
            }
            /*
            std::vector<int> dofs_markers() const {
                std::vector<int> dofs_ = dofs();
                std::vector<int> dofs_markers_(dofs_.size());
                for (int i = 0, n = dofs_.size(); i < n; ++i) { dofs_markers_[i] = dof_handler_->dof_marker(dofs_[i]); }
            return dofs_markers_;
            }
            */
           // da aggiustare
            BinaryVector<Dynamic> boundary_dofs() const {
                std::vector<int> dofs_ = dofs();
                BinaryVector<Dynamic> boundary(dofs_.size());
                int i = 0;
                for (int dof : dofs_) {
                    if (dof_handler_->is_dof_on_boundary(dof)) boundary.set(i);
                    ++i;
                }
                return boundary;
            }
    };

    // constructor
    DofHandler() = default;
    DofHandler(const MeshType& mesh) : mesh_(std::addressof(mesh)) {
        order_ = mesh_->order();
        n_dofs_per_cell_ = 1;
        for(int i = 0; i < local_dim; i++) n_dofs_per_cell_ *= order_[i] + 1;
        int n_cells = mesh_->n_cells();
        n_dofs_ = (mesh_->basis()).size();
        const auto& dims = mesh_->n_control_points();

        dofs_.resize(n_cells, n_dofs_per_cell_);

        for (int cell_id = 0; cell_id < n_cells; ++cell_id) {
            auto local_dof_multi_indices = active_dofs(cell_id);  // list of [i,j]
            for (int k = 0; k < local_dof_multi_indices.size(); ++k) {
                dofs_(cell_id, k) = local_dof_multi_indices[k];  // flatten [i,j] → scalar
            }
        }

        // Initialize boundary dofs
        boundary_dofs_.resize(n_dofs_);

        for(int id = 0; id < n_dofs_; id++) {
            auto multi_index = unflatten(id);
            for(int d = 0; d < local_dim; d++) {
                if(multi_index[d] == 0 || multi_index[d] == dims[d] - 1) {
                    boundary_dofs_.set(id);
                    break;
                }
            }
        }

        dofs_markers_ = mesh.nodes_markers();
     }

     // dimension n_dofs_ x local_dim: parametric coordinates of each dof
     Eigen::Matrix<double, Dynamic, local_dim> dof_coords() const{
        Eigen::Matrix<double, Dynamic, local_dim> coords(n_dofs_, local_dim);
        auto nurb = mesh_->basis()[0]; // take a nurb
        std::array<std::vector<double>, local_dim> knot_coords;

        // Extract 1D knot positions for each parametric direction
        for (int i = 0; i < local_dim; i++) {
            auto basis = nurb.spline_basis()[i];
            for (const auto& b : basis) {
                knot_coords[i].push_back(b.knot());
            }
        }

        // Carry-on logic: Cartesian product of knot coordinates
        std::vector<int> idx(local_dim, 0);
        while (true) {
            Eigen::Matrix<double, local_dim, 1> coord;
            for (int d = 0; d < local_dim; ++d) {
                coord(d) = knot_coords[d][idx[d]];
            }
            coords.row(flatten(idx)) = coord;

            // Increment multi-index
            int d = local_dim - 1;
            while (d >= 0) {
                idx[d]++;
                if (idx[d] < static_cast<int>(knot_coords[d].size())) break;
                idx[d] = 0;
                --d;
            }
            if (d < 0) break;
        }

     }
     
 

    // getters
    const MeshType* mesh() const {return mesh_;}
    CellType cell(int id) const { return CellType(id, this); }
    int n_dofs() const { return n_dofs_; }
    int n_dofs_per_cell() const { return n_dofs_per_cell_; }
    bool is_dof_on_boundary(int i) const { return boundary_dofs_[i]; }
    const std::vector<int>& dofs_markers() const { return dofs_markers_; }
    int dof_marker(int dof) const { return dofs_markers_[dof]; }
    int n_boundary_dofs() const { return boundary_dofs_.count(); }
    int n_boundary_dofs(int marker) const {
        int i = 0, sum = 0;
        for (int dof_marker : dofs_markers_) { sum += (dof_marker == marker && boundary_dofs_[i++]) ? 1 : 0; }
        return sum;
    }
    std::vector<int> filter_dofs_by_marker(int marker) const {
        std::vector<int> result;
        for (int i = 0; i < n_dofs_; ++i) {
            if (dofs_markers_[i] == marker) result.push_back(i);
        }
        return result;
    }

    // iterate over geometric cells coupled with dofs, possibly filtered by marker
    class cell_iterator :  public internals::filtering_iterator<cell_iterator, CellType> {
        using Base = internals::filtering_iterator<cell_iterator, CellType>;
        using Base::index_;
        friend Base;
        const DofHandler* dof_handler_;
        int marker_;

        cell_iterator& operator()(int i){
            Base::val_ = dof_handler_->cell(i);
            return *this;
        }

        public:
        cell_iterator() = default;
        cell_iterator(int index, const DofHandler* dof_handler, const BinaryVector<Dynamic>& filter, int marker):
            Base(index, 0, dof_handler->mesh()->n_cells(), filter),
            dof_handler_(dof_handler), 
            marker_(marker){
                for (; index_ < Base::end_ && !filter[index_]; ++index_);
                if (index_ != Base::end_) { operator()(index_); }
            }
            cell_iterator(int index, const DofHandler* dof_handler, int marker) :
            cell_iterator(
              index, dof_handler,
              marker == TriangulationAll ?
                BinaryVector<Dynamic>::Ones(dof_handler->mesh()->n_cells()) :   // apply no filter
                make_binary_vector(
                  dof_handler->mesh()->cells_markers().begin(),
                  dof_handler->mesh()->cells_markers().end(), marker),
              marker) { }
        int marker() const { return marker_; }

        
    }; 

    cell_iterator cells_begin(int marker = TriangulationAll) const {
        const std::vector<int>& cells_markers = mesh_->cells_markers();
        fdapde_assert(marker == TriangulationAll || (marker >= 0 && cells_markers.size() != 0));
        return cell_iterator(0, this, marker);
    }
    cell_iterator cells_end(int marker = TriangulationAll) const {
        fdapde_assert(marker == TriangulationAll || (marker >= 0 && mesh_->cells_markers().size() != 0));
        return cell_iterator(mesh_->n_cells(), this, marker);
    }


    class BoundaryDofType {
        int id_;
        const DofHandler* dof_handler_;
       public:
        BoundaryDofType() = default;
        BoundaryDofType(int id, const DofHandler* dof_handler) : id_(id), dof_handler_(dof_handler) { }
        int id() const { return id_; }
        int marker() const { return dof_handler_->dofs_markers_[id_]; }
        Eigen::Matrix<double, local_dim, 1> coord() const {
            return dof_handler_->dofs_coords_[id_];
        }
    };

    class boundary_dofs_iterator : public internals::filtering_iterator<boundary_dofs_iterator, BoundaryDofType> {
        using Base = internals::filtering_iterator<boundary_dofs_iterator, BoundaryDofType>;
        using Base::index_;
        friend Base;
        const DofHandler* dof_handler_;
        int marker_;
        boundary_dofs_iterator& operator()(int i) {
            Base::val_ = BoundaryDofType(i, dof_handler_);
            return *this;
        }
       public:
        boundary_dofs_iterator(
          int index, const DofHandler* dof_handler, const BinaryVector<Dynamic>& filter, int marker) :
            Base(index, 0, dof_handler->n_dofs(), filter), dof_handler_(dof_handler), marker_(marker) {
            for (; index_ < Base::end_ && !filter[index_]; ++index_);
            if (index_ != Base::end_) { operator()(index_); }
        }
        // filter boundary dofs by marker
        boundary_dofs_iterator(int index, const DofHandler* dof_handler, int marker) :
            boundary_dofs_iterator(
              index, dof_handler,
              marker == BoundaryAll ? dof_handler->boundary_dofs_ :
                                      dof_handler->boundary_dofs_ &
                                        make_binary_vector(
                                          dof_handler->dofs_markers_.begin(), dof_handler->dofs_markers_.end(), marker),
              marker) { }
        int marker() const { return marker_; }
    };
    boundary_dofs_iterator boundary_dofs_begin(int marker = BoundaryAll) const {
        return boundary_dofs_iterator(0, this, marker);
    }
    boundary_dofs_iterator boundary_dofs_end(int marker = BoundaryAll) const {
        return boundary_dofs_iterator(n_dofs_, this, marker);
    }


    // In any given knot span [u_i, u_{i+1}) at most p+1 basis functions are non zero, namely N_{i-p,p}, ..., N_{i,p}
    // (property P2.2, pag 55, Piegl, L., & Tiller, W. (2012). The NURBS book. Springer Science & Business Media.)
    // Evaluation of the non zero basis functions in the a given knot span
    // voglio gli ID, non i punti std::vector<int>
    std::vector<int> active_dofs(int id) const { // id is the cell id
        std::vector<int> dofs;
        auto multi_index = mesh_->cell_multi_index(id);
        std::array<std::vector<double>,local_dim> param_nodes = mesh_->param_nodes();
        Eigen::Matrix<double,local_dim,1> u;
        for(int i = 0; i < local_dim; i++) u(i) = param_nodes[i][multi_index[i]];
        auto nurb = mesh_->basis()[0];
        auto spline_basis = nurb.spline_basis();

        std::vector<std::vector<int>> span_indices(local_dim);
        std::cout << "Cell ID: " << id << std::endl;
        std::cout << "u: " << u.transpose() << std::endl;

        for(int i = 0; i < local_dim; i++) {
            auto basis = spline_basis[i];
            int span = basis->find_span(u(i));
            int p = order_[i];
            for(int j = 0; j <= p; j++) {
                span_indices[i].push_back(span - p + j); 
            }
        }

        // Print the span indices
        std::cout << "Span indices: ";
        for (int i = 0; i < local_dim; ++i) {
            std::cout << "[";
            for (const auto& index : span_indices[i]) {
                std::cout << index << " ";
            }
            std::cout << "] ";
        }
        std::cout << std::endl;

        // Carry-on logic: Cartesian product of all local spans
        std::vector<int> idx(local_dim, 0);
        while (true) {
            std::array<int,local_dim> dof_index;
            for (int d = 0; d < local_dim; ++d) {
                dof_index[d] = span_indices[d][idx[d]];
            }
            dofs.push_back(flatten(dof_index));

            // Increment multi-index
            int d = local_dim - 1;
            while (d >= 0) {
                idx[d]++;
                if (idx[d] < span_indices[d].size()) break;
                idx[d] = 0;
                --d;
            }
            if (d < 0) break;
        }

        // print the active dofs
        std::cout << "Active dofs: ";
        for (int i = 0; i < dofs.size(); i++) {
            std::cout << dofs[i] << " ";
        }
        std::cout << std::endl;

        return dofs;
    }


    private:
    Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor> dofs_; // dofs active on cell: each row = global DOFs on one cell ...
    BinaryVector<Dynamic> boundary_dofs_; // boundary dofs
    int n_dofs_per_cell_ = 0, n_dofs_ = 0;
    std::vector<int> dofs_markers_; // dofs markers
    const MeshType* mesh_;
    std::array<int,local_dim> order_;

};

}


#endif // __FDAPDE_NURBS_DOF_HANDLER_H__