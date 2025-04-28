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

    public:
    int flatten(const std::array<int, local_dim>& multi_idx) const {
        int id = 0;
        int stride = 1;
        for (int d = 0; d < local_dim ; d++) {
            id += multi_idx[d] * stride;
            stride *= dims_[d];
        }
        return id;
    }
    std::array<int,local_dim> unflatten(int id) const {
        std::array<int,local_dim> multi_idx;
        for (int d = 0; d < local_dim ; d++)  {
            multi_idx[d] = id % dims_[d];
            id /= dims_[d];
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
    DofHandler(const MeshType& mesh) : mesh_(std::addressof(mesh)), dof_constraints_(*this) {
        order_ = mesh_->order();
        basis_pde_ = mesh_->basis_pde();
        n_dofs_per_cell_ = 1;
        for(int i = 0; i < local_dim; i++) n_dofs_per_cell_ *= order_[i] + 1;
        int n_cells = mesh_->n_cells();
        n_dofs_ = (basis_pde_).size();
        //const auto& dims = mesh_->n_control_points();
        for (int d = 0; d < local_dim; ++d) {
            dims_[d] = basis_pde_[0].spline_basis()[d]->n_knots() - order_[d] - 1;
            std::cout << "Dimension " << d << ": " << dims_[d] << std::endl;
        }
        

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
                if((multi_index[d] == 0 || multi_index[d] == dims_[d] - 1) && (!mesh_->is_periodic(d))) {
                    boundary_dofs_.set(id);
                    break;
                }
            }
        }

        // for the moment unmarked dofs
        dofs_markers_ = std::vector<int>(n_dofs_, Unmarked);

        //dofs_markers_ = mesh.nodes_markers();
        dof_map_.resize(n_dofs_);
        for (int i = 0; i < n_dofs_; ++i) {
            dof_map_[i] = i;
        }


        build_periodic_dof_map2();
        /*
        // print the periodic dof map
        std::cout << "Periodic DOF map: " << std::endl;
        for (int i = 0; i < n_dofs_; ++i) {
            std::cout << "DOF " << i << " maps to " << dof_map_[i] << std::endl;

         }
            */


         int next_index = 0;

        for (int i = 0; i < dof_map_.size(); ++i) {
            int mapped = dof_map_[i];
            if (compressed_map_.find(mapped) == compressed_map_.end()) {
                compressed_map_[mapped] = next_index++;
            }
        }

        reduced_dof_map_.resize(n_dofs_);
        for (int i = 0; i < n_dofs_; ++i) {
            reduced_dof_map_[i] = compressed_map_[dof_map_[i]];
        }

        // Print the reduced dof map
        //std::cout << "Reduced DOF map: " << std::endl;
        //for (int i = 0; i <n_dofs_; ++i) {
          //  std::cout << "DOF " << i << " maps to " << reduced_dof_map_[i] << std::endl;
        //}
        // Print the number of mapped dofs
        //std::cout << "Number of mapped DOFs: " << n_mapped_dofs_ << std::endl;
        // Print the number of dofs


    }

    

     template <typename SystemMatrix, typename SystemRhs>
     void enforce_constraints(SystemMatrix&& A, SystemRhs&& b) const {
         dof_constraints_.enforce_constraints(std::forward<SystemMatrix>(A), std::forward<SystemRhs>(b));
     }

     void enforce_constraints(Eigen::SparseMatrix<double>& A) const {
        dof_constraints_.enforce_constraints(A);
    }

    void enforce_constraints(Eigen::Matrix<double, Dynamic, 1>& b) const {
        dof_constraints_.enforce_constraints(b);
    }

    void set_hom_dirichlet_constraint(int marker = BoundaryAll) {
        dof_constraints_.set_hom_dirichlet_constraint(marker);

        
    }

    void set_periodic_constraint() { // non so se funziona, only C0 continuity
        std::set<int> constrained_dofs;
        std::set<int> master_set;
        for (int d = 0; d < local_dim; ++d) {
            if (!mesh_->is_periodic(d)) continue;
            std::cout << "Periodic constraint in dimension " << d << std::endl;
        
            int p = order_[d]; //order_[p];      // spline degree
            int n = dims_[d];       // number of DOFs in this direction
        
            std::array<int, local_dim> idx_min, idx_max;
            
            // Loop over all combinations in other dimensions
            std::vector<int> dofs_min, dofs_max;
            std::array<int, local_dim> multi;
            int master_north_pole = -1;
            int master_south_pole = -1;
            for (int i = 0; i < n_dofs_; ++i) {
                
                // For a sphere: v = 0 is south pole, v = dims_[1]-1 is north pole
                if (unflatten(i)[1] == dims_[1] - 1 && master_north_pole == -1) {
                    master_north_pole = i;
                    master_set.insert(i);
                    continue;  // Skip this DoF from periodic constraint setup
                }
                if (unflatten(i)[1] == 0 && master_south_pole == -1) {
                    master_south_pole = i;
                    master_set.insert(i);
                    continue;  // Same here
                }
                    

                multi = unflatten(i);
                if (multi[d] < order_[d]) {
                    auto mapped = multi;
                    mapped[d] += n - order_[d];
                    int mapped_id = flatten(mapped);
                    dofs_min.push_back(i);
                    dofs_max.push_back(mapped_id);
                    master_set.insert(i);
                }
                // TEMPORARY: only for sphere: constraint on the dofs at the poles, put the other dofs equal to the first one
            }

            // print dofs_min and dofs_max
            for(int i = 0; i < dofs_min.size(); ++i) {
                std::cout << "dofs_min: " << dofs_min[i] << " dofs_max: " << dofs_max[i] << std::endl;
            }
            
            
            std::cout<<"Master north pole: "<<master_north_pole<<std::endl;
            std::cout<<"Master south pole: "<<master_south_pole<<std::endl;
            
            
            for(int i = 0; i < n_dofs_; ++i) {
                // v = dims_[1] - 1 is north pole, v = 0 is south pole
                if(unflatten(i)[1] == dims_[1] - 1 && i != master_north_pole) { // north pole
                    dof_constraints_.set_master_slave_constraint(master_north_pole, i);
                    constrained_dofs.insert(i);
                }
                if(unflatten(i)[1] == 0 && i != master_south_pole) { // south pole
                    dof_constraints_.set_master_slave_constraint(master_south_pole, i);
                    constrained_dofs.insert(i);
                }
            }
                
                
                
                
            //dof_constraints_.set_master_slave_constraint(master_south_pole, master_north_pole);

            std::cout << "number of Periodic constraint: " << dofs_min.size() << " dofs" << std::endl;
        
            // Enforce the periodic constraint: max ↔ min
            for (int i = 0; i < dofs_min.size(); ++i) {
                int master = dofs_min[i];
                int slave  = dofs_max[i];
                std::cout << "Periodic constraint: " << master << " ↔ " << slave << std::endl;
                for(int j = 0; j < local_dim; ++j) {
                    std::cout << "multi: " << unflatten(master)[j] << " ↔ " << unflatten(slave)[j] << std::endl;
                }
                if(constrained_dofs.find(slave) != constrained_dofs.end()) {
                    std::cout << "Skipped constraint: " << slave << " already constrained.\n";
                }
                if(constrained_dofs.find(slave) == constrained_dofs.end()) {
                    dof_constraints_.set_master_slave_constraint(master, slave);
                    constrained_dofs.insert(slave);
                }
                //dof_constraints_.set_master_slave_constraint(master, slave); // enforce slave = master
            }


            
        }
    }
    void set_periodic_constraint2() { // C0 periodic continuity
        std::set<int> constrained_dofs;

        // Handle corner wrap (only once, when all dimensions are periodic)
        bool all_periodic = true;
        for (int d = 0; d < local_dim; ++d) {
            if (!mesh_->is_periodic(d)) {
                all_periodic = false;
                break;
            }
        }

        if (all_periodic) {
            int p0 = order_[0];
            int p1 = order_[1];
            int n0 = dims_[0];
            int n1 = dims_[1];
            int master = flatten({0, 0});
            for (int i = 0; i < 1; ++i) {
                for (int j = 0; j < 1; ++j) {
                    if (i == 0 && j == 0) continue;
                    int slave = flatten({i + n0 - p0, j + n1 - p1});
                    if (constrained_dofs.find(slave) == constrained_dofs.end()) {
                        dof_constraints_.set_master_slave_constraint(master, slave);
                        constrained_dofs.insert(slave);
                    }
                }
            }

            // Edge wrap along top (in j)
            for (int j = 0; j < p1; ++j) {
                for (int i = 0; i < n0; ++i) {
                    int master = flatten({i, j});
                    int slave = flatten({i, j + n1 - p1});
                    if (constrained_dofs.find(slave) == constrained_dofs.end()) {
                        dof_constraints_.set_master_slave_constraint(master, slave);
                        constrained_dofs.insert(slave);
                    }
                }
            }

            // Edge wrap along right (in i)
            for (int i = 0; i < p0; ++i) {
                for (int j = 0; j < n1; ++j) {
                    int master = flatten({i, j});
                    int slave = flatten({i + n0 - p0, j});
                    if (constrained_dofs.find(slave) == constrained_dofs.end()) {
                        dof_constraints_.set_master_slave_constraint(master, slave);
                        constrained_dofs.insert(slave);
                    }
                }
            }
        }

        for (int d = 0; d < local_dim; ++d) {
            if (!mesh_->is_periodic(d)) continue;
            std::cout << "Periodic constraint in dimension " << d << std::endl;

            int p = order_[d];
            int n = dims_[d];

            // Edge wrapping (1D per dimension)
            for (int i = 0; i < n_dofs_; ++i) {
                auto idx = unflatten(i);
                if (idx[d] < p) {
                    auto mapped = idx;
                    mapped[d] = idx[d] + (n - p);
                    int master = flatten(idx);
                    int slave = flatten(mapped);
                    if (constrained_dofs.find(slave) != constrained_dofs.end()) {
                        std::cout << "Skipped constraint: " << slave << " already constrained.\n";
                    }
                    if (constrained_dofs.find(slave) == constrained_dofs.end()) {
                        dof_constraints_.set_master_slave_constraint(master, slave);
                        constrained_dofs.insert(slave);
                    }
                }
            }
        }


    }



    void set_periodic_constraint3() {
        std::set<int> constrained_dofs;
        for (int i = 0; i < n_dofs_; ++i) {
            int mapped = dof_map_[i];
            if (mapped != i && constrained_dofs.find(i) == constrained_dofs.end()) {
                std::cout << "Periodic constraint: " << i << " ↔ " << mapped << std::endl;
                dof_constraints_.set_master_slave_constraint(mapped, i);
                constrained_dofs.insert(i);
            }
        }
    }
    
    
    
    void get_boundary_dofs_for_dimension(int dim, bool min_side, std::vector<int>& dofs) const {
        for (int i = 0; i < n_dofs_; ++i) {
            if(boundary_dofs_[i]) {
                auto multi_index = unflatten(i);
                if ((min_side && multi_index[dim] == 0) || (!min_side && multi_index[dim] == dims_[dim] - 1)) {
                    dofs.push_back(i);
                }
            }
        }
    }


     // dimension n_dofs_ x local_dim: parametric coordinates of each dof
     Eigen::Matrix<double, Dynamic, local_dim> dof_coords() const{
        Eigen::Matrix<double, Dynamic, local_dim> coords(n_dofs_, local_dim);
        auto nurb = basis_pde_[0]; // take a nurb
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
     
    void build_periodic_dof_map() {
        // Initialize all entries with identity mapping
        dof_map_.resize(n_dofs_);
        std::iota(dof_map_.begin(), dof_map_.end(), 0);  // identity map

        for (int d = 0; d < local_dim; ++d) {
            if (!mesh_->is_periodic(d)) continue;

            int n = dims_[d];

            for (int i = 0; i < n_dofs_; ++i) {
                auto idx = unflatten(i);
                if (idx[d] == n - 1) {
                    auto wrapped = idx;
                    wrapped[d] = 0;
                    int target = flatten(wrapped);
                    dof_map_[i] = target;
                }
            }
        }
        

        // Handle poles
        int master_south = -1, master_north = -1;
        for (int i = 0; i < n_dofs_; ++i) {
            auto idx = unflatten(i);
            if (idx[1] == 0) {
                if (master_south == -1) master_south = i;
                dof_map_[i] = master_south;
            } else if (idx[1] == dims_[1] - 1) {
                if (master_north == -1) master_north = i;
                dof_map_[i] = master_north;
            }
        }
            

        // print the dof_map
        std::cout << "DOF map: " << std::endl;
        for (int i = 0; i < n_dofs_; ++i) {
            std::cout << "DOF " << i << " maps to " << dof_map_[i] << std::endl;
        }

        std::unordered_set<int> unique_dofs;
        for (int i = 0; i < n_dofs_; ++i)
            unique_dofs.insert(dof_map_[i]);

        n_mapped_dofs_ = unique_dofs.size();
    }
    
    
    void build_periodic_dof_map2() {
        for (int d = 0; d < local_dim; ++d) {
            if (!mesh_->is_periodic(d)) continue;
            int n = dims_[d];

            for (int i = 0; i < n_dofs_; ++i) {
                auto multi = unflatten(i);
                if (multi[d] < order_[d]) {
                    auto mapped = multi;
                    mapped[d] += n - order_[d];
                    dof_map_[i] = flatten(mapped);
                }
            }
        }
        if (mesh_->is_periodic(0) && mesh_->is_periodic(1)) {
            int p0 = order_[0];
            int p1 = order_[1];
            int n0 = dims_[0];
            int n1 = dims_[1];
            
            for (int i = 0; i < p0; ++i) {
                for (int j = 0; j < p1; ++j) {
                    int master = flatten({i, j});
                    int slave  = flatten({i + n0 - p0, j + n1 - p1});
                    dof_map_[slave] = master;
                }
            }
        }


        
        // Handle poles (collapse DOFs at the north and south poles)
        /*
        int master_south = -1, master_north = -1;
        for (int i = 0; i < n_dofs_; ++i) {
            auto idx = unflatten(i);
            if (idx[1] == 0) {
                if (master_south == -1) master_south = i;
                dof_map_[i] = master_south;
            } else if (idx[1] == dims_[1] - 1) {
                if (master_north == -1) master_north = i;
                dof_map_[i] = master_north;
            }
        }
        */
            
            

        std::unordered_set<int> unique_dofs;
        for (int i = 0; i < n_dofs_; ++i)
            unique_dofs.insert(dof_map_[i]);
        n_mapped_dofs_ = unique_dofs.size();



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

    int n_mapped_dofs() const { return n_mapped_dofs_; }
    const std::vector<int>& reduced_dof_map() const { return reduced_dof_map_; }
    std::array<int, local_dim> dims() const { return dims_; }

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
        BoundaryDofType(int id, const DofHandler* dof_handler) : id_(id), dof_handler_(dof_handler) { 
        }
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
          Base(index, 0, dof_handler->n_dofs(), filter), dof_handler_(dof_handler), marker_(marker)
             {
            for (; index_ < Base::end_ && !filter[index_]; ++index_);
            if (index_ != Base::end_) { operator()(index_); }
        }
        // filter boundary dofs by marker
        boundary_dofs_iterator(int index, const DofHandler* dof_handler, int marker) :
            boundary_dofs_iterator(
              index, dof_handler,
              dof_handler->boundary_dofs_,
              marker) { 
              }
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
        auto nurb = basis_pde_[0];
        auto spline_basis = nurb.spline_basis();

        std::vector<std::vector<int>> span_indices(local_dim);
        //std::cout << "Cell ID: " << id << std::endl;
        //std::cout << "u: " << u.transpose() << std::endl;

        for(int i = 0; i < local_dim; i++) {
            auto basis = spline_basis[i];
            int span = basis->find_span(u(i));
            int p = order_[i];
            for(int j = 0; j <= p; j++) {
                span_indices[i].push_back(span - p + j); 
            }
        }
        /*
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
        */

        // Carry-on logic: Cartesian product of all local spans
        std::vector<int> idx(local_dim, 0);
        while (true) {
            std::array<int,local_dim> dof_index;
            for (int d = 0; d < local_dim; ++d) {
                dof_index[d] = span_indices[d][idx[d]];
            }
            /*
            std::cout << "Inserting dof_index: ";
            for (int d = 0; d < local_dim; ++d) {
                std::cout << dof_index[d] << " ";
            }
            std::cout<<"as dof: "<<flatten(dof_index) << std::endl;
            */
            
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
        /*
        // print the active dofs
        std::cout << "Active dofs: ";
        for (int i = 0; i < dofs.size(); i++) {
            std::cout << dofs[i] << " ";
        }
        std::cout << std::endl;
        */

        return dofs;
    }

    std::vector<int> dof_map() const {
        return dof_map_;
    }


    private:
    std::vector<int> dof_map_;  // dof_map_[original_dof] = reduced_dof
    int n_mapped_dofs_;         // Number of unique DOFs after collapsing periodic ones
    std::unordered_map<int, int> compressed_map_; // map old index → compressed index
    std::vector<int> reduced_dof_map_; // reduced_dof_map_[compressed_dof] = original_dof
    
    
    NurbsBasis<local_dim> basis_pde_;
    std::array<int,local_dim> dims_; // number of control points in each direction
    Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor> dofs_; // dofs active on cell: each row = global DOFs on one cell ...
    BinaryVector<Dynamic> boundary_dofs_; // boundary dofs
    int n_dofs_per_cell_ = 0, n_dofs_ = 0;
    std::vector<int> dofs_markers_; // dofs markers
    const MeshType* mesh_;
    std::array<int,local_dim> order_;

    DofConstraints<DofHandler> dof_constraints_;

};

}


#endif // __FDAPDE_NURBS_DOF_HANDLER_H__