#ifndef __ISO_MESH_H__
#define __ISO_MESH_H__


#include "nurbs.h"
#include "iso_square.h"


namespace fdapde {

// public iterator types, problema per la doppia definizion in triangulation.h CHIEDERE
/*
template <typename MeshType> struct CellIterator : public MeshType::cell_iterator {
    CellIterator(int index, const MeshType* mesh) : MeshType::cell_iterator(index, mesh) { }
    CellIterator(int index, const MeshType* mesh, int marker) :
        MeshType::cell_iterator(index, mesh, marker) { }
}; */
/*
template <typename MeshType> struct BoundaryIterator : public MeshType::boundary_iterator {
    BoundaryIterator(int index, const MeshType* mesh) : MeshType::boundary_iterator(index, mesh) { }
    BoundaryIterator(int index, const MeshType* mesh, int marker) :
        MeshType::boundary_iterator(index, mesh, marker) { }
};*/

template <int LocalDim, int EmbedDim> class IsoMesh;

template <int LocalDim, int EmbedDim, typename Derived> class IsoMeshBase{ //typname Derived
    public:
    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;
    static constexpr int n_nodes_per_cell = 1<<LocalDim;
    static constexpr int n_neighbors_per_cell = 2*LocalDim;
    static constexpr bool is_manifold = !(local_dim == embed_dim);

    using CellType = std::conditional_t<
      local_dim == 1, IsoSegment<Derived>, std::conditional_t<local_dim == 2, IsoSquare<Derived>, IsoCube<Derived>>>;
    using MeshType = Derived;
    

    class NodeType {
        int id_;
        const MeshType* mesh_;
        public:
        NodeType() = default;
        NodeType(int id, const MeshType* mesh) : id_(id), mesh_(mesh) { }
        int id() const { return id_; }
        //Eigen::Matrix<double, embed_dim, 1> coords() const { return mesh_->node(id_); }
        //std::vector<int> patch() const { return mesh_->node_patch(id_); }         // cells having this node as vertex
        //std::vector<int> one_ring() const { return mesh_->node_one_ring(id_); }   // directly connected nodes
    };
    

    IsoMeshBase() = default;

    IsoMeshBase(std::array<std::vector<double>,LocalDim> & knots, MdArray<double,full_dynamic_extent_t<LocalDim>> & weights, 
         MdArray<double,full_dynamic_extent_t<LocalDim+1>> & control_points, int order, int flags=0): knots_(knots), weights_(weights), control_points_(control_points), flags_(flags){
            basis_ = NurbsBasis<LocalDim>(knots, weights, order);
            // compute number of cells and nodes
            n_cells_ = 1; // is_cell cambialo in cells, ovunque
            n_nodes_ = 1; //
            for(int i=0; i < LocalDim; ++i){
                n_cells_ *= knots_[i].size() - 1;
                n_nodes_ *= knots_[i].size();
            }
            // compute neighbors and cells, da capire se tenerli
            //neighbors_ = neighbors();
            //cells_ = cells();
            //param_nodes_ = parametric_nodes();

        };


    // getters
    const NurbsBasis<LocalDim>& basis() const { return basis_; }
    const MdArray<double, full_dynamic_extent_t<LocalDim+1>>& control_points() const { return control_points_; }
    const std::array<std::vector<double>,LocalDim>& knots() const { return knots_; }
    int n_cells() const { return n_cells_; }
    int n_nodes() const { return n_nodes_; }

    protected:
    int compute_stride_(int dim, bool is_cell) const {
        return is_cell ? this->knots_[dim].size() - 1 : this->knots_[dim].size();   
    }

    // Algo A4.3 from NURBS book pag. 103, evaluation of a NURBS curve
    Eigen::Matrix<double, EmbedDim, 1> eval_param_(const std::array<double, LocalDim>& u) const {

        std::vector<std::vector<double>> basis_eval(LocalDim);
        std::array<int,LocalDim> spans= {0};
        auto order= this->basis_[0].order();
        auto nurb = this->basis_[0];
        double total_weight = 0.0;
        
        for(int i = 0; i < LocalDim; i++){
            basis_eval[i] = nurb.spline_basis()[i]->evaluate_basis(u[i], false); // evaluate basis functions, padding = false
            spans[i] = std::distance(knots_[i].begin(), std::upper_bound(knots_[i].begin(), knots_[i].end(), u[i])) - 1 + order;
        }

        Eigen::Matrix<double, EmbedDim, 1> Sw = Eigen::Matrix<double, EmbedDim, 1>::Zero();

        std::vector<int> index(LocalDim,0);
        bool done =  false;
        while(!done){
            double eval = 1.0;

            std::array<int,LocalDim> full_indices;
            for(int i = 0; i < LocalDim; i++){
                eval *= basis_eval[i][index[i]];
                full_indices[i] = spans[i] - order + index[i];
            }

            Eigen::Matrix<double, EmbedDim, 1> cp;
            for(int i = 0; i < EmbedDim; i++){
                const auto cp_slice = this->control_points_.template slice<LocalDim>(i);
                cp(i) = cp_slice(full_indices);
            }

            double w = weights_(index);  // Retrieve weight from weight array
            cp *= w;  // Apply weight to control point
            Sw += eval * cp; // Accumulate weighted control point
            
            total_weight += eval * w ;  // Accumulate total weight

            for (int d = LocalDim - 1; d >= 0; d--) {
            if (++index[d] > order) {
                index[d] = 0;
                if (d == 0) done = true;
            } else 
                break;
            }
        }
        return Sw/total_weight;
    }

    // Algo A4.3 from NURBS book pag. 103, evaluation of the physical derivative a NURBS curve
    Eigen::Matrix<double, EmbedDim, LocalDim, Eigen::RowMajor> eval_param_derivative_(const std::array<double, LocalDim>& u) const {
        std::vector<std::vector<double>> basis_eval(LocalDim);
        std::vector<std::vector<double>> basis_deriv_eval(LocalDim);
        std::array<int, LocalDim> spans = {0};
        auto order = this->basis_[0].order();
        auto nurb = this->basis_[0];
        
        for (int i = 0; i < LocalDim; i++) {
            basis_eval[i] = nurb.spline_basis()[i]->evaluate_basis(u[i], false);
            basis_deriv_eval[i] = nurb.spline_basis()[i]->evaluate_der_basis(u[i],1,false);
            spans[i] = std::distance(knots_[i].begin(), std::upper_bound(knots_[i].begin(), knots_[i].end(), u[i])) - 1 + order;
        }
        
        Eigen::Matrix<double, EmbedDim, LocalDim, Eigen::RowMajor> dSw = Eigen::Matrix<double, EmbedDim, LocalDim>::Zero();
        Eigen::Matrix<double, EmbedDim, 1> Sw = Eigen::Matrix<double, EmbedDim, 1>::Zero();
        Eigen::Matrix<double, LocalDim, 1> dW = Eigen::Matrix<double, LocalDim, 1>::Zero();
        double total_weight = 0.0;
        
        std::vector<int> index(LocalDim, 0);
        bool done = false;
        while (!done) {
            double eval = 1.0;
            std::array<double, LocalDim> eval_der = {0.0};        
            std::array<int, LocalDim> full_indices;

            for (int i = 0; i < LocalDim; i++) {
                eval *= basis_eval[i][index[i]];
                eval_der[i] = basis_deriv_eval[i][index[i]];
                full_indices[i] = spans[i] - order + index[i];
            }

            
            Eigen::Matrix<double, EmbedDim, 1> cp;
            for (int i = 0; i < EmbedDim; i++) {
                const auto cp_slice = this->control_points_.template slice<LocalDim>(i);
                cp(i) = cp_slice(full_indices);
            }
            
            double w = weights_(index);
            cp *= w;
            Sw += eval * cp;
            total_weight += eval * w;
            
            for (int j = 0; j < LocalDim; j++) {
                double w_temp = w*eval_der[j];
                Eigen::Matrix<double, EmbedDim, 1> cp_temp = eval_der[j] *cp;
                for(int i = 0; i < LocalDim; i++){
                    if(i != j){
                        w_temp *= basis_eval[i][index[i]];
                        cp_temp *= basis_eval[i][index[i]];
                    }
                }
                dW(j) += w_temp;
                dSw.col(j) += cp_temp;
            }
            
            for (int d = LocalDim - 1; d >= 0; d--) {
                if (++index[d] > order) {
                    index[d] = 0;
                    if (d == 0) done = true;
                } else 
                    break;
            }
        }
        
        // Apply derivative of the ratio d/du (Sw / total_weight)
        for (int j = 0; j < LocalDim; j++) {
            dSw.col(j) = (dSw.col(j) - (Sw * dW(j) / total_weight)) / total_weight;
        }
        
        return dSw;
    }  

    // Compute the id of a cell (or a node if is_cell is false) from the multi-index
    int compute_id_(const std::array<int, LocalDim>& multi_index, bool is_cell=true) const {
        int id = multi_index[LocalDim - 1]; // Start with the last index
        for (int i = LocalDim - 1; i > 0; --i) { // Iterate backward
            id = id * compute_stride_(i-1,is_cell) + multi_index[i - 1];
        }
        return id;
    }

    // Compute the multi-index of a cell (or a node if is_cell is false) from the id
    std::array<int, LocalDim> compute_multi_index_(const int& id, bool is_cell=true) const {
        auto tempID = id;
        std::array<int, LocalDim> multi_index;
        for (int i = LocalDim; i > 0; --i) {
            if (i == 1) {
                multi_index[LocalDim-1] = tempID; // The last index is the remaining id
            } else {
                int stride = compute_stride_(i - 2, is_cell); // Use (i-2) to match dimensions
                multi_index[LocalDim - i] = tempID % stride; // Extract the current index
                tempID /= stride; // Update id for the next dimension
            }
        }
        return multi_index;
    }
    public:
    // Compute the vertices of a square (cell) given its id, fai MdArray , sono metodi pirvati
    std::array<std::array<int, LocalDim>,2> compute_lr_vertices_(const int& id)const {
        auto multi_index = compute_multi_index_(id);
        std::array<std::array<int, LocalDim>,2> vertices;
        for(int i = 0; i < LocalDim; ++i){
            vertices[0][i] = multi_index[i];
            vertices[1][i] = multi_index[i] + 1;
        }
        return vertices;
    }

    public:

    Eigen::Matrix<double, Dynamic, LocalDim, Eigen::RowMajor> parametric_nodes() const {
        Eigen::Matrix<double, Dynamic, LocalDim> nodes;
        nodes.resize(n_nodes_, LocalDim);
        for (int i = 0; i < n_nodes_; ++i) {
            auto multi_index = compute_multi_index_(i, false);
            for (int j = 0; j < LocalDim; ++j) {
                nodes(i, j) = knots_[j][multi_index[j]];
            }
        }
        return nodes;
    }

    // Compute all cells of the mesh (dimension #cells x LocalDim ) // is_cell diventa cells, Eigen::Matrix<int, LocalDim, Dynamic>
    Eigen::Matrix<int, Dynamic, LocalDim, Eigen::RowMajor> cells() const {
        Eigen::Matrix<int, Dynamic, LocalDim> cells;
        cells.resize(n_cells_,LocalDim);
        for(int i = 0; i < n_cells_; ++i){
            cells[i] = compute_multi_index_(i);
        }
        return cells;
    }
    
    // check if a cell/node is on the boundary, usa is_node_on_boundary, is_cell_on_boundary
    bool is_node_on_boundary(const int& id) const {
        auto multi_index = compute_multi_index_(id,false);
        for(int i = 0; i < LocalDim; ++i){
            if(multi_index[i] == 0 || multi_index[i] == compute_stride_(i,false) - 1){
                return true;
            }
        }
        return false;
    }
    bool is_cell_on_boundary(const int& id) const {
        auto multi_index = compute_multi_index_(id);
        for(int i = 0; i < LocalDim; ++i){
            if(multi_index[i] == 0 || multi_index[i] == compute_stride_(i,true) - 1){
                return true;
            }
        }
        return false;
    }

    //guarda triangulation.h per vedere come fare, riga 83
    
    // Compute the id neighbors of a cell given its id (dimension #cells x 2M(at most))
    //   - Colonne `2*j` : Voisin dans la direction négative pour la dimension `j`.
    //   - Colonne `2*j + 1` : Voisin dans la direction positive pour la dimension `j`.
    // - Si un voisin n'existe pas dans une direction, la valeur est -1.
    Eigen::Matrix<int, Dynamic, 2 * LocalDim, Eigen::RowMajor> neighbors() const {
        // Initialize the neighbors matrix with dynamic rows and 2 * LocalDim columns
        Eigen::Matrix<int, Eigen::Dynamic, 2 * LocalDim, Eigen::RowMajor> neighbors;
        neighbors.resize(n_cells_, 2 * LocalDim);

        // Loop through each cell
        for (int i = 0; i < n_cells_; ++i) {
            // Compute the multi-index for the current cell
            auto multi_index = compute_multi_index_(i);

            // Loop through each dimension
            for (int j = 0; j < LocalDim; ++j) {
                // Handle the neighbor in the negative direction
                if (multi_index[j] > 0) {
                    auto updated_index = multi_index;
                    updated_index[j] -= 1; // Update the index in the negative direction
                    neighbors(i, 2 * j) = compute_id_(updated_index);
                } else {
                    neighbors(i, 2 * j) = -1; // No neighbor in the negative direction
                }

                // Handle the neighbor in the positive direction
                if (multi_index[j] < compute_stride_(j, true) - 1) {
                    auto updated_index = multi_index;
                    updated_index[j] += 1; // Update the index in the positive direction
                    neighbors(i, 2 * j + 1) = compute_id_(updated_index);
                } else {
                    neighbors(i, 2 * j + 1) = -1; // No neighbor in the positive direction
                }
            }
        }
        // order the neighbors rows in increasing order, but put the -1 at the end
        // ONLY FOR THE TESTS
        /*
        for(int i = 0; i < n_cells_; ++i){
            std::vector<int> row(neighbors.row(i).data(), neighbors.row(i).data() + neighbors.cols());
            std::sort(row.begin(), row.end(), [](int a, int b) { return a == -1 ? false : b == -1 ? true : a < b; });
            std::copy(row.begin(), row.end(), neighbors.row(i).data());
        }
        */

        return neighbors;
    }
    
    // capire se mettere una parametric_measure

    // Access individual cells, Create a cell object given its id
    // cell(id) -> IsoSquare
    // da cambiare una volta che hai IsoSegment, IsoSquare, IsoCube

    // Make eval_param and eval_param_derivative public
    Eigen::Matrix<double, EmbedDim, 1> eval_param(const std::array<double, LocalDim>& u) const {
        return eval_param_(u);
    }

    Eigen::Matrix<double, EmbedDim, LocalDim, Eigen::RowMajor> eval_param_derivative(const std::array<double, LocalDim>& u) const {
        return eval_param_derivative_(u);
    }
    

    class cell_iterator: public internals::filtering_iterator<cell_iterator,const CellType*>   { 
    // cell_iterator, guarda triangulation.h al rigo 96, uniforma a quella implementazione
    // public internals::filtering_iterator<cell_iterator, const CellType*>
    // call operator a riga 102
    private:
        using Base = internals::filtering_iterator<cell_iterator, const CellType*>;
        using Base::index_;
        friend Base;
        const Derived* mesh_;
        int marker_;

        mutable CellType cell_;  // Store the actual cell (mutable for const correctness)

        cell_iterator& operator()(int i) {
            cell_ = mesh_->cell(i);  // Store the value in a variable
            Base::val_ = &cell_;  // Take the address of the variable
            return *this;
        }

    public:
        using MeshType = Derived;
        cell_iterator() = default;

        cell_iterator(int index, const Derived* mesh, const BinaryVector<Dynamic>& filter, int marker):
            Base(index, 0, mesh->n_cells_, filter), mesh_(mesh), marker_(marker) {
            for (; index_ < Base::end_ && !filter[index_]; ++index_);
            if (index_ != Base::end_) { operator()(index_); }
        }

        cell_iterator(int index, const Derived* mesh, int marker):
            cell_iterator(index, mesh, marker == TriangulationAll ? BinaryVector<Dynamic>::Ones(mesh->n_cells_) :
                make_binary_vector(mesh->cells_markers_.begin(), mesh->cells_markers_.end(), marker), marker) { }

        int marker() const { return marker_; }

    };

    CellIterator<Derived> cells_begin(int marker = TriangulationAll) const {
        fdapde_assert(marker == TriangulationAll || (marker >= 0 && cells_markers_.size() != 0));
        return CellIterator<Derived>(0, static_cast<const Derived*>(this), marker);
    }
    CellIterator<Derived> cells_end(int marker = TriangulationAll) const {
        fdapde_assert(marker == TriangulationAll || (marker >= 0 && cells_markers_.size() != 0));
        return CellIterator<Derived>(n_cells_, static_cast<const Derived*>(this), marker);
    }

    // set cells markers
    template <typename Lambda> void mark_cells(int marker, Lambda&& lambda)
        requires(requires(Lambda lambda, CellType c) {
            { lambda(c) } -> std::same_as<bool>;
        }) {
        fdapde_assert(marker >= 0);
        cells_markers_.resize(n_cells_);
        for (cell_iterator it = cells_begin(); it != cells_end(); ++it) {
            cells_markers_[it->id()] = lambda(*it) ? marker : Unmarked;
        }
    }
    template <int Rows, typename XprType> void mark_cells(const BinMtxBase<Rows, 1, XprType>& mask) {
        fdapde_assert(mask.rows() == n_cells_);
        cells_markers_.resize(n_cells_);
        for (cell_iterator it = cells_begin(); it != cells_end(); ++it) {
            cells_markers_[it->id()] = mask[it->id()] ? 1 : 0;
        }
    }

    template <typename Iterator> void mark_cells(Iterator first, Iterator last) {
        fdapde_static_assert(
          std::is_convertible_v<typename Iterator::value_type FDAPDE_COMMA int>, INVALID_ITERATOR_RANGE);
        int n_markers = std::distance(first, last);
        bool all_markers_positive = std::all_of(first, last, [](auto marker) { return marker >= 0; });
        fdapde_assert(n_markers == n_cells_ && all_markers_positive);
        cells_markers_.resize(n_cells_, Unmarked);
        for (int i = 0; i < n_cells_; ++i) { cells_markers_[i] = *(first + i); }
    }
    void mark_cells(int marker) {   // marks all cells with m
        fdapde_assert(marker >= 0);
        cells_markers_.resize(n_cells_);
	std::for_each(cells_markers_.begin(), cells_markers_.end(), [marker](int& marker_) { marker_ = marker; });
    }
    void clear_cell_markers() {
        std::for_each(cells_markers_.begin(), cells_markers_.end(), [](int& marker) { marker = Unmarked; });
    }
    
    protected:
        std::array<std::vector<double>,LocalDim> knots_; // knots in each direction
        MdArray<double,full_dynamic_extent_t<LocalDim>> weights_; // weights in each direction
        MdArray<double,full_dynamic_extent_t<LocalDim+1>> control_points_; // control points in each direction
        NurbsBasis<LocalDim> basis_; // basis of the mesh

        // other attributes to be computed
        //Eigen::Matrix<int, Dynamic, 2 * LocalDim, Eigen::RowMajor> neighbors_ {};  // neighbors of each cell
        //Eigen::Matrix<int, Dynamic, LocalDim, Eigen::RowMajor>  cells_ {};  // cells of the mesh
        //Eigen::Matrix<double, Dynamic, LocalDim> param_nodes_ {};  // nodes of the mesh


        //BinaryVector<Dynamic> boundary_markers_ {};   // j-th element is 1 \iff node j is on boundary

        int n_nodes_ = 0, n_cells_ = 0;
        int flags_ = 0;
        std::vector<int> cells_markers_ {};   // marker associated to i-th cell
        std::vector<int> nodes_markers_ {};   // marker associated to i-th node


};



template <int N> class IsoMesh<2, N>: public IsoMeshBase<2, N, IsoMesh<2, N>> {
    fdapde_static_assert(N == 2 || N == 3, THIS_CLASS_IS_FOR_2D_OR_3D_MESHES_ONLY);
    public:
    using Base = IsoMeshBase<2, N, IsoMesh<2, N>>;
    static constexpr int n_nodes_per_edge = 2;
    static constexpr int n_edges_per_cell = 4;
    static constexpr int n_faces_per_edge = 2;
    //using EdgeType = typename Base::CellType::EdgeType; aspettiamo le classi IsoSegment, IsoSquare, IsoCube
    //using LocationPolicy = TreeSearch<IsoMesh<2, N>>;
    //using Base::cells_; //questo manca perché è calcolato a richiesta
    using Base::embed_dim;
    using Base::local_dim;
    using Base::n_cells_;
    using Base::n_nodes_per_cell;
    //static constexpr auto edge_pattern =
    //  Matrix<int, binomial_coefficient(n_nodes_per_cell, n_nodes_per_edge), n_nodes_per_edge>(
    //    combinations(n_nodes_per_edge, n_nodes_per_cell)); // ???

    IsoMesh() = default;

    IsoMesh(std::array<std::vector<double>, 2>& knots, MdArray<double, MdExtents<Dynamic, Dynamic>>& weights,
         MdArray<double, MdExtents<Dynamic,Dynamic,Dynamic>>& control_points, int order, int flags = 0) :
        Base(knots, weights, control_points, order, flags) { }

    const typename Base::CellType& cell(int id) const {
        cell_ = typename Base::CellType(id, this);
        return cell_;
    }

    protected:
    std::vector<int> edges_ {};                        // nodes (as row indexes in nodes_ matrix) composing each edge
    std::vector<int> edge_to_cells_ {};                // for each edge, the ids of adjacent cells
    Eigen::Matrix<int, Dynamic, Dynamic, Eigen::RowMajor> cell_to_edges_ {};   // ids of edges composing each cell
    BinaryVector<Dynamic> boundary_edges_ {};   // j-th element is 1 \iff edge j is on boundary
    std::vector<int> edges_markers_ {};
    int n_edges_ = 0;
    //mutable std::optional<LocationPolicy> location_policy_ {};
    // cell caching
    std::vector<typename Base::CellType> cell_cache_;
    mutable typename Base::CellType cell_;   // used in case cell caching is off

};


}; // namespace fdapde

#endif // __ISO_MESH_H__