#ifndef __ISO_MESH_H__
#define __ISO_MESH_H__


#include "nurbs.h"
#include "iso_square.h"

namespace fdapde {

template<int M, int N>
class IsoMesh {

    std::array<std::vector<double>,M> knots_; // knots in each direction
    MdArray<double,full_dynamic_extent_t<M>> weights_; // weights in each direction
    MdArray<double,full_dynamic_extent_t<M+1>> control_points_; // control points in each direction

    int n_elements_; // number of elements
    int n_nodes_; // number of nodes

    NurbsBasis<M> basis_; // basis of the mesh

    //DMatrix<std::size_t> boundary_dof_ {}; // dofs of boundary of the mesh ( DofHandler.h will be implemented later)

    // Compute stride based on dimension and type (element or node)
    std::size_t compute_stride(std::size_t dim, bool element) const {
        return element ? this->knots_[dim].size() - 1 : this->knots_[dim].size();   
    }   

    public:

    // Evaluate the physical coordinates of a point u given its parametric coordinates with index sliceIdx
    inline double eval_param(const std::array<double, M>& u, int sliceIdx) const {
        const auto CpSlice = this->control_points_.template slice<M>(sliceIdx);
        double x = 0.0;
        for (const auto& nurb : this->basis_) {
            x += nurb(u) * CpSlice(nurb.index());
        }
        return x;
    }

    // evaluate the derivative of the physical coordinates of a point given its parametric coordinates with index sliceIdx and derivative_index as the derivative index
    inline double eval_param_derivative(const std::array<double, M>& u, int sliceIdx, std::size_t derivative_index) const {
        const auto CpSlice = this->control_points_.template slice<M>(sliceIdx);
        double x = 0.0;
        for (const auto& nurb : this->basis_) {
            x += nurb.derive(derivative_index)(u) * CpSlice(nurb.index());
        }
        return x;
    }


    // Compute the ID of a cell (or a node if element is false) from the multi-index
    std::size_t compute_ID(const std::array<std::size_t, M>& eMultiIndex, bool element=true) const {
        std::size_t ID = eMultiIndex[M - 1]; // Start with the last index
        for (std::size_t i = M - 1; i > 0; --i) { // Iterate backward
            ID = ID * compute_stride(i-1,element) + eMultiIndex[i - 1];
        }
        return ID;
    }

    public:

    // Constructor
    IsoMesh(std::array<std::vector<double>,M> & knots, MdArray<double,full_dynamic_extent_t<M>> & weights, int order,
         MdArray<double,full_dynamic_extent_t<M+1>> & control_points): knots_(knots), weights_(weights), control_points_(control_points){

            basis_ = NurbsBasis<M>(knots, weights, order);

            // compute number of elements and nodes
            n_elements_ = 1;
            n_nodes_ = 1;
            for(std::size_t i=0; i < M; ++i){
                n_elements_ *= knots_[i].size() - 1;
                n_nodes_ *= knots_[i].size();
            }
        };



    // Compute the multi-index of a cell (or a node if element is false) from the ID
    std::array<std::size_t, M> compute_multiIndex(const std::size_t& ID, bool element=true) const {
        auto tempID = ID;
        std::array<std::size_t, M> eMultiIndex;
        for (std::size_t i = M; i > 0; --i) {
            if (i == 1) {
                eMultiIndex[M-1] = tempID; // The last index is the remaining ID
            } else {
                std::size_t stride = compute_stride(i - 2, element); // Use (i-2) to match dimensions
                eMultiIndex[M - i] = tempID % stride; // Extract the current index
                tempID /= stride; // Update ID for the next dimension
            }
        }
        return eMultiIndex;
    }

    // Compute the vertices of a square (cell) given its ID
    std::array<std::array<int, M>,2> compute_lr_vertices(const std::size_t& ID)const {
        auto eMultiIndex = compute_multiIndex(ID);
        std::array<std::array<int, M>,2> vertices;
        for(std::size_t i = 0; i < M; ++i){
            vertices[0][i] = eMultiIndex[i];
            vertices[1][i] = eMultiIndex[i] + 1;
        }
        return vertices;
    }

    // check if a cell/node is on the boundary
    bool is_boundary(const std::size_t& ID, bool element=true) const {
        auto eMultiIndex = compute_multiIndex(ID,element);
        for(std::size_t i = 0; i < M; ++i){
            if(eMultiIndex[i] == 0 || eMultiIndex[i] == compute_stride(i,element) - 1){
                return true;
            }
        }
        return false;
    }
    
    // Compute all nodes of the mesh (dimension #nodes x M )
    std::vector<std::array<double,M>> compute_nodes() const {
        std::vector<std::array<double,M>> nodes;
        // print number of nodes
        nodes.resize(n_nodes_);
        for(std::size_t i = 0; i < n_nodes_; ++i){
            auto eMultiIndex = compute_multiIndex(i, false);
            for(std::size_t j = 0; j < M; ++j){
                nodes[i][j] = knots_[j][eMultiIndex[j]];
            }
        }
        return nodes;
    }

    // Compute all elements of the mesh (dimension #elements x M )
    std::vector<std::array<std::size_t,M>> compute_elements() const {
        std::vector<std::array<std::size_t,M>> elements;
        elements.resize(n_elements_);
        for(std::size_t i = 0; i < n_elements_; ++i){
            elements[i] = compute_multiIndex(i);
        }
        return elements;
    }

    // Compute the ID neighbors of a cell given its ID (dimension #elements x 2M(at most))
    std::vector<std::size_t> get_neighbors_ID(const std::size_t& ID) const {
        auto eMultiIndex = compute_multiIndex(ID);
        std::vector<std::size_t> neighbors;
        bool element = true; // neighbors are elements

        for(std::size_t i = 0; i < M; ++i){
            auto neighbor = eMultiIndex;
            if(eMultiIndex[i] > 0){
                neighbor[i] = eMultiIndex[i] - 1;
                neighbors.push_back(compute_ID(neighbor));
            }
            neighbor = eMultiIndex;
            if(eMultiIndex[i] < compute_stride(i,element) - 1 ){
                neighbor[i] = eMultiIndex[i] + 1;
                neighbors.push_back(compute_ID(neighbor));
            }
        }
        std::sort(neighbors.begin(), neighbors.end());
        return neighbors;
    }

    // Compute the coordinates of the physical coords using eval_param method of the mesh
    // point in Parametric Space (M dim) -> point in Physical Space (M+1 dim)
    std::array<double,N> compute_physical_coords(const std::array<double, M> & p) const {
        std::array<double, N> x;
        // note that the control points are stored in a tensor of dimension M+1, we need to slice it
        for(int i = 0; i < N; ++i){
            x[i] = eval_param(p,i);
        }
        return x;
    }

    // Access individual cells, Create a cell object given its ID
    IsoSquare<M,N> getCell(const std::size_t& ID) const {
        return IsoSquare<M,N>(ID, this);
    }


    class CellIterator {
    private:
        const IsoMesh* parentMesh_;
        std::size_t currentIndex_;
        IsoSquare<M,N> currentCell_; 
    
    public:

        CellIterator(const IsoMesh* mesh, size_t index): 
            parentMesh_(mesh), currentIndex_(index), currentCell_(mesh->getCell(index)) {}

        IsoSquare<M,N> operator*() const { return currentCell_;}

        CellIterator& operator++() {
            currentIndex_++;
            currentCell_ = parentMesh_->getCell(currentIndex_); 
            return *this;
        }

        // Confronto
        bool operator==(const CellIterator& other) const { return currentIndex_ == other.currentIndex_ ;}

        bool operator!=(const CellIterator& other) const { return !(*this == other);}
    };

    // Metodi per ottenere gli iteratori
    CellIterator beginCells() const { return CellIterator(this, 0);}
    CellIterator endCells() const { return CellIterator(this, this->n_elements_);}

    class BoundaryIterator {
        private:
            const IsoMesh* parentMesh_;
            std::size_t currentIndex_;
            IsoSquare<M,N> currentCell_;
        public:
            BoundaryIterator(const IsoMesh* mesh, std::size_t index)
                : parentMesh_(mesh), currentIndex_(index), currentCell_(mesh->getCell(index)) {
                while (currentIndex_ < parentMesh_->n_elements_ && !parentMesh_->is_boundary(currentIndex_)) {
                    currentIndex_++;
                }
            }
            IsoSquare<M,N> operator*() const { return currentCell_; }
            BoundaryIterator& operator++() {
                do {
                    currentIndex_++;
                } while (currentIndex_ < parentMesh_->n_elements_ && !parentMesh_->is_boundary(currentIndex_));
                if (currentIndex_ < parentMesh_->n_elements_) {
                    currentCell_ = parentMesh_->getCell(currentIndex_);
                }
                return *this;
            }
            bool operator==(const BoundaryIterator& other) const { return currentIndex_ == other.currentIndex_; }
            bool operator!=(const BoundaryIterator& other) const { return !(*this == other); }
        };

    BoundaryIterator beginBoundaryCells() const { return BoundaryIterator(this, 0);}
    BoundaryIterator endBoundaryCells() const { return BoundaryIterator(this, n_elements_); }


    // some getters
    const NurbsBasis<M>& basis() const { return basis_; }
    const MdArray<double, full_dynamic_extent_t<M+1>>& control_points() const { return control_points_; }
    const std::array<std::vector<double>,M>& knots() const { return knots_; }

    std::size_t n_elements() const { return n_elements_; }
    std::size_t n_nodes() const { return n_nodes_; }
};
}; // namespace fdapde

#endif // __ISO_MESH_H__