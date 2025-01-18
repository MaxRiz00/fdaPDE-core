#ifndef __ISO_MESH_H__
#define __ISO_MESH_H__

#include "nurbs_basis.h"
#include "iso_square.h"
#include "linear_algebra/mdarray.h"

namespace fdapde {

template<int M, int N>
class IsoMesh {

    //template<typename T> using DMatrix = MdArray<T, Dynamic, Dynamic>;

    std::array<std::vector<double>,M> knots_; // knots in each direction
    MdArray<double,full_dynamic_extent_t<M>> weights_; // weights in each direction
    MdArray<double,full_dynamic_extent_t<M+1>> control_points_; // control points in each direction

    int n_elements_; // number of elements
    int n_nodes_; // number of nodes

    
    NurbsBasis<M> basis_; // basis of the mesh


    //std::vector<ElementIGA<M,N>> elements_cache_; // elements of the mesh

    //DMatrix<double> nodes_ {}; // nodes of the mesh, size: r by M, where r is the product of the numbers of unique knots along each dimension
    //std::vector<std::size_t> boundary_ {}; // : r by 1 contains a boolean value for each node, true if the node is on the boundary of Ω
    //DMatrix<std::size_t> elements_ {}; // elements of the mesh
    //DMatrix<std::size_t> boundary_dof_ {}; // dofs of boundary of the mesh

    // implementation of the private methods for the mesh parametrization

    std::size_t compute_stride(std::size_t dim, bool element) const {
        return element ? this->knots_[dim].size() - 1 : this->knots_[dim].size();   
    }   


    inline double eval_param(const std::array<double, M>& u, int sliceIdx) const {
        auto CpSlice = this->control_points_.template slice<M>(sliceIdx);
        double x = 0.0;
        for (const auto& nurb : this->nurbs_basis_) {
            x += nurb(x) * CpSlice(nurb.index());
        }
        return x;
    }

    inline double eval_param_derivative(const std::array<double, M>& u, int sliceIdx, std::size_t derivative_index) const {
        auto CpSlice = this->control_points_.template slice<M>(sliceIdx);
        double x = 0.0;
        for (const auto& nurb : this->nurbs_basis_) {
            x += nurb.derive(derivative_index)(x) * CpSlice(nurb.index());
        }
        return x;
    }


    // HERE

    public:

    IsoMesh(std::array<std::vector<double>,M> & knots, MdArray<double,full_dynamic_extent_t<M>> & weights, int order,
        MdArray<double,full_dynamic_extent_t<M+1>> & control_points): knots_(knots), weights_(weights), control_points_(control_points){

            basis_ = NurbsBasis<M>(knots, weights, order);

            // compute numeber of elements and nodes
            n_elements_ = 1;
            n_nodes_ = 1;

            for(std::size_t i=0; i < M; ++i){
                n_elements_ *= knots_[i].size() - 1;
                n_nodes_ *= knots_[i].size();
            }

            //std::cout<<"n_elements: "<<n_elements_<<std::endl;

        };

    
    // function to compute the ID of a cell from the multi-index

    // function to compute the ID of a cell from the multi-index

    std::size_t compute_el_ID(const std::array<std::size_t, M>& eMultiIndex) const {
        std::size_t ID = eMultiIndex[M - 1]; // Start with the last index
        for (std::size_t i = M - 1; i > 0; --i) { // Iterate backward
            ID = ID * (this->knots_[i - 1].size() - 1) + eMultiIndex[i - 1];
        }
        return ID;
    }

    std::array<std::size_t, M> compute_multiIndex(std::size_t ID, bool element=true) const {
        std::array<std::size_t, M> eMultiIndex;
        for (std::size_t i = M; i > 0; --i) {
            if (i == 1) {
                eMultiIndex[M-1] = ID; // The last index is the remaining ID
            } else {
                std::size_t stride = element ? this->knots_[i - 2].size() - 1 : this->knots_[i - 2].size(); // Use (i-2) to match dimensions
                eMultiIndex[M - i] = ID % stride; // Extract the current index
                ID /= stride; // Update ID for the next dimension
            }
        }
        return eMultiIndex;
    }

    std::array<std::array<int, M>,2> compute_param_el_vertices(const std::size_t ID)const {
        auto eMultiIndex = compute_multiIndex(ID);
        std::array<std::array<int, M>,2> vertices;
        for(std::size_t i = 0; i < M; ++i){
            vertices[0][i] = eMultiIndex[i];
            vertices[1][i] = eMultiIndex[i] + 1;
        }
        return vertices;
    }

    bool is_boundary(const std::size_t ID) const {
        auto eMultiIndex = compute_multiIndex(ID);
        for(std::size_t i = 0; i < M; ++i){
            if(eMultiIndex[i] == 0 || eMultiIndex[i] == knots_[i].size() - 2){
                return true;
            }
        }
        return false;
    }

    bool is_node_boundary(const std::size_t ID) const {
        auto eMultiIndex = compute_multiIndex(ID, false);
        for(std::size_t i = 0; i < M; ++i){
            if(eMultiIndex[i] == 0 || eMultiIndex[i] == knots_[i].size() - 1){
                return true;
            }
        }
        return false;
    }

    // compute nodes, tensor for each node, compute the multi-index of the node, dimension r (number of nodes) by M

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

    // compute elements

    std::vector<std::array<std::size_t,M>> compute_elements() const {
        std::vector<std::array<std::size_t,M>> elements;
        elements.resize(n_elements_);
        for(std::size_t i = 0; i < n_elements_; ++i){
            elements[i] = compute_multiIndex(i);
        }
        return elements;
    }

    // get neighbors

    std::vector<std::size_t> get_neighbors_ID(const std::size_t ID) const {
        // compute the multi-index of the element
        auto eMultiIndex = compute_multiIndex(ID);
        std::vector<std::size_t> neighbors;

        // Size of the neighbors vector is at most 2M
        // each element's neighbors are elements who have the same indexes along each dimension
        // except for one dimension where the index is the previous or the next one (if they exist)
        // w.r.t. the considered element
        // to simplify the computation, we add the "next" element along a direction to the current one's neighbor list
        // while simultaneously adding the current one to the "next" one's list

        // loop over all the dimensions
        for(std::size_t i = 0; i < M; ++i){
            // create a copy of the multi-index
            auto neighbor = eMultiIndex;
            // check if the element is not on the boundary
            if(eMultiIndex[i] > 0){
                // add the previous element along the i-th direction
                neighbor[i] = eMultiIndex[i] - 1;
                neighbors.push_back(compute_el_ID(neighbor));
            }
            // create a copy of the multi-index
            neighbor = eMultiIndex;
            // check if the element is not on the boundary
            if(eMultiIndex[i] < knots_[i].size() - 2){
                // add the next element along the i-th direction
                neighbor[i] = eMultiIndex[i] + 1;
                neighbors.push_back(compute_el_ID(neighbor));
            }
        }
        // sort the neighbors
        std::sort(neighbors.begin(), neighbors.end());
        return neighbors;
    }

    // Da mettere coordinate fisiche, problema con MdArray che son supporta somme di array
    // This function computes the phisical coordinates of a point given its parametric coordinates
    /*
    std::array<double, M+1> compute_physical_coords(const std::array<double, M> & p) const {
        std::array<double, M+1> x;

        for (const auto& nurb : this->nurbs_basis_) {
            x += nurb(x) * control_points_(nurb.index()); // devi fare una slice e poi sommare ... Da rivedere 
        }
        return x;
    }
    */

    // getCell function
    
    IsoSquare<M,N> getCell(std::size_t index) const {
        // compute the multi-index of the element
        auto eMultiIndex = compute_multiIndex(index);
        return IsoSquare<M,N>(eMultiIndex, this);
    }

    // some getters

    const NurbsBasis<M>& basis() const { return basis_; }
    const MdArray<double, full_dynamic_extent_t<M+1>>& control_points() const { return control_points_; }
    const std::array<std::vector<double>,M>& knots() const { return knots_; }

    public:

    class CellIterator {
    public:

        CellIterator(const IsoMesh* mesh, size_t index): parentMesh_(mesh), currentIndex_(index), currentCell_(mesh->getCell(index)) {}

        // Dereferenziazione
        IsoSquare<M,N> operator*() const {
            return currentCell_;
        }

        // Incremento
        CellIterator& operator++() {
            currentIndex_++;
            currentCell_ = parentMesh_->getCell(currentIndex_); // Update cached cell
            return *this;
        }

        // Confronto
        bool operator==(const CellIterator& other) const {
            return currentIndex_ == other.currentIndex_ ;
        }

        bool operator!=(const CellIterator& other) const {
            return !(currentIndex_ == other.currentIndex_ );
        }

    private:
        const IsoMesh* parentMesh_;
        size_t currentIndex_;
        IsoSquare<M,N> currentCell_; // Cache the current cell
    };

    // Metodi per ottenere gli iteratori
    CellIterator beginCells() const {
        return CellIterator(this, 0);
    }

    CellIterator endCells() const {
        return CellIterator(this, this->n_elements_ - 1);
    }


    // some getters
    std::size_t n_elements() const { return n_elements_; }
    std::size_t n_nodes() const { return n_nodes_; }

};



}; // namespace fdapde

#endif // __ISO_MESH_H__