#ifndef __FDAPDE_ISO_SQUARE_H__
#define __FDAPDE_ISO_SQUARE_H__

#include "header_check.h"

namespace fdapde {

/**
 * @brief 2D parametric square element embedded in physical space.
 * Specialization of IsoCell for squares (rectangles)
 * 
 * @tparam MeshType Parent mesh type (must have local_dim = 1)
 */
template <typename MeshType> class IsoSquare: public IsoCell<MeshType::local_dim, MeshType::embed_dim>{
    fdapde_static_assert(MeshType::local_dim == 2, THIS_CLASS_IS_FOR_2D_MESHES_ONLY);
    using Base = IsoCell<MeshType::local_dim, MeshType::embed_dim>;
    public:
    // === Constructors === //
    IsoSquare() = default;
    IsoSquare(int id, const MeshType* mesh) : 
        id_(id), mesh_(mesh), boundary_(false)  {
        boundary_ = mesh_->is_cell_on_boundary(id_);
        //std::tie(this->left_coords_, this->right_coords_) = mesh_->compute_lr_vertices(id_);
        this->left_coords_ = mesh_->compute_lr_vertices(id_).first;
        this->right_coords_ = mesh_->compute_lr_vertices(id_).second;

        // print left_coords_ << this->left_coords_ << std::endl;
        // print right_coords_ << this->right_coords_ << std::endl;
        // print total number of cells
        //std::cout<<"Total number of cells: " << mesh_->n_cells() << std::endl;
        //std::cout<<"Element ID: " << id_ << std::endl;
        //std::cout << "Left coords: " << this->left_coords_ << std::endl;
        //std::cout << "Right coords: " << this->right_coords_ << std::endl;
    }

    // === Edge Type === //
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

        /**
         * @brief Evaluate the n linspaced physical points for the edge. Only for plotting purposes.
         * 
         * @param n Number of points to evaluate
         * @return Physical coordinates in embedding space
         */
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

    // === Public Member Functions === //

    /**
     * @brief Evaluate the physical point corresponding to a reference coordinate.
     * 
     * @param p Point in reference domain [-1, 1]
     * @return Physical coordinate in embedding space
     */
    Eigen::Matrix<double, MeshType::embed_dim, 1> parametrization(const Eigen::Matrix<double, MeshType::local_dim,1>& p, bool param = false) const {
        if (param)
            return mesh_->eval_param(p);
        else
            return mesh_->eval_param(this->affine_map(p));
    }

    /**
     * @brief Evaluate the Jacobian of the mapping at a reference point.
     * 
     * @param p Point in reference domain [-1, 1]
     * @return First derivative (Jacobian matrix)
     */
    Eigen::Matrix<double, MeshType::embed_dim, MeshType::local_dim> parametrization_gradient(const Eigen::Matrix<double, MeshType::local_dim,1>& p, bool param=false) const {
        //std::cout<<"Param"<<param<<std::endl;
        if(param){
            return mesh_->eval_param_derivatives(p,false).first_derivative;
        }
        else
            return mesh_->eval_param_derivatives(this->affine_map(p),false).first_derivative;
    }


    MdArray<double, MdExtents<MeshType::embed_dim, MeshType::local_dim, MeshType::local_dim>> parametrization_hessian(const Eigen::Matrix<double, MeshType::local_dim,1>& p, bool param=false) const {
        //std::cout<<"Param"<<param<<std::endl;
        if(param){
            return *(mesh_->eval_param_derivatives(p,true).second_derivative);
        }
        else
            return *(mesh_->eval_param_derivatives(this->affine_map(p),true).second_derivative);
    }

    /**
     * @brief Compute the metric tensor at a reference point. Computes Fᵀ·F where F is the Jacobian of the mapping.
     * 
     * @param p Point in reference domain [-1, 1]
     * @return Symmetric metric tensor matrix
     */
    Eigen::Matrix<double, MeshType::local_dim, MeshType::local_dim> metric_tensor(const Eigen::Matrix<double, MeshType::local_dim,1>& p, bool param=false) const {
        auto F = parametrization_gradient(p,param);
        return F.transpose() * F; 
    }

    /**
     * @brief Compute the square root of the determinant of the metric tensor.
     * 
     * @param p Point in reference domain [-1, 1]
     * @return Determinant of the metric tensor (metric scaling factor)
     */
    double metric_determinant(const Eigen::Matrix<double, MeshType::local_dim,1>& p, bool param=false) const {
        return std::sqrt(metric_tensor(p,param).determinant()); 
    }

    /**
     * @brief Evaluate a regular grid of physical points over the square element. Only for plotting purposes.
     * 
     * Performs a 2D tensor-product linear interpolation in parametric space,
     * maps each point to physical space, and stores the result.
     * 
     * @param n Number of evaluation points per parametric direction (produces n x n grid)
     * @return 3D array (n x n x embed_dim) of physical coordinates
     */
    MdArray<double, full_dynamic_extent_t<MeshType::local_dim + 1>> linspace_evaluation(int n, MdArray<double, full_dynamic_extent_t<MeshType::local_dim + 1>>& parametric_points ) const {
        MdArray<double, full_dynamic_extent_t<MeshType::local_dim + 1>> res(n, n, MeshType::embed_dim);
        auto param_nodes = mesh_->parametric_nodes();
        auto left_coords = this->left_coords_;
        auto right_coords = this->right_coords_;

        // store the parametric points
        parametric_points.resize(n, n, MeshType::local_dim);

        // compute the step
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Eigen::Matrix<double, MeshType::local_dim, 1> p;
                auto t1 = static_cast<double>(i) / (n - 1);
                auto t2 = static_cast<double>(j) / (n - 1);

                p(0) = (1 - t1) * left_coords(0) + t1 * right_coords(0);  // Linear interpolation
                p(1) = (1 - t2) * left_coords(1) + t2 * right_coords(1);  // Linear interpolation

                for(int k = 0; k < MeshType::local_dim; k++){
                    parametric_points(i, j, k) = p(k);
                }

                auto param = mesh_->eval_param(p);
                for(int k = 0; k < MeshType::embed_dim; k++){
                    res(i, j, k) = param(k);
                }
            }
        }
        return res;

    }


    // === Getters === // 
    int id() const { return id_; }
    Eigen::Matrix<int, 1, 2 * MeshType::local_dim> neighbors() const { return mesh_->neighbors().row(id_); }
    Eigen::Matrix<int, 1, MeshType::local_dim> node_ids() const { return mesh_->cells().row(id_); }
    bool on_boundary() const { return boundary_; }
    operator bool() const { return mesh_ != nullptr; }
    EdgeType edge(int n){
        fdapde_assert(n<this->n_edges);
        return EdgeType(mesh_->cell_to_edes()(id_,n), mesh_);
    }
    int marker() const {return mesh_->cell_markers().size() ? mesh_->cells_markers()[id_] : Unmarked; }

    // === Edge Iterators === //
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
    int id_ = 0;                     ///< id of the square element
    const MeshType* mesh_ = nullptr; ///< pointer to the parent mesh
    bool boundary_ = false;          ///< true if square element is on the boundary
};
    
    
}; // namespace fdapde

#endif // __FDAPDE_ISO_SQUARE_H__