#ifndef __ISO_SQUARE_H__
#define __ISO_SQUARE_H__

#include <array>
#include <cmath>
#include <Eigen/Dense>



namespace fdapde {
        
// Forward declaration of IsoMesh
template<int M, int N> class IsoMesh;

template<int M, int N>
class IsoSquare{

    const IsoMesh<M,N>* mesh_; // mesh pointer
    std::size_t ID_; // element ID
    std::array<int, M> left_coords_; // coordinates of the left corner of the element
    std::array<int, M> right_coords_; // coordinates of the right corner of the element
    bool boundary_; // true if the element has at least one vertex on the boundary

    public:

    // Default constructor
    IsoSquare() = default;

    // Constructor with element ID and parent mesh pointer
    IsoSquare(const std::size_t& ID,const IsoMesh<M,N>* mesh): mesh_(mesh), ID_(ID){ 

        left_coords_ = mesh_->compute_lr_vertices(ID_)[0];
        right_coords_ = mesh_->compute_lr_vertices(ID_)[1];

        boundary_ = mesh_->is_boundary(ID_);

    };


    // Affine map from reference domain [-1, 1]^M to parametric domain [left_coords, right_coords]
    std::array<double, M> affine_map(const std::array<double, M> & p) const {
            std::array<double, M> x;
            for(std::size_t i = 0; i < M; ++i){
                x[i] = 0.5*(right_coords_[i] + left_coords_[i] + (right_coords_[i] - left_coords_[i]) * p[i]);
            }
            return x;
        }

    
    // All these methods want a point p in \Omega, [-1,1]^M, and return the corresponding value in the physical domain
    // parametrization gradient F, puo' magari essere resa piu' efficiente

    // also parametrization
    std::array<double,N> parametrization(const std::array<double, M>& p) const {
        // compute x, with the affine map
        auto x = affine_map(p);
        // now x is in [left_coords, right_coords]
        std::array<double,N> F;
        for(int i=0;i<N;i++){
            F[i] = mesh_->eval_param(x,i);
        }
        return F;
    }

    Eigen::Matrix<double,N,M> parametrization_gradient(const std::array<double, M>& p) const {
        // compute x, with the affine map
        auto x = affine_map(p);
        Eigen::Matrix<double,N,M> F;
        for(int i=0;i<N;i++){
            for(int j=0;j<M;j++){
                F(i,j) = mesh_->eval_param_derivative(x,i,j);
            }
        }
        return F;
    }

    // Metric tensor F^T * F
    Eigen::Matrix<double,M,M> metric_tensor(const std::array<double, M>& p) const {
        auto x = affine_map(p);
        auto F = parametrization_gradient(x);
        return F.transpose() * F; 
    }

    // metric determinant sqrt(det(F^T * F))
    double metric_determinant(const std::array<double, M>& p) const {
        auto x = affine_map(p);
        return std::sqrt(metric_tensor(x).determinant()); 
    }

    // Compute the neighbors of the current element
    std::vector<std::size_t> get_neighbors_ID() const {
        return mesh_->get_neighbors_ID(ID_);
    }


    // get the physical coordinates of vertices of a cell (dimension 2^M x N)
    std::vector<std::array<double,N>> compute_physical_vertices() const {
    
    // Number of vertices (2^M)
    std::size_t num_vertices = 1 << M; // 2^M vertices for an M-dimensional cell

    // Store the physical coordinates of each vertex
    std::vector<std::array<double, N>> physical_vertices(num_vertices);

    for (std::size_t vertex_idx = 0; vertex_idx < num_vertices; ++vertex_idx) {
        // Compute the parametric coordinates for this vertex
        std::array<double, M> parametric_coords;
        for (std::size_t dim = 0; dim < M; ++dim) {
            // Use binary representation of `vertex_idx` to determine left or right
            parametric_coords[dim] = (vertex_idx & (1 << dim)) ? right_coords_[dim] : left_coords_[dim];
        }

        // Store the physical coordinates
        physical_vertices[vertex_idx] =  mesh_->compute_physical_coords(parametric_coords);

    }
    return physical_vertices;
    }

    // compute parametric measure

    double parametric_measure() const {
        double measure = 1.0;
        for(std::size_t i = 0; i < M; ++i){
            measure *= right_coords_[i] - left_coords_[i];
        }
        return measure/(1<<M);
    }


    // Getters
    std::size_t ID() const { return ID_; }
    const std::array<int,M> & left_coords() const { return left_coords_; }
    const std::array<int,M> & right_coords() const { return right_coords_; }
    bool is_boundary() const { return boundary_; }

    };

}; // namespace fdapde

#endif // __ISO_SQUARE_H__
