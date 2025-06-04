#include "isogeometric.h"
#include "exact_solution.h"
#include "../helpers.h"
#include "utils/utils.h"

using namespace fdapde;
int main() {
    typedef Eigen::Matrix<double, 2, 1> Point2D;

    // Define degree-2 full NURBS circle
    std::vector<double> knots = {
        0.0, 0.0, 0.0,
        0.25, 0.25,
        0.5, 0.5,
        0.75, 0.75,
        1.0, 1.0, 1.0
    };

    MdArray<double, MdExtents<Dynamic>> weights(9);
    MdArray<double, MdExtents<Dynamic, Dynamic>> control_points(9, 2);

    double w = std::sqrt(2) / 2;

    weights(0) = 1.0;
    weights(1) = w;
    weights(2) = 1.0;
    weights(3) = w;
    weights(4) = 1.0;
    weights(5) = w;
    weights(6) = 1.0;
    weights(7) = w;
    weights(8) = 1.0;

    control_points(0, 0) = 1.0;  control_points(0, 1) = 0.0;
    control_points(1, 0) = 1.0;  control_points(1, 1) = 1.0;
    control_points(2, 0) = 0.0;  control_points(2, 1) = 1.0;
    control_points(3, 0) = -1.0; control_points(3, 1) = 1.0;
    control_points(4, 0) = -1.0; control_points(4, 1) = 0.0;
    control_points(5, 0) = -1.0; control_points(5, 1) = -1.0;
    control_points(6, 0) = 0.0;  control_points(6, 1) = -1.0;
    control_points(7, 0) = 1.0;  control_points(7, 1) = -1.0;
    control_points(8, 0) = 1.0;  control_points(8, 1) = 0.0;

    IsoMeshData<1> curve(knots, weights, control_points, 2);

    // Degree elevation (increase to degree 3)
    IsoMeshData<1> elevated = iso_algorithms::degree_elevation(curve, 3);
    std::cout << "Original degree: " << curve.degree[0] << ", Elevated degree: " << elevated.degree[0] << std::endl;

    // Evaluation parameter values
    std::vector<double> test_params = {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};

    // Print degree-elevated structure
    std::cout << "\n[Degree Elevation] Knots: ";
    for (const auto& knot : elevated.knots[0]) std::cout << knot << " ";
    std::cout << "\nWeights: ";
    for (int i = 0; i < elevated.weights.extent(0); ++i) std::cout << elevated.weights(i) << " ";
    std::cout << "\nControl points:\n";
    for (int i = 0; i < elevated.control_points.extent(0); ++i) {
        for (int j = 0; j < elevated.control_points.extent(1); ++j) {
            std::cout << elevated.control_points(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Build meshes
    IsoMesh<1,2> new_mesh(elevated.knots, elevated.weights, elevated.control_points, elevated.degree);
    IsoMesh<1,2> old_mesh(curve.knots, curve.weights, curve.control_points, curve.degree);

    // Elevate IsoMesh object (p-refinement)
    IsoMesh<1,2> p_refined_mesh(curve.knots, curve.weights, curve.control_points, curve.degree);
    p_refined_mesh.elevate_degree({1});
    std::cout << "\n[P-Refinement] Degree: " << p_refined_mesh.degree()[0] << std::endl;

    std::cout << "\n[P-Refinement] Knots: ";
    for (const auto& knot : p_refined_mesh.knots()[0]) std::cout << knot << " ";
    std::cout << "\nWeights: ";
    for (int i = 0; i < p_refined_mesh.weights().extent(0); ++i) std::cout << p_refined_mesh.weights()(i) << " ";
    std::cout << "\nControl points:\n";
    for (int i = 0; i < p_refined_mesh.control_points().extent(0); ++i) {
        for (int j = 0; j < p_refined_mesh.control_points().extent(1); ++j) {
            std::cout << p_refined_mesh.control_points()(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Compare geometry at test points
    Eigen::Matrix<double,1,1> u_param;
    std::cout << "\n[Comparison] Degree-elevated vs p-refined at sample points:\n";
    for (double u : test_params) {
        u_param(0) = u;
        Point2D pt_elev = new_mesh.eval_param(u_param);
        Point2D pt_pref = p_refined_mesh.eval_param(u_param);

        for (int d = 0; d < 2; ++d) {
            double diff = std::abs(pt_elev(d) - pt_pref(d));
            if (diff > 1e-8) {
                std::cout << "u = " << u << ", dim = " << d
                          << ", elev = " << pt_elev(d)
                          << ", pref = " << pt_pref(d)
                          << ", |diff| = " << diff << std::endl;
            }
        }
    }

    return 0;
}