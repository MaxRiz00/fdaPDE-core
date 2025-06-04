#include "isogeometric.h"
#include "exact_solution.h"
#include "../helpers.h"
#include "utils/utils.h"

using namespace fdapde;
int main() {
    typedef Eigen::Matrix<double, 2, 1> Point2D;

    auto sphere = IsoMesh<2, 3>::sphere();

    auto sphere_ec = IsoMesh<2, 3>::sphere();
    sphere_ec.elevate_degree({2, 2});

    // Define test parameters  as vector of some 2d points between 0 and 1
    std::vector<Point2D> test_params = {
        {0.0, 0.0}, {0.125, 0.125}, {0.25, 0.25}, {0.375, 0.375},
        {0.5, 0.5}, {0.625, 0.625}, {0.75, 0.75}, {0.875, 0.875}, {1.0, 1.0}
    };


    // Compare geometry at test points
    std::cout << "\n[Comparison] Degree-elevated vs p-refined at sample points:\n";
    for (auto u : test_params) {
        auto  pt_elev = sphere.eval_param(u);
        auto  pt_pref = sphere_ec.eval_param(u);

        for (int d = 0; d < 3; ++d) {
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