#ifndef __ISO_ALGORITHMS_H__
#define __ISO_ALGORITHMS_H__

#include "header_check.h"

namespace fdapde {

// Algorithm to refine a knot vector
template <int EmbedDim>
void knot_refinement(int degree,
                    const std::vector<double>& old_knots,
                    const std::vector<double>& insert_knots,
                    std::vector<double>& new_knots,
                    const std::vector<std::array<double,EmbedDim+1>>& old_cpw,
                    std::vector<std::array<double,EmbedDim+1>>& new_cpw) {
    
    // get the number of control points
    int n = old_cpw.size() - 1;
    int cp_size = EmbedDim+1;
    int m = degree + n + 1;

    std::cout<<"n = "<<n<<std::endl;
    std::cout<<"m = "<<m<<std::endl;

    int r = insert_knots.size() - 1;

    std::cout<<"Old knots size: "<<old_knots.size()<<std::endl;
    std::cout<<"Insert knots size: "<<insert_knots.size()<<std::endl;

    // get the span
    auto old_basis  = BSplineBasis(old_knots, degree);
    int a = old_basis.find_span(insert_knots[0]);
    int b = old_basis.find_span(insert_knots[r]) + 1;

    std::cout<<"b = "<<b<<std::endl;

    // get the new control points
    for(int j=0; j<=a-degree; j++) {
        for(int i=0; i<cp_size; i++) {
            new_cpw[j][i] = old_cpw[j][i];
        }
    }

    for(int j=b-1; j<=n; j++) {
        for(int i=0; i<cp_size; i++) {
            new_cpw[j+r+1][i] = old_cpw[j][i];
        }
    }
    // get the new knots
    
    for(int j=0; j<=a; j++) {
        new_knots[j] = old_knots[j];}
    for(int j=b+degree; j<=m+1; j++) new_knots[j+r+1] = old_knots[j]; // ?????????

    // get the new control points
    int i = b + degree - 1;
    int k = b + degree + r;

    std::cout<<"k = "<<k<<std::endl;


    for(int j=r; j>=0; j--) {
        while(insert_knots[j] <= old_knots[i] && i > a) {
            for(int l=0; l<cp_size; l++) {
                new_cpw[k-degree-1][l] = old_cpw[i-degree-1][l];
            }
            new_knots[k] = old_knots[i];
            k = k - 1;
            i = i - 1;
        }
        for(int l =0;l<cp_size;l++) {
            new_cpw[k-degree-1][l] = new_cpw[k-degree][l];
        }
        for(int l = 1; l<=degree; l++) {
            int ind = k-degree+l;
            double alpha = new_knots[k+l] - insert_knots[j];
            if(alpha == 0.0) {
                for(int m=0; m<cp_size; m++) {
                    new_cpw[ind-1][m] = new_cpw[ind][m];
                }
            } else {
                alpha = alpha / (new_knots[k+l] - old_knots[i-degree+l]);
                for(int m=0; m<cp_size; m++) {
                    new_cpw[ind-1][m] = alpha * new_cpw[ind-1][m] + (1.0 - alpha) * new_cpw[ind][m];
                }
            }
        }
        new_knots[k] = insert_knots[j];
        k = k - 1;  
    }

}

} // namespace fdapde

#endif //__ISO_ALGORITHMS_H__