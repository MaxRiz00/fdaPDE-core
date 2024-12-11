#include <iostream>
#include <vector>
#include "utils/symbols.h"
#include "fields/spline.h"
# include "splines/spline_basis.h"


int main() {

        // Define a knot vector
        DVector<double> knots(11); // Specify the size
        knots << 0,0,0,1,2,3,4,4,5,5,5;

        // Instantiate a cubic spline (degree 3) centered at index 2
        int degree = 2;
        int knotIndex = 2;
        fdapde::Spline spline(knots, knotIndex, degree);

        // Evaluate the spline at a given point
        double point = 1.5; // Parameter value to evaluate the spline
        double splineValue = spline(point);
        std::cout << "Spline value at " << point << ": " << splineValue << std::endl;

        int n=2;

        // // Compute the gradient (derivative)
           auto gradient = spline.gradient(n);
           double derivativeValue = gradient(point);
           std::cout << "Spline "<<n<<"-th derivative at " << point << ": " << derivativeValue << std::endl;


        fdapde::SplineBasis spline_base(knots, degree);

        for (int i =0 ;i<spline_base.size();i++){
            std::cout<< spline_base[i](point)<<std::endl;
        }

    return 0;
}