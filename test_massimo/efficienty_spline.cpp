#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <random>
#include <fstream>
#include "utils/symbols.h"
#include "fields/spline.h"
#include "splines_old/basis/old_spline.h"


template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}


int main() {

        // Define a knot vector
        std::vector<double> knots = {0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
        // Instantiate a cubic spline (degree 3) centered at index 2
        const int degree = 4;
        int knotIndex = 4;

        int num_int = 1000;


        auto linspaced = linspace<double>(knots[0], knots[knots.size()-1], num_int );

        constexpr int n=1;

        //fdapde::Spline spline(knots, knotIndex, degree); // Spline nuova
        //auto gradient = spline.gradient(n);
        //fdapde::core::OldSpline<degree> spline(knots,knotIndex); //spline vecchia

        std::vector<int> times;


        for (int i=0; i<linspaced.size(); ++i){
            // per ogni replica valuta la spline 10000 volte in un punto p
            // vantaggio computazionale

            const double p = linspaced[i];
            std::array<double, 1> p_array = {p};
            
            //std::cout<<"Point: "<<p<<std::endl;

            auto t1 = std::chrono::high_resolution_clock::now();

            for(int j = 0; j < 10000; ++j) { 
                //spline.derive<n>()(p); // valuta spline in un punto
                spline(p_array);
            }

            auto t2 = std::chrono::high_resolution_clock::now();
            auto ms_int = duration_cast<std::chrono::milliseconds>(t2 - t1);
            //std::cout<<"Time elapsed: "<<ms_int;
            times.push_back(ms_int.count());
        }

        // Compute max
        double max_time = *std::max_element(times.begin(), times.end());

        // Compute min
        double min_time = *std::min_element(times.begin(), times.end());

        // Compute mean
        double mean_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();

        // Display results
        std::cout << "Max time: " << max_time << std::endl;
        std::cout << "Min time: " << min_time << std::endl;
        std::cout << "Mean time: " << mean_time << std::endl;



        // Save results to a CSV file
        std::ofstream file("results_1der_oldSpline.csv");
        if (file.is_open()) {
            // Write header
            file << "Iteration,Time (ms)" << std::endl;

            // Write times
            for (size_t i = 0; i < times.size(); ++i) {
                file << i + 1 << "," << times[i] << std::endl;
            }

            file.close();
            std::cout << "Results saved to results.csv" << std::endl;
        } else {
            std::cerr << "Unable to open file for writing!" << std::endl;
        }

    return 0;
}