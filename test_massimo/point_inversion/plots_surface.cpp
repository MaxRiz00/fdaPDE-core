#include "isogeometric.h"
#include "helpers.h"

using namespace fdapde;

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> array;
    double step = (b - a) / (n - 1);

    while(a <= b) {
        array.push_back(a);
        a += step;
    }

    return array;
}




int main(){

    using vector_t = Eigen::Matrix<double, Dynamic, 1>;
    using matrix_t = Eigen::SparseMatrix<double>;

    constexpr int M = 3;

    using Vec = Eigen::Matrix<double, M, 1>;
    using Fun = std::function<double(const Vec&)>;
    
    std::string folder = "volume3/";

    std::string path = "../../plots/data/" + folder + "/";

    std::vector<int> ref = {0,1,2,3,4};

    std::ofstream csv_file("../results/timings_nref_twisted.csv");
    csv_file << "n_pts,time_sec\n";  // header

    int N = 46656;

    for (const auto& r : ref){
    
        auto mesh = load_mesh<3,3>(path);
        //auto mesh = IsoMesh<2,3>::torus();
        
        //if(r != 0) 
        mesh.refine_knots({ 1, r,r});

        std::cout<<"n_cells: "<<mesh.n_cells()<<std::endl;

        //int points_per_dir = static_cast<int>(std::round(std::pow(10, r / 2.0)));
        int points_per_dir = std::cbrt(N);

        // linspace of the parametric domain
        std::array<std::vector<double>, 3> test_u;
        test_u[0] = linspace(mesh.knots()[0][0] + 1e-5, mesh.knots()[0][mesh.knots()[0].size()-1]-1e-5, points_per_dir);
        test_u[1] = linspace(mesh.knots()[1][0]+ 1e-5, mesh.knots()[1][mesh.knots()[1].size()-1]-1e-5, points_per_dir);
        test_u[2] = linspace(mesh.knots()[2][0]+ 1e-5, mesh.knots()[2][mesh.knots()[2].size()-1]-1e-5, points_per_dir);

        std::vector<Vec> test_P;
        std::vector<Eigen::Matrix<double, 3, 1>> test_n;
        std::cout<<"# test points: "<<test_u[0].size()*test_u[1].size()<<std::endl;
        

        for(auto u : test_u[0]){
            for(auto v : test_u[1]){
                for(auto w : test_u[2]){
                    test_n.push_back({u, v,w});
                    test_P.push_back(mesh.eval_param({u, v,w}));
                }
            }
        }

        std::cout<<"# test points: "<<test_P.size()<<std::endl;

        // tic
        
        auto start = std::chrono::high_resolution_clock::now();

        Vec P ;
        double t1_sum = 0, t2_sum = 0;
        for(int h=0;h<test_P.size();h++){
            double t1, t2;
            P = test_P[h];
            auto u = test_n[h];
            auto u_ = mesh.invert_point(P,t1,t2);
            //t1_sum += t1;
            //t2_sum += t2;
            //std::cout<<"P: "<<P.transpose()<<" u: "<<u.transpose()<<" u_: "<<u_.transpose()<<std::endl;
        }

        auto end = std::chrono::high_resolution_clock::now();

        // time in seconds double
        auto elapsed_time = std::chrono::duration<double>(end - start).count();
        std::cout << "Total time: " << elapsed_time << " seconds" << std::endl;

        //csv_file << test_P.size() << "," << elapsed_time << "\n";
        csv_file << mesh.n_cells() << "," << elapsed_time << "\n";

        //std::cout<<"t1 mean: "<<t1_sum/test_P.size()<<" t2 mean: "<<t2_sum/test_P.size()<<std::endl;

    }

    csv_file.close();
        
        

    return 0;
}