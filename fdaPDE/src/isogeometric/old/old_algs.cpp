    // evaluate the derivative of the physical coordinates of a point given its parametric coordinates with index sliceIdx and derivative_index as the derivative index
    // puoi ritornare Eigen::Matrix<double, N, M>, direttamente lo jacobiano
    Eigen::Matrix<double, N, M, Eigen::RowMajor> eval_param_derivative_(const std::array<double, M>& u) const {
        Eigen::Matrix<double, N, M, Eigen::RowMajor> x = Eigen::Matrix<double, N, M>::Zero();
        for(int i = 0;i<N;i++){
            for(int j = 0; j<M; j++){
                const auto cp_slice = this->control_points_.template slice<M>(i);
                for (const auto& nurb : this->basis_) {
                    x(i,j) += nurb.derive(j)(u) * cp_slice(nurb.index());
                }
            }
        }
        return x;
    }


    // mettili privati e basta, puoi outputtare  Eigen::Matrix<double, N, 1> per fare i calcoli
    // Evaluate the physical coordinates of a point u given its parametric coordinates with index sliceIdx
    // falla u che diventa in [-1,1]^M, portando la affine_map in questa classe
    Eigen::Matrix<double, N, 1> eval_param_(const std::array<double, M>& u) const {
        Eigen::Matrix<double, N, 1> x = Eigen::Matrix<double, N, 1>::Zero();
        x.setZero();
        for(int i = 0;i<N;i++){
            const auto cp_slice = this->control_points_.template slice<M>(i); // cp_slice , slice_idx
            for (const auto& nurb : this->basis_) {
                x(i) += nurb(u) * cp_slice(nurb.index());  
            }
        }
        return x;
    }