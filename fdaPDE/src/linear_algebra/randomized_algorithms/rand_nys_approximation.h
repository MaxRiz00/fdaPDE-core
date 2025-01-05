//
// Created by Marco Galliani on 02/10/24.
//

#ifndef NYSTROM_APPROXIMATION_H
#define NYSTROM_APPROXIMATION_H

#include <algorithm>
#include <random>
#include <unordered_set>
#include <Eigen/Cholesky>
#include <Eigen/QR>

#include "utils/random_utils.h"
#include "fdaPDE/utils/symbols.h"

namespace fdapde{
namespace core{

template<typename MatrixType, IterationPolicy ItPolicy> class NysIterations;

//Nystrom Low-rank Approximator
template<typename MatrixType,IterationPolicy ItPolicy>
class NystromApproximator{
private:
    //method parameters
    int block_sz_;
    double shift_;
    double tolerance_=1e-3;
    int seed_=fdapde::random_seed;
    //Storage of the decomposition
    DMatrix<double> Z_; //factor matrix: A + shift*I = Z*Z^T
    //Iterator
    friend NysIterations<MatrixType,ItPolicy>;
    NysIterations<MatrixType,ItPolicy> InitNysIterations(const MatrixType &A){
        return NysIterations<MatrixType,ItPolicy>(A, this);
    };
public:
    //constructors
    NystromApproximator() = default;
    NystromApproximator(int block_sz, int seed=fdapde::random_seed) : block_sz_(block_sz), seed_(seed){}
    //computation
    NystromApproximator compute(const MatrixType &A){
        for(auto rf_it=InitNysIterations(A); !rf_it.stop(); ++rf_it){}
        return *this;
    }
    //setters
    inline void setBlockSize(int block_sz){ block_sz_ = block_sz;}
    inline void setSeed(int seed){ seed_ = seed;}
    inline void setTol(double tol){ tolerance_ = tol;}

    //getters
    inline DMatrix<double>& matrixF(){ return Z_;}
    inline double shift(){ return shift_;}
    //freeing the space occupied by the storage of the decompositions
    void flush(){
        Z_.resize(0,0);
    }
};

//NysBKI: Nystrom Block Krylov Iterations
//assumption -> A is sdp
template<typename MatrixType>
class NysIterations<MatrixType,IterationPolicy::BlockKrylovIterations>{
private:
    int index_,maxIter_;
    DMatrix<double> A_; //matrix to be decomposed
    double norm_A_;
    //Storage of the decomposition
    DMatrix<double> X_; //range and corange
    DMatrix<double> T_; //core matrix
    //Method parameters
    size_t block_sz_;
    double shift_;
    //Pointer to the calling Range Finder
    NystromApproximator<MatrixType,IterationPolicy::BlockKrylovIterations> *nys_approx_;
public:
    NysIterations(const MatrixType &A, NystromApproximator<MatrixType,IterationPolicy::BlockKrylovIterations>* nys_approx) :
    A_(A),nys_approx_(nys_approx){
        index_=0;
        block_sz_ = nys_approx_->block_sz_;
        //init matrices dimensions
        maxIter_ = std::ceil((double)A.rows()/(double)block_sz_); //A.rows() is equal to A.cols()
        auto maxMatDim = maxIter_*block_sz_; //A.rows() is equal to A.cols()
        //The Krylov subspace is used to approximate the range of A, it is stored in Q_
        X_.resize(A_.rows(), maxMatDim);
        //Residual matrix
        T_.resize(A.rows(), maxMatDim);
        //Init T_ and shift_ parameter
        T_.leftCols(block_sz_) = GaussianMatrix(A.rows(), block_sz_, nys_approx_->seed_);
        shift_ = DOUBLE_TOLERANCE*A_.trace();
        //stopping criterion
        norm_A_ = A_.norm();
    }
    bool stop(){
        if(index_==0) return false;
        double squared_reconstr_err = std::pow(norm_A_,2)-std::pow((X_.leftCols(index_*block_sz_)*T_.leftCols(index_*block_sz_).transpose()).norm(),2);
        return squared_reconstr_err < std::pow(norm_A_*nys_approx_->tolerance_,2) || index_>=maxIter_;
    }
    NysIterations& operator++(){

      // another block grahm schimit here
      
        Eigen::HouseholderQR<DMatrix<double>> qr;
        X_.middleCols(index_*block_sz_,block_sz_) = T_.middleCols(std::max(index_-1,0)*block_sz_,block_sz_);
        //Block Gram-Schmidt
        X_.middleCols(index_*block_sz_,block_sz_) =
                (DMatrix<double>::Identity(A_.rows(),A_.rows()) - X_.leftCols(index_*block_sz_) * X_.leftCols(index_*block_sz_).transpose()) * X_.middleCols(index_*block_sz_,block_sz_);
        X_.middleCols(index_*block_sz_,block_sz_) =
                (DMatrix<double>::Identity(A_.rows(),A_.rows()) - X_.leftCols(index_*block_sz_) * X_.leftCols(index_*block_sz_).transpose()) * X_.middleCols(index_*block_sz_,block_sz_);
        qr.compute(X_.middleCols(index_*block_sz_, block_sz_));
        X_.middleCols(index_*block_sz_, block_sz_) = qr.householderQ() * DMatrix<double>::Identity(A_.rows(),block_sz_);
        //Updating the residual matrix
        T_.middleCols(index_*block_sz_,block_sz_) = (A_+shift_*DMatrix<double>::Identity(A_.rows(),A_.cols()))*X_.middleCols(index_*block_sz_, block_sz_);

        index_++;
        return *this;
    }
    ~NysIterations(){
        Eigen::LLT<DMatrix<double>> chol = Eigen::LLT<DMatrix<double>>(X_.leftCols(index_*block_sz_).transpose()*T_.leftCols(index_*block_sz_));
        nys_approx_->Z_ = chol.matrixU().solve<Eigen::OnTheRight>(T_.leftCols(index_*block_sz_));
        nys_approx_->shift_ = shift_; //passing the shift parameter (computed on A)
    }
};


/*
//NysSI: Nystrom Subspace Iterations
template<typename MatrixType, StoppingPolicy StopPolicy>
class NysIterations<MatrixType,IterationPolicy::SubspaceIterations,StopPolicy>{
private:
    int index_,maxIter_;
    DMatrix<double> A_; //matrix to be decomposed
    //Storage of the decomposition
    DMatrix<double> X_; //range
    DMatrix<double> T_; //residual matrix
    //Structured storage
    DMatrix<double> Z_; //factor: A + shift*I = Z*Z^T, alternative storage
    //Method parameters
    size_t block_sz_;
    double shift_;
    //Pointer to the calling Range Finder
    NystromApproximator<MatrixType,IterationPolicy::SubspaceIterations,StopPolicy> *nys_approx_;
public:
    NysIterations(const MatrixType &A, NystromApproximator<MatrixType,IterationPolicy::SubspaceIterations,StopPolicy>* nys_approx) :
    A_(A),nys_approx_(nys_approx){
        index_=0;
        block_sz_ = nys_approx_->block_sz_;
        //Init X_ and shift_ parameter
        T_ = GaussianMatrix(A.cols(), block_sz_, nys_approx_->seed_);
        shift_ = std::numeric_limits<float>::denorm_min()*A_.trace();
    }
    bool stop(){
        if(index_==0) return false;
        if constexpr(StopPolicy==StoppingPolicy::ReconstructionAccuracy){
            return (A_-X_*T_.transpose()).norm() < nys_approx_->tolerance_;
        }
    }
    NysIterations& operator++(){
        this->index_++;
        Eigen::HouseholderQR<decltype(T_)> qr(T_);
        X_ = qr.householderQ()*DMatrix<double>::Identity(A_.cols(), block_sz_);
        T_ = A_ * X_;
        return *this;
    }
    ~NysIterations(){
        T_ = T_ + shift_*X_;
        auto C = Eigen::LLT<DMatrix<double>>(X_.transpose()*T_).matrixL();
        nys_approx_->Z_ = C.solve<Eigen::OnTheRight>(T_);
        nys_approx_->shift_ = shift_; //passing the shift parameter (computed on A)
    }
};
*/

}//core
}//fdapde

#endif //NYSTROM_APPROXIMATION_H
