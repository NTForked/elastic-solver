#include "editing_energy.h"

#include <iostream>

#include "util.h"

using namespace std;
using namespace Eigen;

namespace cj { namespace elastic {

ReducedEditEnergy::ReducedEditEnergy(const size_t nbr_frame,
                                     const size_t reduced_dim,
                                     const double time_step,
                                     const double alpha,
                                     const double beta,
                                     const VectorXd &lambda)
    : nbr_frame_(nbr_frame),
      reduced_dim_(reduced_dim),
      h_(time_step),
      alpha_(alpha),
      beta_(beta),
      Lam_(lambda) {
    assert(lambda.rows() == reduced_dim_);
    Id_.setOnes(reduced_dim_, 1);
    Damp_ = alpha_*Id_ + beta_*Lam_;
}

size_t ReducedEditEnergy::Nx() const {
    return reduced_dim_ * nbr_frame_;
}

int ReducedEditEnergy::Val(const double *x, double *val) const {
    Map<const MatrixXd> X(x, reduced_dim_, nbr_frame_);
    for (size_t i = 1; i <= nbr_frame_-2; ++i) {
        *val += 0.5*h_*( 1.0/(h_*h_)*(X.col(i+1)-2*X.col(i)+X.col(i-1))
                         +1.0/h_*Damp_.asDiagonal()*(X.col(i+1)-X.col(i))
                         +Lam_.asDiagonal()*X.col(i)
                       ).squaredNorm();
    }
    return 0;
}

int ReducedEditEnergy::Gra(const double *x, double *gra) const {
    Map<const MatrixXd> X(x, reduced_dim_, nbr_frame_);
    Map<MatrixXd> grad(gra, reduced_dim_, nbr_frame_);
    for (size_t i = 1; i <= nbr_frame_-2; ++i) {
        VectorXd gra = h_*( 1.0/(h_*h_)*(X.col(i+1)-2*X.col(i)+X.col(i-1))
                            +1.0/h_*Damp_.asDiagonal()*(X.col(i+1)-X.col(i))
                            +Lam_.asDiagonal()*X.col(i)
                          );
        grad.col(i+1) += (1.0/(h_*h_)*Id_+1.0/h_*Damp_).asDiagonal() * gra;
        grad.col(i) += (-2.0/(h_*h_)*Id_ - 1.0/h_*Damp_ + Lam_).asDiagonal() * gra;
        grad.col(i-1) += 1.0/(h_*h_) * gra;
    }
    return 0;
}

int ReducedEditEnergy::Hes(const double *x, SparseMatrix<double> *hes) const {
    vector<Triplet<double>> trips;
    for (size_t i = 1; i <= nbr_frame_-2; ++i) {
        // i-1, i -1, i, i+1
        // i, i-1, i, i+1
        //i+1, i-1, i, i+1
    }
    return 0;
}

}}
