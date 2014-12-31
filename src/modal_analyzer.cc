#include "modal_analyzer.h"

#include <iostream>
#include "arpaca.h"
#include "util.h"

using namespace std;
using namespace Eigen;

namespace cj { namespace elastic {

ModalAnalyzer::ModalAnalyzer(const SparseMatrix<double> &K,
                             const SparseMatrix<double> &M,
                             const size_t nbr_basis,
                             const unordered_set<size_t> &fix)
    : K_(K), M_(M), nbr_(nbr_basis) {
    // construct L^{-1}
    invL_.resize(M_.rows());
    for (size_t i = 0; i < invL_.rows(); ++i) {
        double diag = 1.0 / sqrt(M_.coeff(i, i));
        invL_.diagonal()[i] = diag;
    }
    // construct flag
    flag_.resize(M_.cols());
    size_t j = 0;
    for (size_t i = 0; i < flag_.size(); ++i) {
        if ( fix.find(i) != fix.end() )
            flag_[i] = -1;
        else
            flag_[i] = j++;
    }
}

int ModalAnalyzer::Compute() {
    // build $L^{-1}KL^{-1}$ and then
    // solve $L^{-1}KL^{-1}y = \lambda y$
    SparseMatrix<double> LHS = K_;
    for (size_t j = 0; j < LHS.cols(); ++j) {
        for (SparseMatrix<double>::InnerIterator it(LHS, j); it; ++it) {
            it.valueRef() *= invL_.diagonal()[it.row()];
            it.valueRef() *= invL_.diagonal()[it.col()];
        }
    }
    RemoveSparseRowCol(LHS, flag_);

    arpaca::SymmetricEigenSolver<double> sol
            = arpaca::Solve(LHS, nbr_, arpaca::MAGNITUDE_SMALLEST);
    printf("[INFO] arpack %d iter, %d converged, %s\n",
           sol.num_actual_iterations(), sol.num_converged_eigenvalues(), sol.GetInfo());
    modes_ = sol.eigenvectors();
    freqs_ = sol.eigenvalues();

    // extend modes
    MatrixXd U(K_.rows(), modes_.cols());
    for (size_t i = 0; i < U.rows(); ++i) {
        if ( flag_[i] == -1 )
            U.row(i).setZero();
        else
            U.row(i) = modes_.row(flag_[i]);
    }
    modes_ = U;
    /// here we can check $\psi_i^T*M*\psi_j=\delta_{ij}$
    /// and $\psi_i^T*K*\psi_j=\lambda_{i}\delta_{ij}$
    for (size_t i = 0; i < modes_.cols(); ++i) {
        modes_.col(i) = invL_ * modes_.col(i).eval();
    }
    return 0;
}

MatrixXd ModalAnalyzer::get_modes() const {
    return modes_;
}

MatrixXd& ModalAnalyzer::get_modes() {
    return modes_;
}

VectorXd ModalAnalyzer::get_freqs() const {
    return freqs_;
}

VectorXd& ModalAnalyzer::get_freqs() {
    return freqs_;
}

}}
