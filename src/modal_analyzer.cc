//#include "modal_analyzer.h"

//#include <iostream>
//#include "arpaca.h"
//#include "ArpackSelfAdjointEigenSolver.h"

//using namespace std;
//using namespace Eigen;

//namespace cj { namespace elastic {

//ModalAnalyzer::ModalAnalyzer(const SparseMatrix<double> &K,
//                             const SparseMatrix<double> &M,
//                             const size_t nbr_basis)
//    : K_(K), M_(M), nbr_(nbr_basis) { }

//int ModalAnalyzer::Compute() {
//    ArpackGeneralizedSelfAdjointEigenSolver<SparseMatrix<double>>
//            sol(K_, M_, nbr_, "SM");
//    modes_ = sol.eigenvectors();
//    freqs_ = sol.eigenvalues();
//    return 0;
//}

//MatrixXd ModalAnalyzer::get_modes() {
//    return modes_;
//}

//const MatrixXd& ModalAnalyzer::get_modes() const {
//    return modes_;
//}

//VectorXd ModalAnalyzer::get_freqs() {
//    return freqs_;
//}

//const VectorXd& ModalAnalyzer::get_freqs() const {
//    return freqs_;
//}

//}}
