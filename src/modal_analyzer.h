#ifndef MODEL_ANALYZER_H
#define MODEL_ANALYZER_H

#include <Eigen/Sparse>

namespace cj { namespace elastic {

class ModalAnalyzer {
public:
    ModalAnalyzer(const Eigen::SparseMatrix<double> &K,
                  const Eigen::SparseMatrix<double> &M,
                  const size_t nbr_basis);
    int Compute();
    Eigen::MatrixXd get_modes();
    const Eigen::MatrixXd& get_modes() const;
    Eigen::VectorXd get_freqs();
    const Eigen::VectorXd& get_freqs() const;
    void Display();
private :
    const Eigen::SparseMatrix<double> &K_, &M_;
    const size_t nbr_;
    Eigen::DiagonalMatrix<double, -1> invL_;
    Eigen::MatrixXd modes_;
    Eigen::VectorXd freqs_;
};

}}
#endif
