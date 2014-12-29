#ifndef MODEL_ANALYZER_H
#define MODEL_ANALYZER_H

#include <unordered_set>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

class ModalAnalyzer {
public:
    ModalAnalyzer(const Eigen::SparseMatrix<double> &K,
                  const Eigen::SparseMatrix<double> &M,
                  const size_t nbr_basis,
                  const std::unordered_set<size_t> &fix);
    int Compute();
    Eigen::MatrixXd get_modes() const;
    Eigen::MatrixXd& get_modes();
    Eigen::VectorXd get_freqs() const;  // rvalue
    Eigen::VectorXd& get_freqs();       // lvalue
private :
    const Eigen::SparseMatrix<double> &K_, &M_;
    const size_t nbr_;
    std::vector<size_t> flag_;
    Eigen::DiagonalMatrix<double, -1> invL_;
    Eigen::MatrixXd modes_;
    Eigen::VectorXd freqs_;
};

}}
#endif
