#ifndef ELASTIC_UTIL_H
#define ELASTIC_UTIL_H

#include <unordered_set>
#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace cj { namespace elastic {

template <typename T>
int RemoveRowCol(const std::unordered_set<size_t> &idx,
                 Eigen::SparseMatrix<T> &A,
                 zjucad::matrix::matrix<size_t> &l2g) {
    using Eigen::Triplet;
    using zjucad::matrix::matrix;
    matrix<size_t> g2l(A.cols());
    l2g.resize(A.cols() - idx.size());
    for (size_t i = 0, j = 0; i < g2l.size(); ++i) {
        if ( idx.find(i) == idx.end() ) {
            g2l[i] = j;
            l2g[j] = i;
            ++j;
        } else {
            g2l[i] = -1;
        }
    }
    std::vector<Triplet<T>> trips;
    for (size_t j = 0; j < A.cols(); ++j) {
        for (size_t cnt = A.outerIndexPtr()[j];
             cnt < A.outerIndexPtr()[j + 1]; ++cnt) {
            size_t i = A.innerIndexPtr()[cnt];
            if ( g2l[i] != -1 && g2l[j] != -1 )
                trips.push_back(Triplet<T>(g2l[i], g2l[j], A.valuePtr()[cnt]));
        }
    }
    const size_t dim = l2g.size();
    A.resize(dim, dim);
    A.setZero();
    A.reserve(trips.size());
    A.setFromTriplets(trips.begin(), trips.end());
    return 0;
}


/// Generalized eigenvalue solver for
/// vibration modes: Kx = \lambda Mx
int SolveModalBasis(const Eigen::SparseMatrix<double> &K,
                    const Eigen::DiagonalMatrix<double, -1> &M,
                    const std::unordered_set<size_t> &fix,
                    const size_t nr,
                    Eigen::MatrixXd &U,
                    Eigen::VectorXd &lambda);

int SolveMassPCA();

}}
#endif
