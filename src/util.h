#ifndef ELASTIC_UTIL_H
#define ELASTIC_UTIL_H

#include <unordered_set>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

template <typename T>
int RemoveSparseRowCol(Eigen::SparseMatrix<T> &A,
                       const std::vector<size_t> g2l) {
    size_t n = 0;
    for (size_t i = 0; i < g2l.size(); ++i) {
        if ( g2l[i] != -1 )
            ++n;
    }
    std::vector<Eigen::Triplet<T>> trips;
    for (size_t j = 0; j < A.cols(); ++j) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it) {
            if ( g2l[it.row()] != -1 && g2l[it.col()] != -1 )
                trips.push_back(Eigen::Triplet<T>(g2l[it.row()], g2l[it.col()], it.value()));
        }
    }
    A.resize(n, n);
    A.reserve(trips.size());
    A.setFromTriplets(trips.begin(), trips.end());
    return 0;
}

template <class Mat>
bool isSymmetric(const Mat &A) {
    Mat AT = A.transpose();
    if ( (AT - A).squaredNorm() < 1e-20 ) {
        return true;
    }
    return false;
}

}}
#endif
