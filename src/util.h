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

template <typename T>
Eigen::Matrix<T, 3, 3> make_skew_symm(const T *element) {
    Eigen::Matrix<T, 3, 3> mat;
    mat.setZero();
    mat(1, 0) = element[0];
    mat(0, 1) = -element[0];
    mat(2, 0) = element[1];
    mat(0, 2) = -element[1];
    mat(2, 1) = element[2];
    mat(1, 2) = -element[2];
    return mat;
}

template <typename T>
Eigen::Matrix<T, 3, 3> make_symm(const T *element) {
    Eigen::Matrix<T, 3, 3> mat;
    mat(0, 0) = element[0];
    mat(1, 0) = element[1];
    mat(2, 0) = element[2];
    mat(0, 1) = mat(1, 0);
    mat(1, 1) = element[3];
    mat(2, 1) = element[4];
    mat(0, 2) = mat(2, 0);
    mat(1, 2) = mat(2, 1);
    mat(2, 2) = element[5];
    return mat;
}

// in next four overload arithmetic binary opertions,
// the operands require same size at run time. The
// function itself does not check for this fact, which
// needs to be guaranteed by users manually when invoking them

template <typename T>
inline Eigen::DiagonalMatrix<T, -1> operator +(const Eigen::DiagonalMatrix<T, -1> &a,
                                               const Eigen::DiagonalMatrix<T, -1> &b) {
    return Eigen::DiagonalMatrix<T, -1>(a.diagonal()+b.diagonal());
}

template <typename T>
inline Eigen::DiagonalMatrix<T, -1> operator -(const Eigen::DiagonalMatrix<T, -1> &a,
                                               const Eigen::DiagonalMatrix<T, -1> &b) {
    return Eigen::DiagonalMatrix<T, -1>(a.diagonal()-b.diagonal());
}

template <typename T>
inline Eigen::DiagonalMatrix<T, -1> operator *(const Eigen::DiagonalMatrix<T, -1> &a,
                                               const Eigen::DiagonalMatrix<T, -1> &b) {
    return Eigen::DiagonalMatrix<T, -1>(a.diagonal().cwiseProduct(b.diagonal()));
}

template <typename T>
inline Eigen::DiagonalMatrix<T, -1> operator /(const Eigen::DiagonalMatrix<T, -1> &a,
                                               const Eigen::DiagonalMatrix<T, -1> &b) {
    return Eigen::DiagonalMatrix<T, -1>(a.diagonal().cwiseQuotient(b.diagonal()));
}

}}
#endif
