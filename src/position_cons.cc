#include "position_cons.h"

#include <zjucad/matrix/itr_matrix.h>

using namespace std;
using namespace zjucad::matrix;

PositionCons::PositionCons(const vector<size_t> &idx,
                           const matrix<double> &uc,
                           const double w)
    : idx_(idx), uc_(uc), w_(w) { }

size_t PositionCons::Nx() const {
    return uc_.size();
}

size_t PositionCons::Nf() const {
    return idx_.size() * 3;
}

int PositionCons::Val(const double *x, double *val) const {
    itr_matrix<const double *> dx(3, uc_.size(2), x);
    itr_matrix<double *> v(3 * idx_.size(), 1, val);
    for (size_t i = 0; i < idx_.size(); ++i) {
        v(colon(3 * i, 3 * i + 2)) = w_ * (dx(colon(), idx_[i]) - uc_(colon(), idx_[i]));
    }
    return 0;
}

int PositionCons::Jac(const double *x, Eigen::SparseMatrix<double> *jac) const {
    vector<Eigen::Triplet<double>> trips;
    for (size_t i = 0; i < idx_.size(); ++i) {
        trips.push_back(Eigen::Triplet<double>(3 * i + 0, 3 * idx_[i] + 0, w_));
        trips.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * idx_[i] + 1, w_));
        trips.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * idx_[i] + 2, w_));
    }
    jac->resize(3 * idx_.size(), uc_.size());
    jac->reserve(trips.size());
    jac->setFromTriplets(trips.begin(), trips.end());
    jac->makeCompressed();
    return 0;
}

