#include "warping.h"

#include <zjucad/matrix/itr_matrix.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace cj { namespace elastic {

WarpingEnergy::WarpingEnergy(const matrix_t &tets,
                             const matrix_d &nods,
                             const matrix_d &G,
                             const matrix_d &tetRS,
                             const matrix_d &vol)
    : tets_(tets), nods_(nods), G_(G), tetRS_(tetRS), vol_(vol) { }

size_t WarpingEnergy::Nx() const {
    return nods_.size();
}

int WarpingEnergy::Val(const double *x, double *val) const {
    itr_matrix<const double*> X(3, Nx() / 3, x);
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> U = X(colon(), tets_(colon(), i));
        Matrix3d RS = make_skew_symm(&tetRS_(0, i)).exp()
                * (Matrix3d::Identity() + make_symm(&tetRS_(3, i)))
                - Matrix3d::Identity();
        double v = 0;
        axb_energy_(&v, &U[0], &G_(0, i), RS.data(), &vol_[i]);
        *val += v;
    }
    return 0;
}

int WarpingEnergy::Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, Nx() / 3, x);
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> U = X(colon(), tets_(colon(), i));
        Matrix3d RS = make_skew_symm(&tetRS_(0, i)).exp()
                * (Matrix3d::Identity() + make_symm(&tetRS_(3, i)))
                - Matrix3d::Identity();
        matrix<double> g(12);
        axb_energy_jac_(&g[0], &U[0], &G_(0, i), RS.data(), &vol_[i]);
        for (size_t j = 0; i < 12; ++j) {
            gra[3 * tets_(j / 3, i) + j % 3] += g[j];
        }
    }
    return 0;
}

int WarpingEnergy::Hes(const double *x, SparseMatrix<double> *hes) const {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> H(12, 12);
        axb_energy_hes_(&H[0], nullptr, &G_(0, i), nullptr, &vol_[i]);
        for (size_t p = 0; p < 12; ++p) {
            for (size_t q = 0; q < 12; ++q) {
                if ( H(p, q) != 0.0 ) {
                    size_t I = 3 * tets_(p / 3, i) + p % 3;
                    size_t J = 3 * tets_(q / 3, i) + q % 3;
                    trips.push_back(Triplet<double>(I, J, H(p, q)));
                }
            }
        }
    }
    hes->resize(Nx(), Nx());
    hes->reserve(trips.size());
    hes->setFromTriplets(trips.begin(), trips.end());
    return 0;
}

FixMassCenter::FixMassCenter(const matrix_t &tets,
                             const matrix_d &nods,
                             const matrix_d &vols)
    : tets_(tets), nods_(nods), vols_(vols) {
    wgt_ = zeros<double>(nods_.size(2), 1);
    for (size_t i = 0; i < vols_.size(2); ++i)
        wgt_(tets_(colon(), i)) += 0.25 * vols_[i] * ones<double>(4, 1);
    double frac = zjucad::matrix::sum(wgt_);
    wgt_ /= frac;
}

size_t FixMassCenter::Nx() const {
    return nods_.size();
}

size_t FixMassCenter::Nf() const {
    return 3;
}

int FixMassCenter::Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx() / 3, x);
    itr_matrix<double *> V(3, Nf() / 3, val);
    V += X * wgt_;
    return 0;
}

int FixMassCenter::Jac(const double *x, SparseMatrix<double> *jac) const {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < Nx() / 3; ++i) {
        trips.push_back(Triplet<double>(0, 3 * i + 0, wgt_[i]));
        trips.push_back(Triplet<double>(1, 3 * i + 1, wgt_[i]));
        trips.push_back(Triplet<double>(2, 3 * i + 2, wgt_[i]));
    }
    jac->resize(Nf(), Nx());
    jac->reserve(trips.size());
    jac->setFromTriplets(trips.begin(), trips.end());
    return 0;
}

}}
