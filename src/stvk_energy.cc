#include "stvk_energy.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>
#include "stvk_tet.h"

using namespace std;
using namespace zjucad::matrix;

StVKEnergy::StVKEnergy(const zjucad::matrix::matrix<size_t> &tets,
                       const zjucad::matrix::matrix<double> &nods,
                       const double lambda,
                       const double miu,
                       const double w)
    : tets_(tets), nods_(nods), lambda_(lambda), miu_(miu), w_(w) {
    Dm_.resize(tets_.size(2));
    volume_.resize(tets_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> tet_nods = nods_(colon(), tets_(colon(), i));
        Dm_[i] = tet_nods(colon(), colon(1, 3)) - tet_nods(colon(), 0) * ones<double>(1, 3);
        volume_[i] = fabs(det(Dm_[i])) / 6.0;
        if ( inv(Dm_[i]) ) {
            cerr << "[INFO] degenerated tet.\n";
        }
    }
}

size_t StVKEnergy::Nx() const {
    return nods_.size();
}

int StVKEnergy::Val(const double *x, double *val) const {
    itr_matrix<const double *> dx(3, nods_.size(2), x);
    matrix<double> nods = nods_  + dx;
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
        double v = 0;
        stvk_tet_(&v, &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
        *val += v;
    }
    *val *= w_;
    return 0;
}

int StVKEnergy::Gra(const double *x, double *gra) const {
    itr_matrix<const double *> dx(3, nods_.size(2), x);
    matrix<double> nods = nods_ + dx;
    itr_matrix<double *> g(nods_.size(), 1, gra);
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
        matrix<double> g_(12);
        stvk_tet_jac_(&g_[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
        for (size_t k = 0; k < 12; ++k)
            g[tets_(k / 3, i) * 3 + k % 3] += g_[k];
    }
    g *= w_;
    return 0;
}

int StVKEnergy::Hes(const double *x, Eigen::SparseMatrix<double> *hes) const  {
    itr_matrix<const double *> dx(3, nods_.size(2), x);
    matrix<double> nods = nods_ + dx;
    vector<Eigen::Triplet<double>> trips;
    for (size_t i = 0;  i < tets_.size(2); ++i) {
        matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
        matrix<double> H(12, 12);
        stvk_tet_hes_(&H[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
        for (size_t p = 0; p < 12; ++p) {
            for (size_t q = 0; q < 12; ++q) {
                size_t I = tets_(p / 3, i) * 3 + p % 3;
                size_t J = tets_(q / 3, i) * 3 + q % 3;
                trips.push_back(Eigen::Triplet<double>(I, J, w_ * H(p, q)));
            }
        }
    }
    hes->resize(nods_.size(), nods_.size());
    hes->reserve(trips.size());
    hes->setFromTriplets(trips.begin(), trips.end());
    hes->makeCompressed();
    return 0;
}
