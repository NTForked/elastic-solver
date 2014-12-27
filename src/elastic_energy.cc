#include "elastic_energy.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/lapack.h>
#include "types.h"
#include "math_function.h"
#include "constitutive.h"

using namespace zjucad::matrix;

namespace cj { namespace elastic {

// \miu and \lambda are lame parameters

class MaterialEnergy : public Energy {
public :
    MaterialEnergy(const zjucad::matrix::matrix<size_t> &tets,
                   const zjucad::matrix::matrix<double> &nods,
                   const double lambda,
                   const double miu,
                   const MaterialType &type,
                   const double w = 1)
        : tets_(tets), nods_(nods), lambda_(lambda), miu_(miu), type_(type), w_(w) {
        // precompute volume and [x1 - x0, x2 - x0, x3 - x0]^-1
        Dm_.resize(tets_.size(2));
        volume_.resize(tets_.size(2));
#pragma omp parallel for
        for (size_t i = 0; i < tets_.size(2); ++i) {
            matrix<double> tet_nods = nods_(colon(), tets_(colon(), i));
            Dm_[i] = tet_nods(colon(), colon(1, 3)) - tet_nods(colon(), 0) * ones<double>(1, 3);
            matrix<double> cache = Dm_[i];
            volume_[i] = fabs(det(cache)) / 6.0;   // det will change parameter because of in place LU decomposition
            if ( inv(Dm_[i]) ) {
                std::cerr << "[INFO] degenerated tet.\n";
            }
        }
        //-------set lambda1 and miu1 for Fung model------//
        const double E = 5 * 1e4;
        const double v = 0.45;
        lambda1_ = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
        miu1_ = E / (2.0 * (1.0 + v));
        c_ = 10;
        //------------------------------------------------//
    }

    ~MaterialEnergy() { }

    size_t Nx() const {
        return nods_.size();
    }

    int Val(const double *x, double *val) const {
        itr_matrix<const double *> dx(3, nods_.size(2), x);
        matrix<double> nods = nods_ + dx;
        for (size_t i = 0; i < tets_.size(2); ++i) {
            matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
            double v = 0;
            switch ( type_ ) {
                case MaterialType::LINEAR : {
                    linear_elastic_tet_(&v, &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::STVK : {
                    stvk_tet_(&v, &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::COROTATED : {
                    break;
                }
                case MaterialType::NEOHOOKEAN : {
                    neohookean_tet_(&v, &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::FUNG : {
                    fung_tet_(&v, &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_, &lambda1_, &miu1_, &c_);
                    break;
                }
            }
            *val += v;
        }
        *val *= w_;
        return 0;
    }

    int Gra(const double *x, double *gra) const {
        itr_matrix<const double *> dx(3, nods_.size(2), x);
        matrix<double> nods = nods_ + dx;
        itr_matrix<double *> g(nods_.size(), 1, gra);
#pragma omp parallel for
        for (size_t i = 0; i < tets_.size(2); ++i) {
            matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
            matrix<double> g_(12);
            switch ( type_ ) {
                case MaterialType::LINEAR : {
                    linear_elastic_tet_jac_(&g_[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::STVK : {
                    stvk_tet_jac_(&g_[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::COROTATED : {
                    break;
                }
                case MaterialType::NEOHOOKEAN : {
                    neohookean_tet_jac_(&g_[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::FUNG : {
                    fung_tet_jac_(&g_[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_, &lambda1_, &miu1_, &c_);
                    break;
                }
            }
#pragma omp critical
            for (size_t k = 0; k < 12; ++k)
                g[tets_(k / 3, i) * 3 + k % 3] += g_[k];
        }
        g *= w_;
        return 0;
    }

    int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const {
        itr_matrix<const double *> dx(3, nods_.size(2), x);
        matrix<double> nods = nods_ + dx;
        std::vector<Eigen::Triplet<double>> trips;
#pragma omp parallel for
        for (size_t i = 0;  i < tets_.size(2); ++i) {
            matrix<double> tet_nods = nods(colon(), tets_(colon(), i));
            matrix<double> H(12, 12);
            switch ( type_ ) {
                case MaterialType::LINEAR : {
                    linear_elastic_tet_hes_(&H[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::STVK : {
                    stvk_tet_hes_(&H[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::COROTATED : {
                    break;
                }
                case MaterialType::NEOHOOKEAN : {
                    neohookean_tet_hes_(&H[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_);
                    break;
                }
                case MaterialType::FUNG : {
                    fung_tet_hes_(&H[0], &tet_nods[0], &Dm_[i][0], &volume_[i], &lambda_, &miu_, &lambda1_, &miu1_, &c_);
                    break;
                }
            }
#pragma omp critical
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

private :
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;
    MaterialType type_;
    const double lambda_;
    const double miu_;
    const double w_;
    zjucad::matrix::matrix<zjucad::matrix::matrix<double>> Dm_;
    zjucad::matrix::matrix<double> volume_;
    double lambda1_;
    double miu1_;
    double c_;
};

Energy* BuildElasticEnergy(const matrix<size_t> &tets,
                           const matrix<double> &nods,
                           const double lambda,
                           const double miu,
                           const std::string &type,
                           const double w) {
    std::map<std::string, MaterialType> option;
    option["linear"]     = MaterialType::LINEAR;
    option["stvk"]       = MaterialType::STVK;
    option["corotated"]  = MaterialType::COROTATED;
    option["neohookean"] = MaterialType::NEOHOOKEAN;
    option["fung"]       = MaterialType::FUNG;

    auto it = option.find(type);
    if ( it != option.end() ) {
        std::cerr << "[INFO] " << type << " elastic model.\n";
        return new MaterialEnergy(tets, nods, lambda, miu, it->second, w);
    } else {
        return nullptr;
    }
}

}}
