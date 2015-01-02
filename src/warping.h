#ifndef WARPING_ENERGY_H
#define WARPING_ENERGY_H

#include <zjucad/matrix/matrix.h>

#include "math_function.h"

extern "C" {

void axb_energy_(double *val,
                 const double *X,
                 const double *A,
                 const double *B,
                 const double *w);

void axb_energy_jac_(double *gra,
                     const double *X,
                     const double *A,
                     const double *B,
                     const double *w);

void axb_energy_hes_(double *hes,
                     const double *X,
                     const double *A,
                     const double *B,
                     const double *w);
}

namespace cj { namespace elastic {

class WarpingEnergy : public Energy {
public:
    typedef zjucad::matrix::matrix<size_t> matrix_t;
    typedef zjucad::matrix::matrix<double> matrix_d;
    WarpingEnergy(const matrix_t &tets,
                  const matrix_d &nods,
                  const matrix_d &G,
                  const matrix_d &tetRS,
                  const matrix_d &vol);
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;
private:
    const matrix_t &tets_;
    const matrix_d &nods_;
    const matrix_d &G_;
    const matrix_d &tetRS_;
    const matrix_d &vol_;
};

class FixMassCenter : public Constraint {
public:
    typedef zjucad::matrix::matrix<size_t> matrix_t;
    typedef zjucad::matrix::matrix<double> matrix_d;
    FixMassCenter(const matrix_t &tets,
                  const matrix_d &nods,
                  const matrix_d &vols);
    size_t Nx() const;
    size_t Nf() const;
    int Val(const double *x, double *val) const;
    int Jac(const double *x, Eigen::SparseMatrix<double> *jac) const;
private :
    const matrix_t &tets_;
    const matrix_d &nods_;
    const matrix_d &vols_;
    matrix_d wgt_;
};

}}
#endif
