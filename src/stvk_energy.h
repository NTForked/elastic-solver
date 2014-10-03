#ifndef __ELASTIC_STVK_ENERGY_H__
#define __ELASTIC_STVK_ENERGY_H__

#include <zjucad/matrix/matrix.h>
#include "math_function.h"

namespace cj { namespace elastic {

// \miu and \lambda are lame parameters
// V * (\miu * E : E + 0.5 * \lambda * tr^2(E))

class StVKEnergy : public Energy {
public :
    StVKEnergy(const zjucad::matrix::matrix<size_t> &tets,
               const zjucad::matrix::matrix<double> &nods,
               const double lambda,
               const double miu,
               const double w = 1);
    ~StVKEnergy() { }
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;

private :
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;
    const double lambda_;
    const double miu_;
    const double w_;
    zjucad::matrix::matrix<zjucad::matrix::matrix<double>> Dm_;
    zjucad::matrix::matrix<double> volume_;
};

}}

#endif
