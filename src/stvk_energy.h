#ifndef __CJ_STVK_ENERGY_H__
#define __CJ_STVK_ENERGY_H__

#include <zjucad/matrix/matrix.h>
#include "math_function.h"

class StVKEnergy : public Energy {
public :
    StVKEnergy(const zjucad::matrix::matrix<size_t>  &tets,
                            const zjucad::matrix::matrix<double> &nods,
                            const double lambda,
                            const double miu,
                            const double w = 1);
    virtual ~StVKEnergy();
    virtual size_t Nx() const;
    virtual int Val(const double *x, double *val) const;
    virtual int Gra(const double *x, double *gra) const;
    virtual int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;

private :
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;
    const double lambda_;
    const double miu_;
    const double w_;
    zjucad::matrix::matrix<zjucad::matrix::matrix<double>> Dm_;
    zjucad::matrix::matrix<double> volume_;
};

#endif
