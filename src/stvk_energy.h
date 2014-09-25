#ifndef __CJ_STVK_ENERGY_H__
#define __CJ_STVK_ENERGY_H__

#include <zjucad/matrix/matrix.h>
#include "math_function.h"

class StVKEnergy : public Energy {
public :
    StVKEnergy();
    virtual ~StVKEnergy();

    virtual size_t Nx() const;
    virtual size_t Nf() const;
    virtual int Val(const double *x, double *val) const;
    virtual int Gra(const double *x, double *gra) const;
    virtual int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;

private :
    const zjucad::matrix::matrix<double> &tets;
    const zjucad::matrix::matrix<double> &nods;
};
















#endif
