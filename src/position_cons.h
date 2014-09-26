#ifndef __CJ_POSITION_CONSTRAINT_H__
#define __CJ_POSITION_CONSTRAINT_H__

#include <zjucad/matrix/matrix.h>
#include "math_function.h"

class  PositionCons : public Constraint {
public :
    PositionCons(const std::vector<size_t> &idx, const double w = 1);
    virtual ~PositionCons();
    virtual size_t Nx() const;
    virtual size_t Nf() const;
    virtual int Val(const double *x, double *val) const;
    virtual int Jac(const double *x, Eigen::SparseMatrix<double> *hes) const;
private :
    const std::vector<size_t> idx_;
    zjucad::matrix::matrix<double> uc_;
    const double w_;
};

#endif
