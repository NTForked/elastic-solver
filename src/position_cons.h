#ifndef __CJ_POSITION_CONSTRAINT_H__
#define __CJ_POSITION_CONSTRAINT_H__

#include <zjucad/matrix/matrix.h>
#include "math_function.h"

namespace cj { namespace elastic {

class PositionCons : public Constraint {
public :
    PositionCons(const std::vector<size_t> &idx,
                 const zjucad::matrix::matrix<double> &uc,
                 const double w = 1);
    ~PositionCons() { }
    size_t Nx() const;
    size_t Nf() const;
    int Val(const double *x, double *val) const;
    int Jac(const double *x, Eigen::SparseMatrix<double> *jac) const;
private :
    const std::vector<size_t> idx_;
    const zjucad::matrix::matrix<double> uc_;
    const double w_;
};

}}
#endif
