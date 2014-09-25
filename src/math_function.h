#ifndef __CJ_MATH_FUNCTION_H__
#define __CJ_MATH_FUNCTION_H__

#include <Eigen/Sparse>

class Energy {
public :
    Energy();
    virtual ~Energy();

    virtual size_t Nx() const = 0;
    virtual int Val(const double *x, double *val) const = 0;
    virtual int Gra(const double *x, double *gra) const = 0;
    virtual int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const = 0;
};

class Constraint {
public :
    Constraint();
    virtual ~Constraint();

    virtual size_t Nx() const = 0;
    virtual size_t Nf() const = 0;
    virtual int Val(const double *x, double *val) const = 0;
    virtual int Jac(const double *x, Eigen::SparseMatrix<double> *hes) const = 0;
};






#endif
