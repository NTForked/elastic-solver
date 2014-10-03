#ifndef __ELASTIC_MATH_FUNCTION_H__
#define __ELASTIC_MATH_FUNCTION_H__

#include <Eigen/Sparse>


namespace cj { namespace elastic {

// interface of energy functional
class Energy {
public :
    virtual ~Energy() { }
    virtual size_t Nx() const = 0;
    virtual int Val(const double *x, double *val) const = 0;
    virtual int Gra(const double *x, double *gra) const = 0;
    virtual int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const = 0;
};

// interface of constraint map
class Constraint {
public :
    virtual ~Constraint() { }
    virtual size_t Nx() const = 0;
    virtual size_t Nf() const = 0;
    virtual int Val(const double *x, double *val) const = 0;
    virtual int Jac(const double *x, Eigen::SparseMatrix<double> *jac) const = 0;
};

}}
#endif
