#ifndef __ELASTIC_MATH_FUNCTION_H__
#define __ELASTIC_MATH_FUNCTION_H__

#include <memory>
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

class SumEnergy : public Energy {
public :
    typedef std::shared_ptr<std::vector<std::shared_ptr<Energy>>> EnergyBuffer;
    SumEnergy(EnergyBuffer &pool);
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;
private :
    const EnergyBuffer pool_;
};

}}
#endif
