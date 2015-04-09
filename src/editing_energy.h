#ifndef EDITING_ENERGY_H
#define EDITING_ENERGY_H

#include "math_function.h"

namespace cj { namespace elastic {

class ReducedEditEnergy : Energy
{
public:
    ReducedEditEnergy(const size_t nbr_frame,
                      const size_t reduced_dim,
                      const double time_step,
                      const double alpha,
                      const double beta,
                      const Eigen::VectorXd &lambda);
    size_t Nx() const;
    int Val(const double *x, double *val) const;
    int Gra(const double *x, double *gra) const;
    int Hes(const double *x, Eigen::SparseMatrix<double> *hes) const;
private:
    const size_t nbr_frame_;
    const size_t reduced_dim_;
    const double h_;
    const double alpha_;
    const double beta_;
    Eigen::VectorXd lambda_;
    Eigen::DiagonalMatrix<double, -1> Lam_;
    Eigen::DiagonalMatrix<double, -1> Id_;
};

class Manipulation
{

};

}}
#endif
