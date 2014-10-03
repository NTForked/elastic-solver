#ifndef __ELASTIC_MASS_MATRIX_H__
#define __ELASTIC_MASS_MATRIX_H__

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

class MassMatrix {
public :
    MassMatrix(const zjucad::matrix::matrix<size_t> &tets,
               const zjucad::matrix::matrix<double> &nods,
               const double density);
    int Compute(Eigen::SparseMatrix<double> &M, bool lumped = false);
private :
    void Lump(Eigen::SparseMatrix<double> &M);
    void ComputeCompactMass();

    const double rho_;
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;
    Eigen::SparseMatrix<double> compact_mass_mat;
};

}}
#endif
