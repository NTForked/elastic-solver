#ifndef __CJ_MASS_MATRIX_H__
#define __CJ_MASS_MATRIX_H__

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

class MassMatrix {
public :
    MassMatrix(const zjucad::matrix::matrix<size_t> &tets,
               const zjucad::matrix::matrix<double> &nods,
               const double density);
    int ComputeLumpedMassMatrix(Eigen::SparseMatrix<double> &M);
    int Lump();
private :
    const double rho_;
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;
};

}}
#endif
