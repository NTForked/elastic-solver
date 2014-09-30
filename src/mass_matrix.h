#ifndef __CJ_MASS_MATRIX_H__
#define __CJ_MASS_MATRIX_H__

#include <Eigen/Sparse>

class MassMatrix {
public :
    MassMatrix();
    int compute(Eigen::SparseMatrix<double> &M);
    int lump();
private :
};






#endif
