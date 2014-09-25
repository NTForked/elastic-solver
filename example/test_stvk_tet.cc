#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/lapack.h>
#include <Eigen/Dense>
#include "src/stvk_tet.h"

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
    matrix<double> nods(3, 4), D(3, 3);
    nods = zeros<double>(3, 4);
    nods(colon(), colon(1, 3)) = identity_matrix<double>(3);

    cerr << "nodes: \n" << nods << endl;

    D = nods(colon(), colon(1, 3)) - nods(colon(), 0) * ones<double>(1, 3);
    double volume = fabs(det(D)) / 6.0;

    cerr << "\nD: \n" << D << endl;

    if ( inv(D) ) {
        cerr << "[INFO] inv fail!\n";
        return __LINE__;
    }

    matrix<double> rand_disp = rand(3, 1);
    nods += rand_disp * ones<double>(1, 3);
    double val = 0;
    double lambda = 0.1, miu = 0.2;
    stvk_tet_(&val, &nods[0], &D[0], &volume, &lambda, &miu);

    cout << "[INFO] tet energy value: " << val << endl;
    return 0;
}
