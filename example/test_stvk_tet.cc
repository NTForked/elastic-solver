#include <iostream>
#include <fstream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "src/stvk_tet.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
    matrix<size_t> tets = colon(0, 3);
    matrix<double> nods = zeros<double>(3, 4);
    nods(colon(), colon(1, 3)) = identity_matrix<double>(3);
    {
        ofstream os("../../result/rest.vtk");
        tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
    }
    cerr << "nodes: \n" << nods << endl;

    matrix<double> D = nods(colon(), colon(1, 3)) - nods(colon(), 0) * ones<double>(1, 3);
    double volume = fabs(det(D)) / 6.0;

    if ( inv(D) ) {
        cerr << "[INFO] inv fail!\n";
        return __LINE__;
    }

    // apply transformation
    matrix<double> t = rand(3, 1);
    nods += t * ones<double>(1, 3);
    {
        std::ofstream os("../../result/transform.vtk");
        tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
    }

    // apply a rotation
    const double PI = 3.14159265358979323;
    const double angle = PI / 3;
    matrix<double> R(3, 3);
    R(0, 0) = 1; R(1, 1) = R(2, 2) = cos(angle);
    R(1, 2) = -sin(angle); R(2, 1) = sin(angle);
    nods = temp(R * nods);
    {
        std::ofstream os("../../result/rotation.vtk");
        tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
    }

    // apply scale
//    nods(colon(), 2) *= 10;

    double val = 0;
    matrix<double> gra = zeros<double>(12, 1);
    matrix<double> hes = zeros<double>(12, 12);
    double lambda = 0.1, miu = 0.2;
    stvk_tet_(&val, &nods[0], &D[0], &volume, &lambda, &miu);
    stvk_tet_jac_(&gra[0], &nods[0], &D[0], &volume, &lambda, &miu);
    stvk_tet_hes_(&hes[0], &nods[0], &D[0], &volume, &lambda, &miu);

    cout << "[INFO] H - HT: " << norm(hes - trans(hes)) << endl;
    cout << "[INFO] tet energy value: " << val << endl;
    return 0;
}
