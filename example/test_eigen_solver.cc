#include <iostream>
#include <Eigen/Sparse>

#include "src/arpaca.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
    SparseMatrix<double> A, B;
    A.resize(20, 20);
    B.resize(20, 20);
    A.setIdentity();
    B.setIdentity();
    arpaca::SymmetricEigenSolver<double>
            _sol = arpaca::Solve(A, 10, arpaca::MAGNITUDE_SMALLEST);
    cout << _sol.eigenvalues().transpose() << "\n";
    cout << "done\n";
    return 0;
}
