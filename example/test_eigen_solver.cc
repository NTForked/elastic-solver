#include <iostream>
#include <Eigen/Sparse>

#include "src/ArpackSelfAdjointEigenSolver.h"
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
    ArpackGeneralizedSelfAdjointEigenSolver<SparseMatrix<double>>
            sol(A, B, 10, "SM");
    arpaca::SymmetricEigenSolver<double>
            _sol = arpaca::Solve(A, 10, arpaca::MAGNITUDE_SMALLEST);
    cout << sol.eigenvalues().transpose() << "\n";
    cout << _sol.eigenvalues().transpose() << "\n";
    cout << "done\n";
    return 0;
}
