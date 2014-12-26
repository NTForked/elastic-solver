#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include "src/arpaca.h"
#include "src/util.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

#ifdef TEST_ARPACA
int main(int argc, char *argv[])
{
    if ( argc != 3 ) {
        cerr << "# Usage: ./prog rows cols\n";
        return __LINE__;
    }
    MatrixXd A(atoi(argv[1]), atoi(argv[2]));
    A.setOnes();
    MatrixXd ATA = A.transpose() * A;
    SparseMatrix<double> H = ATA.sparseView();

    SelfAdjointEigenSolver<MatrixXd> sl(ATA);
    cout << sl.eigenvalues().head(20).transpose() << endl << endl;

    arpaca::SymmetricEigenSolver<double>
            sol = arpaca::Solve(H, 20, arpaca::ALGEBRAIC_SMALLEST);
    cout << sol.eigenvalues().transpose() << endl << endl;

    cout << "done\n";
    return 0;
}
#else
int main(int argc, char *argv[])
{
//    Matrix4d A;
//    A.setRandom();
//    SparseMatrix<double> B = A.sparseView();
//    cout << B << endl << endl;

//    std::unordered_set<size_t> idx;
//    idx.insert(1);
//    idx.insert(2);
//    idx.insert(3);

//    matrix<size_t> l2g;
//    cj::elastic::RemoveRowCol(idx, B, l2g);

//    cout << B << endl << endl;
//    cout << "done\n" << endl;
    return 0;
}
#endif
