#include <iostream>
#include <unordered_set>
#include <boost/property_tree/ptree.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <zjucad/matrix/io.h>

#include "src/util.h"
#include "src/warping.h"
#include "src/constitutive.h"
#include "src/util.h"
#include "src/editing_energy.h"

#define CALL_SUB_PROG(prog)             \
    int prog(ptree &pt);                \
    if ( strcmp(argv[1], #prog) == 0 )  \
        return prog(pt);

using namespace std;
using boost::property_tree::ptree;
using namespace Eigen;
using namespace cj::elastic;
using namespace zjucad::matrix;

int test_inner_iteration(ptree &pt) {
    SparseMatrix<double> A(3, 3);
    A.setIdentity();
    for (size_t j = 0; j < A.cols(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, j); it; ++it) {
            printf("(%d, %d, %lf)\n", it.row(), it.col(), it.value());
        }
    }
    return 0;
}

int test_remove_sparse(ptree &pt)  {
    MatrixXd B(5, 5);
    B.setRandom();
    SparseMatrix<double> A = B.sparseView();
    cout << A << "\n\n";
    unordered_set<size_t> id;
    id.insert(0);
    id.insert(2);
    std::vector<size_t> g2l(A.cols());
    size_t j = 0;
    for (size_t i = 0; i <g2l.size(); ++i) {
        if ( id.find(i) != id.end() )
            g2l[i] = -1;
        else
            g2l[i] = j++;
    }
    cj::elastic::RemoveSparseRowCol(A, g2l);
    cout << A << "\n";
    return 0;
}

int test_matrix_exponential(ptree &pt) {
    Matrix3d A(3, 3);
    A.setIdentity();
    cout << A.exp() << "\n";
    return 0;
}

int test_axb_energy(ptree &pt) {
    double val = 0;
    double w = 1;
    double X[9] = {0, 1, 1, 1, 1, 0, 1, 1, 0};
    double A[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double B[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    axb_energy_(&val, X, A, B, &w);
    cout << val << "\n";
    return 0;
}

int test_make_matrix(ptree &pt) {
    const double arr[] = {100, 1, 2, 3, 4, 5, 6, 7, 8};
    Matrix3d SH = make_skew_symm<double>(&arr[0]);
    Matrix3d H = make_symm<double>(&arr[3]);
    cout << "skew: \n" << SH << "\n";
    cout << "symm: \n" << H << "\n";
    return 0;
}

int test_matrix_lvalue(ptree &pt) {
    matrix<double> A = rand(5, 1);
    cout << "origin:\n" << A << "\n";
    A(colon(3, 4)) += ones<double>(2, 1);
    cout << "after:\n" << A << "\n";
    return 0;
}

int test_matrix_plus(ptree &pt) {
    matrix<double> I = identity_matrix<double>(3);
    I += 0.5;
    cout << I << "\n";
    return 0;
}

int test_reduced_edit_energy(ptree &pt) {
    VectorXd temp(3);
    temp << 1, 2, 3;
    ReducedEditEnergy e(2, 3, 0.1, 0.0, 0.0, temp);
    return 0;
}

int test_diagonal_arithmetic(ptree &pt) {
    DiagonalMatrix<double, -1> a(5), b(5);
    a.diagonal() << 0, 2, 4, 6, 16;
    b.diagonal() << 2, 2, 2, 2, 2;
    cout << "sum:\n" << (a+b).diagonal() << endl;
    cout << "diff:\n" << (a-b).diagonal() << endl;
    cout << "product:\n" << (a*b).diagonal() << endl;
    cout << "quotient:\n" << (a/b).diagonal() << endl;
    return 0;
}

int main(int argc, char *argv[])
{
    ptree pt;
    CALL_SUB_PROG(test_inner_iteration);
    CALL_SUB_PROG(test_remove_sparse);
    CALL_SUB_PROG(test_matrix_exponential);
    CALL_SUB_PROG(test_axb_energy);
    CALL_SUB_PROG(test_make_matrix);
    CALL_SUB_PROG(test_matrix_lvalue);
    CALL_SUB_PROG(test_matrix_plus);
    CALL_SUB_PROG(test_reduced_edit_energy);
    CALL_SUB_PROG(test_diagonal_arithmetic);
    cout << "no sub program.\n";
    return 0;
}
