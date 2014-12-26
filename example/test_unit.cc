#include <iostream>
#include <boost/property_tree/ptree.hpp>

#include "src/util.h"

#define CALL_SUB_PROG(prog)          \
    int prog(ptree &pt);                             \
    if ( strcmp(argv[1], #prog) == 0 )  \
        return prog(pt);

using namespace std;
using boost::property_tree::ptree;
using namespace Eigen;

int test_inner_iteration(ptree &pt) {
    SparseMatrix<double> A(10, 10);
    A.setIdentity();
    for (size_t j = 0; j < A.cols(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, j); it; ++it) {
            printf("(%d, %d, %lf\n", it.row(), it.col(), it.value());
        }
    }
    return 0;
}

int test_remove_sparse(ptree &pt) {
    return 0;
}

int main(int argc, char *argv[])
{
    ptree pt;
    CALL_SUB_PROG(test_inner_iteration);
    CALL_SUB_PROG(test_remove_sparse);
    cout << "no sub program.\n";
    return 0;
}
