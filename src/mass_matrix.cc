#include "mass_matrix.h"

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace cj { namespace elastic {

MassMatrix::MassMatrix(const matrix<size_t> &tets,
                       const matrix<double> &nods,
                       const double density)
    : tets_(tets), nods_(nods), rho_(density) {
}

void MassMatrix::ComputeCompactMass() {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> tet_nods = nods_(colon(), tets_(colon(), i));
        matrix<double> e = tet_nods(colon(), colon(1, 3))
                         - tet_nods(colon(), 0) * ones<double>(1, 3);
        double dm = rho_ * fabs(det(e)) / 6.0 / 20.0;
        for (size_t p = 0; p < 4; ++p) {
            for (size_t q = p; q < 4; ++q) {
                trips.push_back(Triplet<double>(tets_(p, i), tets_(q, i), dm));
                trips.push_back(Triplet<double>(tets_(q, i), tets_(p, i), dm));
            }
        }
    }
    compact_mass_mat.resize(nods_.size(2), nods_.size(2));
    compact_mass_mat.reserve(trips.size());
    compact_mass_mat.setFromTriplets(trips.begin(), trips.end());
    compact_mass_mat.makeCompressed();
}

int MassMatrix::MassMatrix::Compute(SparseMatrix<double> &M, bool lumped) {
    ComputeCompactMass();
    vector<Triplet<double>> trips;
    size_t dim = compact_mass_mat.cols();
    for (size_t j = 0; j < dim; ++j) {
        for (size_t cnt = compact_mass_mat.outerIndexPtr()[j];
             cnt < compact_mass_mat.outerIndexPtr()[j + 1]; ++cnt) {
            const size_t I = 3 * compact_mass_mat.innerIndexPtr()[cnt];
            const size_t J = 3 * j;
            const double V = compact_mass_mat.valuePtr()[cnt];
            trips.push_back(Triplet<double>(I + 0, J + 0, V));
            trips.push_back(Triplet<double>(I + 1, J + 1, V));
            trips.push_back(Triplet<double>(I + 2, J + 2, V));
        }
    }
    M.resize(dim * 3, dim * 3);
    M.reserve(trips.size());
    M.setFromTriplets(trips.begin(), trips.end());
    if ( lumped )
        Lump(M);
    M.makeCompressed();
}

void MassMatrix::Lump(SparseMatrix<double> &M) {

}

}}
