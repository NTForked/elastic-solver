#ifndef ELASTIC_REDUCED2RS_H
#define ELASTIC_REDUCED2RS_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

extern "C" {

void disp2rs_(double       *val,
              const double *u,
              const double *G);

void disp2rs_jac_(double       *jac,
                  const double *u,
                  const double *G);
}

namespace cj { namespace elastic {

class ReducedToRS {
public :
    typedef zjucad::matrix::matrix<size_t> matrixi_t;
    typedef zjucad::matrix::matrix<double> matrixd_t;
    ReducedToRS(const matrixi_t &tets,
                const matrixd_t &nods,
                const matrixd_t &G,
                const Eigen::MatrixXd &W) {
        Eigen::SparseMatrix<double> Q;
        std::vector<Eigen::Triplet<double>> trips;
#pragma omp parallel for
        for (size_t i = 0; i < tets.size(2); ++i) {
            matrixd_t H(9, 12);
            disp2rs_jac_(&H[0], nullptr, &G(0, i));
            for (size_t p = 0; p < 9; ++p) {
                for (size_t q = 0; q < 12; ++q) {
                    if ( H(p, q) != 0.0 ) {
                        size_t I = 9 * i + p;
                        size_t J = 3 * tets(q / 3, i) + q % 3;
                    #pragma omp critical
                        trips.push_back(Eigen::Triplet<double>(I, J, H(p, q)));
                    }
                }
            }
        }
        Q.resize(9 * tets.size(2), nods.size());
        Q.reserve(trips.size());
        Q.setFromTriplets(trips.begin(), trips.end());
        QW_.resize(Q.rows(), W.cols());
        QW_ = Q * W;
    }
    void operator ()(const Eigen::VectorXd &z,
                     matrixd_t &rs) {
        Eigen::Map<Eigen::VectorXd>(&rs[0], rs.size())
                = QW_ * z;
    }
private :
    Eigen::SparseMatrix<double> QW_;
};
}}

#endif
