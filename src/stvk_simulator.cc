#include "stvk_simulator.h"

#include <iostream>
#include <Eigen/UmfPackSupport>
#include <Eigen/Dense>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include "stvk_energy.h"
#include "position_cons.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace cj { namespace elastic {

StVKSimulator::StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                             const zjucad::matrix::matrix<double> &nods,
                             boost::property_tree::ptree &pt)
    : tets_(tets), nods_(nods), pt_(pt) {
    h_    = pt_.get<double>("stvk.time_step");
    alpha_= pt_.get<double>("stvk.alpha");
    beta_ = pt_.get<double>("stvk.beta");
    disp_ = zeros<double>(3, nods_.size(2));
    fext_ = zeros<double>(3, nods_.size(2));

    // set mass matrix
    double dens = pt_.get<double>("stvk.density");
    vector<Triplet<double>> mass_trip;
    for (size_t i = 0; i < nods_.size(); ++i) {
        mass_trip.push_back(Triplet<double>(i, i, dens));
    }
    M_.resize(nods_.size(), nods_.size());
    M_.reserve(mass_trip.size());
    M_.setFromTriplets(mass_trip.begin(), mass_trip.end());
    M_.makeCompressed();

    // comupte lame first & second parameters
    double E = pt_.get<double>("stvk.Young_modulus");
    double v = pt_.get<double>("stvk.Poisson_ratio");
    double lambda = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
    double miu = E / 2.0 * (1.0 + v);

    x_.resize(nods_.size());
    x_.setZero();
    pe_.reset(new StVKEnergy(tets_, nods_, lambda, miu));
}

void StVKSimulator::SetFixedPoints(const vector<size_t> &idx,
                                   const matrix<double> &uc) {
    cerr << "[INFO] the number of fixed points is: " << idx.size() << endl;
    double position_penalty = pt_.get<double>("stvk.pos_penalty");
    pc_.reset(new PositionCons(idx, uc));
    x_.resize(pe_->Nx() + pc_->Nf());
    x_.setZero();
}

void StVKSimulator::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) += itr_matrix<const double *>(3, 1, force);
}

void StVKSimulator::ClearExternalForce() {
    fext_ = zeros<double>(3, nods_.size(2));
}

int StVKSimulator::Forward() { // Ax = b
    SparseMatrix<double> A;
    VectorXd b;
    AssembleLHS(A);
    AssembleRHS(b);

    UmfPackLU<SparseMatrix<double>> solver;
    solver.compute(A);
    if ( solver.info() != Success ) {
        cerr << "[INFO] decomposition failed.\n";
        return __LINE__;
    }
    x_ = solver.solve(b);
    if ( solver.info() != Success ) {
        cerr << "[INFO] solve failed.\n";
        return __LINE__;
    }
    Map<VectorXd>(&disp_[0], disp_.size()) += h_ * x_.head(nods_.size());
    return 0;
}

matrix<double>& StVKSimulator::disp() {
    return disp_;
}

int StVKSimulator::AssembleLHS(Eigen::SparseMatrix<double> &A) {
    SparseMatrix<double> K, C, L;
    vector<Triplet<double>> trips;
    if ( !pe_.get() || !pc_.get() ) {
        cerr << "[INFO] null pointer.\n";
        return __LINE__;
    }

    pe_->Hes(&disp_[0], &K);
    pc_->Jac(&disp_[0], &C);

    L = (1 + h_ * alpha_) * M_  + h_ * (h_ + beta_) * K;
    C *= h_;
    L.makeCompressed();
    C.makeCompressed();

//    MatrixXd K_(K.rows(), K.cols());
//    for (size_t col = 0; col < L.cols(); ++col) {
//        for (size_t cnt = K.outerIndexPtr()[col]; cnt < K.outerIndexPtr()[col + 1]; ++cnt)
//            K_(K.innerIndexPtr()[cnt], col) = K.valuePtr()[cnt];
//    }
//    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(K_);
//    for (size_t i = 0; i < 10; ++i)
//        cout << solver.eigenvalues()[i] << endl;
//    exit(0);

//    MatrixXd L_(L.rows(), L.cols());
//    for (size_t col = 0; col < L.cols(); ++col) {
//        for (size_t cnt = L.outerIndexPtr()[col]; cnt < L.outerIndexPtr()[col + 1]; ++cnt)
//            L_(L.innerIndexPtr()[cnt], col) = L.valuePtr()[cnt];
//    }
//    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(L_);
//    for (size_t i = 0; i < 10; ++i)
//        cout << solver.eigenvalues()[i] << endl;
//    exit(0);

    size_t dim1 = pe_->Nx();
    size_t dim2 = pc_->Nf();

    // | L  hCT |
    // | hC   0 |
    for (size_t j = 0; j < dim1; ++j) {
        for (size_t cnt = L.outerIndexPtr()[j]; cnt < L.outerIndexPtr()[j + 1]; ++cnt) {
            trips.push_back(Triplet<double>(L.innerIndexPtr()[cnt], j, L.valuePtr()[cnt]));
        }
    }
    for (size_t j = 0; j < dim1; ++j) {
        for (size_t cnt = C.outerIndexPtr()[j]; cnt < C.outerIndexPtr()[j + 1]; ++cnt) {
            trips.push_back(Triplet<double>(dim1 + C.innerIndexPtr()[cnt], j, C.valuePtr()[cnt]));
            trips.push_back(Triplet<double>(j, C.innerIndexPtr()[cnt] + dim1, C.valuePtr()[cnt]));
        }
    }
    A.resize(dim1 + dim2, dim1 + dim2);
    A.reserve(trips.size());
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();

//    MatrixXd A_(A.rows(), A.cols());
//    A.setZero();
//    for (size_t j = 0; j < A.cols(); ++j) {
//        for (size_t cnt = A.outerIndexPtr()[j]; cnt < A.outerIndexPtr()[j + 1]; ++cnt)
//            A_(A.innerIndexPtr()[cnt], j) = A.valuePtr()[cnt];
//    }
//    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(A_);
//    cout << solver.eigenvalues()[0] << endl;
//    exit(0);

    return 0;
}

int StVKSimulator::AssembleRHS(VectorXd &rhs) {
    if ( !pe_.get() || !pc_.get() ) {
        cerr << "[INFO] null pointer.\n";
        return __LINE__;
    }
    size_t dim1 = pe_->Nx();
    size_t dim2 = pc_->Nf();
    rhs.resize(dim1 + dim2);
    VectorXd g(dim1), v(dim2);
    g.setZero(); v.setZero();
    pe_->Gra(&disp_[0], g.data());
    g = -g;
    rhs.head(dim1) = M_ * x_.head(dim1)
            + h_ * ( Map<VectorXd>(&fext_[0], fext_.size()) - g);
    pc_->Val(&disp_[0], v.data());
    rhs.tail(dim2) = -v;
    return 0;
}

}}
