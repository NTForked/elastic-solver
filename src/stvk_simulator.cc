#include "stvk_simulator.h"

#include <iostream>
#include <Eigen/UmfPackSupport>
#include <zjucad/matrix/itr_matrix.h>
#include "stvk_energy.h"
#include "position_cons.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

StVKSimulator::StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                             const zjucad::matrix::matrix<double> &nods,
                             boost::property_tree::ptree &pt)
    : tets_(tets), nods_(nods) {
    h_ =  pt.get<double>("stvk.time_step");
    alpha_ = pt.get<double>("stvk.alpha");
    beta_ = pt.get<double>("stvk.beta");

    disp_ = zeros<double>(3, nods_.size(2));
    fext_ = zeros<double>(3, nods_.size(2));
    x_.resize();

    double lambda = pt.get<double>("stvk.lame_lambda");
    double miu = pt.get<double>("stvk.lame_niu");
    pe_.reset(new StVKEnergy(tets_, nods_, lambda, miu));

    // set mass matrix
    double dens = pt.get<double>("stvk.density");

}

int StVKSimulator::SetFixedPoint(std::vector<size_t> &idx) {
    pc_.reset(new PositionCons());
}

int StVKSimulator::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) += itr_matrix<const double *>(3, 1, force);
}

int StVKSimulator::ClearExternalForce() {
    fext_ = zeros<double>(3, nods_.size(2));
    return 0;
}

int StVKSimulator::Forward() {
    SparseMatrix<double> A;
    VectorXd b;
    AssembleA(A);
    AssembleRhs(b);

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
    disp_ = disp_ + h_ * x_.head();
    return 0;
}

int StVKSimulator::AssembleLHS(Eigen::SparseMatrix<double> &A) {
    SparseMatrix<double> K, C, L;
    pe_->Hes(disp.data(), &K);
    pc_->Jac(disp.data(), &C);
    L = (1 + h_ * alpha_) * M_  + (h_ * h_ + h_ * beta_) * K;
    C *= h_;

    size_t dim1 = pe_->Nx();
    size_t dim2 = pc_->Nf();

    // push back L

    // push back hC and hCT


    A.resize(dim1 + dim2, dim1 + dim2);
    A.reserve(trips.size());
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    return 0;
}

int StVKSimulator::AssembleRHS(VectorXd &rhs) {
    VectorXd g();
    g.setZero();
    pe_->Gra(disp_.data(), g.data());
    rhs.head() = M * x.head() + h_ * fext_ - h_ * g;
    rhs.tail() = -pc_->Val();
    return 0;
}

