#include "stvk_simulator.h"

#include <iostream>
#include <Eigen/UmfPackSupport>
#include <Eigen/Dense>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include "stvk_energy.h"
#include "position_cons.h"
#include "mass_matrix.h"

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

    // set mass matrix, here we want an unlumped matrix
    double dens = pt_.get<double>("stvk.density");
    MassMatrix mass_calculator(tets_, nods_, dens);
    mass_calculator.Compute(M_, false);

    // comupte lame first & second parameters
    // according to Young's modulus and Possion ratio
    double E = pt_.get<double>("stvk.YoungModulus", 2e6);
    double v = pt_.get<double>("stvk.PoissonRatio", 0.45);
    double lambda = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
    double miu = E / (2.0 * (1.0 + v));

    pe_.reset(new StVKEnergy(tets_, nods_, lambda, miu));
    x_.resize(pe_->Nx());
    x_.setZero();
}

void StVKSimulator::SetFixedPoints(const vector<size_t> &idx,
                                   const matrix<double> &uc) {
    cerr << "[INFO] the number of fixed points is: " << idx.size() << endl;
    double pos_penalty = pt_.get<double>("stvk.pos_penalty");
    pc_.reset(new PositionCons(idx, uc));
    x_.resize(pe_->Nx() + pc_->Nf());
    x_.setZero();
}

void StVKSimulator::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) = itr_matrix<const double *>(3, 1, force);
}

void StVKSimulator::ClearExternalForce() {
    fext_ = zeros<double>(3, nods_.size(2));
}

int StVKSimulator::Forward() {
    SparseMatrix<double> A;
    VectorXd             b;

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
    L = (1 + h_ * alpha_) * M_ + h_ * (h_ + beta_) * K;
    L.makeCompressed();

    pc_->Jac(&disp_[0], &C);
    C.makeCompressed();

    size_t dim1 = pe_->Nx();
    size_t dim2 = pc_->Nf();

    // [ L  CT ]|x  | = [Mv + h(fext - f) ]
    // [ C   0 ]|lam|   [  C(uc - u) / h  ]
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
    return 0;
}

int StVKSimulator::AssembleRHS(VectorXd &rhs) {
    if ( !pe_.get() || !pc_.get() ) {
        cerr << "[INFO] null pointer.\n";
        return __LINE__;
    }
    size_t dim1 = pe_->Nx();
    size_t dim2 = pc_->Nf();

    VectorXd f(dim1); f.setZero();
    VectorXd v(dim2); v.setZero();
    pe_->Gra(&disp_[0], f.data());
    pc_->Val(&disp_[0], v.data());

    rhs.resize(dim1 + dim2);
    rhs.head(dim1) = M_ * x_.head(dim1)
            + h_ * (Map<VectorXd>(&fext_[0], fext_.size()) - f);
    rhs.tail(dim2) = -v / h_;
    return 0;
}

}}
