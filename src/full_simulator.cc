#include "full_simulator.h"

#include <iostream>
#include <Eigen/UmfPackSupport>
#include <Eigen/Dense>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
//#include <hjlib/util/hrclock.h>
#include "elastic_energy.h"
#include "position_cons.h"
#include "mass_matrix.h"


using namespace zjucad::matrix;

namespace cj { namespace elastic {

StVKSimulator::StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                             const zjucad::matrix::matrix<double> &nods,
                             boost::property_tree::ptree &pt)
    : tets_(tets), nods_(nods), pt_(pt) {
    h_    = pt_.get<double>("elastic.time_step");
    alpha_= pt_.get<double>("elastic.alpha");
    beta_ = pt_.get<double>("elastic.beta");
    disp_ = zeros<double>(3, nods_.size(2));
    fext_ = zeros<double>(3, nods_.size(2));

    // set mass matrix, here we want an unlumped matrix
    double dens = pt_.get<double>("elastic.density");
    MassMatrix mass_calculator(tets_, nods_, dens);
    mass_calculator.Compute(M_, false);

    // comupte lame first & second parameters
    // according to Young's modulus and Possion ratio
    double E = pt_.get<double>("elastic.YoungModulus", 2e6);
    double v = pt_.get<double>("elastic.PoissonRatio", 0.45);
    double lambda = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
    double miu = E / (2.0 * (1.0 + v));

    std::string model = pt_.get<std::string>("elastic.consitutive_model");

    try {
        pe_.reset(BuildElasticEnergy(tets_, nods_, lambda, miu, model));
        x_.resize(pe_->Nx());
        x_.setZero();
    } catch (...) {
        std::cerr << "[INFO] error: in building elastic energy, and the prgram will exit.\n";
        exit(0);
    }
}

void StVKSimulator::SetFixedPoints(const std::vector<size_t> &idx,
                                   const matrix<double> &uc) {
    std::cerr << "[INFO] the number of fixed points is: " << idx.size() << "\n";
    double pos_penalty = pt_.get<double>("elastic.pos_penalty");
    pc_.reset(new PositionCons(idx, uc));
    x_.resize(pe_->Nx() + pc_->Nf());
    x_.setZero();
}

void StVKSimulator::ClearFixedPoints() {
    if ( pc_.get() )
        pc_.reset();
}

void StVKSimulator::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) = itr_matrix<const double *>(3, 1, force);
}

void StVKSimulator::ClearExternalForce() {
    fext_ = zeros<double>(3, nods_.size(2));
}

int StVKSimulator::Forward() {
    //hj::util::high_resolution_clock hrc;
    //double start = hrc.ms();

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd             b;
    int flagA = AssembleLHS(A);
    int flagB = AssembleRHS(b);
    assert(!flagA && !flagB);

    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if ( solver.info() != Eigen::Success ) {
        std::cerr << "[INFO] decomposition failed.\n";
        return __LINE__;
    }
    x_ = solver.solve(b);
    if ( solver.info() != Eigen::Success ) {
        std::cerr << "[INFO] solve failed.\n";
        return __LINE__;
    }  
    Eigen::Map<Eigen::VectorXd>(&disp_[0], disp_.size()) += h_ * x_.head(nods_.size());
    //std::cerr << "[INFO] time cost: " << hrc.ms() - start << "\n\n";
    return 0;
}

matrix<double>& StVKSimulator::disp() {
    return disp_;
}

int StVKSimulator::AssembleLHS(Eigen::SparseMatrix<double> &A) {
    if ( !pe_.get() ) {
        std::cerr << "[INFO] no energy to drive the deformation!\n";
        return __LINE__;
    }

    Eigen::SparseMatrix<double> K, C, L;
    size_t dim1, dim2;

    dim1 = pe_->Nx();
    pe_->Hes(&disp_[0], &K);
    L = (1 + h_ * alpha_) * M_ + h_ * (h_ + beta_) * K;
    L.makeCompressed();

    if ( pc_.get() ) {
        dim2 = pc_->Nf();
        pc_->Jac(&disp_[0], &C);
        C.makeCompressed();
    } else {
        dim2 = 0;
    }

    // [ L  CT ]|x  | = [Mv + h(fext - f) ]
    // [ C   0 ]|lam|   [  C(uc - u) / h  ]
    std::vector<Eigen::Triplet<double>> trips;
    for (size_t j = 0; j < dim1; ++j) {
        for (size_t cnt = L.outerIndexPtr()[j]; cnt < L.outerIndexPtr()[j + 1]; ++cnt) {
            trips.push_back(Eigen::Triplet<double>(L.innerIndexPtr()[cnt], j, L.valuePtr()[cnt]));
        }
    }
    if ( pc_.get() ) {
        for (size_t j = 0; j < dim1; ++j) {
            for (size_t cnt = C.outerIndexPtr()[j]; cnt < C.outerIndexPtr()[j + 1]; ++cnt) {
                trips.push_back(Eigen::Triplet<double>(dim1 + C.innerIndexPtr()[cnt], j, C.valuePtr()[cnt]));
                trips.push_back(Eigen::Triplet<double>(j, C.innerIndexPtr()[cnt] + dim1, C.valuePtr()[cnt]));
            }
        }
    }
    A.resize(dim1 + dim2, dim1 + dim2);
    A.reserve(trips.size());
    A.setFromTriplets(trips.begin(), trips.end());
    A.makeCompressed();
    return 0;
}

int StVKSimulator::AssembleRHS(Eigen::VectorXd &rhs) {
    if ( !pe_.get() ) {
        std::cerr << "[INFO] no energy to drive the deformation!\n";
        return __LINE__;
    }
    size_t dim1, dim2;

    dim1 = pe_->Nx();
    Eigen::VectorXd f(dim1);
    f.setZero();
    pe_->Gra(&disp_[0], f.data());

    Eigen::VectorXd v;
    if ( pc_.get() ) {
        dim2 = pc_->Nf();
        v.resize(dim2);
        v.setZero();
        pc_->Val(&disp_[0], v.data());
    } else {
        dim2 = 0;
    }

    rhs.resize(dim1 + dim2);
    rhs.head(dim1) = M_ * x_.head(dim1)
            + h_ * (Eigen::Map<Eigen::VectorXd>(&fext_[0], fext_.size()) - f);
    if ( pc_.get() )
        rhs.tail(dim2) = -v / h_;
    return 0;
}

}}
