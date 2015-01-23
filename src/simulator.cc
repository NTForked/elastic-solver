#include "simulator.h"

#include <iostream>
#include <Eigen/Dense>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/util/hrclock.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "elastic_energy.h"
#include "warping.h"
#include "position_cons.h"
#include "mass_matrix.h"
#include "modal_analyzer.h"
#include "util.h"
#include "reduced2rs.h"

using namespace std;
using namespace Eigen;
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
    mass_calculator.Compute(M_, true);

    // comupte lame first & second parameters
    // according to Young's modulus and Possion ratio
    double E = pt_.get<double>("elastic.Young_modulus");
    double v = pt_.get<double>("elastic.Poisson_ratio");
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

int StVKSimulator::Advance() {
    hj::util::high_resolution_clock hrc;
    double start = hrc.ms();

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
    std::cerr << "[INFO] time cost: " << hrc.ms() - start << "\n\n";
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

//------------------------LINEAR REDUCED PART-----------------------------------

LinearReducedSolver::LinearReducedSolver(const matrix<size_t> &tets,
                                         const matrix<double> &nods,
                                         boost::property_tree::ptree &pt)
    : tets_(tets), nods_(nods), pt_(pt) { }

int LinearReducedSolver::Init() {
    h_ = pt_.get<double>("elastic.time_step");
    alpha_ = pt_.get<double>("elastic.alpha");
    beta_ = pt_.get<double>("elastic.beta");
    disp_ = zeros<double>(3, nods_.size(2));
    fext_ = zeros<double>(3, nods_.size(2));

    // set mass matrix, here we want an lumped matrix;
    double dens = pt_.get<double>("elastic.density");
    MassMatrix calculator(tets_, nods_, dens);
    bool lumped = true;
    calculator.Compute(M_, lumped);

    return 0;
}

int LinearReducedSolver::AddElasticEnergy(const double w) {
    const double E = pt_.get<double>("elastic.Young_modulus");
    const double v = pt_.get<double>("elastic.Poisson_ratio");
    const double lambda = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
    const double miu = E / (2.0 * (1.0 + v));

    std::string material = pt_.get<std::string>("elastic.consitutive_model");
    pe_.reset(BuildElasticEnergy(tets_, nods_, lambda, miu, material));
    if ( pe_.get() == nullptr )
        return __LINE__;
    return 0;
}

int LinearReducedSolver::SetPinnedVertices(const vector<size_t> &idx,
                                           const matrix<double> &uc) {
    for (size_t i = 0; i < idx.size(); ++i) {
        fixed_.insert(3 * idx[i] + 0);
        fixed_.insert(3 * idx[i] + 1);
        fixed_.insert(3 * idx[i] + 2);
    }
    return 0;
}

int LinearReducedSolver::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) = itr_matrix<const double *>(3, 1, force);
    return 0;
}

int LinearReducedSolver::ClearExternalForce() {
    fext_ = zeros<double>(3, nods_.size(2));
    return 0;
}

int LinearReducedSolver::Prepare() {
    // prepare for solving in reduced space
    BuildModalBasis(fixed_);
    z_.resize(nbrBasis_);
    z_.setZero();
    dotz_.resize(nbrBasis_);
    dotz_.setZero();

    // prepare for RS coordinates warping
    vols_ = zeros<double>(tets_.size(2), 1);
    tetRS_.resize(9, tets_.size(2));
    G_.resize(9, tets_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> V = nods_(colon(), tets_(colon(), i));
        matrix<double> G = V(colon(), colon(1, 3)) - V(colon(), 0) * ones<double>(1, 3);
        matrix<double> Gbkp = G;
        vols_[i] = fabs(det(Gbkp)) / 6.0;
        if ( inv(G) ) {
            cerr << "[INFO] degenerated tet\n";
        }
        G_(colon(), i) = itr_matrix<const double*>(9, 1, &G[0]);
    }
    pe_warp_.reset(new WarpingEnergy(tets_, nods_, G_, tetRS_, vols_));
    pc_warp_.reset(new PositionCons(vector<size_t>(1, 0), zeros<double>(3, nods_.size(2))));
    SparseMatrix<double> H, J;
    pe_warp_->Hes(nullptr, &H);
    pc_warp_->Jac(nullptr, &J);
    vector<Triplet<double>> trips;
    for (size_t j = 0; j < H.cols(); ++j)
        for (SparseMatrix<double>::InnerIterator it(H, j); it; ++it)
            trips.push_back(Triplet<double>(it.row(), it.col(), it.value()));
    for (size_t j = 0; j < J.cols(); ++j) {
        for (SparseMatrix<double>::InnerIterator it(J, j); it; ++it) {
            trips.push_back(Triplet<double>(pe_warp_->Nx() + it.row(), it.col(), it.value()));
            trips.push_back(Triplet<double>(it.col(), pe_warp_->Nx() + it.row(), it.value()));
        }
    }
    const size_t DIM = pe_warp_->Nx() + pc_warp_->Nf();
    LHS_.resize(DIM, DIM);
    LHS_.reserve(trips.size());
    LHS_.setFromTriplets(trips.begin(), trips.end());
    solver_.compute(LHS_);
    assert(solver_.info() == Success);
    // build the linear map between reduced coordinates
    // and RS coordinates
    reducedtoRS.reset(new ReducedToRS(tets_, nods_, G_, U_));

    printf("[INFO] preparation done\n");
    return 0;
}

int LinearReducedSolver::Advance() {
    // just for linear elasticity, reduced stiffness
    // matrix is constant and simply diagonal
    // no need to evaluate the cubic reduced force
    // and quadratic reduced stiffness matrix as
    // stvk material
    hj::util::high_resolution_clock clk;
    double timestart = clk.ms();
    Map<VectorXd> f(&fext_[0], fext_.size());

    ///> integrate every component separately, implicit euler
#pragma omp parallel for
    for (size_t i = 0; i < z_.size(); ++i) {
        dotz_[i] = (dotz_[i] + h_ * (f.transpose() * U_.col(i) - lambda_[i] * z_[i])) /
                (1.0 + h_ * (beta_ * lambda_[i] + alpha_ + lambda_[i] * h_));
        z_[i] += h_ * dotz_[i];
    }

    printf("\ttime cost: %lf\n", clk.ms() - timestart);
    return 0;
}

matrix<double>& LinearReducedSolver::get_disp() {
#define CONVENTIONAL_LINEAR_ELASTICITY 0
#if CONVENTIONAL_LINEAR_ELASTICITY
    ///> conventional linear elasticity
    Map<VectorXd>(&disp_[0], disp_.size()) = U_ * z_;
#else
    ///> RS-coordinates warping
    ComputeRSCoords(z_);
    RSWarping();
#endif
    return disp_;
#undef CONVENTIONAL_LINEAR_ELASTICITY
}

int LinearReducedSolver::BuildModalBasis(const std::unordered_set<size_t> &fix) {
    Eigen::SparseMatrix<double> K;
    matrix<double> X = zeros<double>(3, nods_.size(2));
    if ( pe_.get() )
        pe_->Hes(&X[0], &K);
    else
        return __LINE__;

    bool flagK = isSymmetric<Eigen::SparseMatrix<double>>(K);
    bool flagM = isSymmetric<Eigen::SparseMatrix<double>>(M_);
    if ( !(flagK && flagM) ) {
        std::cerr << "[INFO] K or M is not symmetric\n";
        return __LINE__;
    }

    nbrBasis_ = pt_.get<size_t>("elastic.basis_number");
    cout << "[INFO] number of modal basis is " << nbrBasis_ << "\n";
    basis_builder_.reset(new ModalAnalyzer(K, M_, nbrBasis_, fix));
    basis_builder_->Compute();
    U_ = basis_builder_->get_modes();
    lambda_ = basis_builder_->get_freqs();

    return 0;
}

int LinearReducedSolver::ComputeRSCoords(const matrix<double> &u) {
#pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
        matrix<double> U = u(colon(), tets_(colon(1, 3), i))
                - u(colon(), tets_(0, i)) * ones<double>(1, 3);
        matrix<double> DF = U * itr_matrix<const double*>(3, 3, &G_(0, i));
        matrix<double> skew = 0.5 * (DF - trans(DF));
        matrix<double> symm = 0.5 * (DF + trans(DF));
        tetRS_(0, i) = skew(1, 0);
        tetRS_(1, i) = skew(2, 0);
        tetRS_(2, i) = skew(2, 1);
        tetRS_(3, i) = symm(0, 0);
        tetRS_(4, i) = symm(1, 0);
        tetRS_(5, i) = symm(2, 0);
        tetRS_(6, i) = symm(1, 1);
        tetRS_(7, i) = symm(2, 1);
        tetRS_(8, i) = symm(2, 2);
    }
    return 0;
}

int LinearReducedSolver::ComputeRSCoords(const VectorXd &z) {
    (*reducedtoRS)(z, tetRS_);
    return 0;
}

int LinearReducedSolver::RSWarping() {
    static const size_t EDIM = pe_warp_->Nx();
    static const size_t CDIM = pc_warp_->Nf();
    VectorXd b(EDIM + CDIM);
    b.setZero();
    if ( pe_warp_.get() )
        pe_warp_->Gra(&disp_[0], &b[0]);
    if ( pc_warp_.get() )
        pc_warp_->Val(&disp_[0], &b[EDIM]);
    b = -b;
    VectorXd solution = solver_.solve(b);
    assert(solver_.info() == Success);
    disp_ += itr_matrix<const double*>(3, disp_.size(2), solution.data());
    return 0;
}

//--------------------------STVK REDUCED PART-----------------------------------

StvkReducedSolver::StvkReducedSolver(const matrixi_t &tets,
                                     const matrixd_t &nods,
                                     boost::property_tree::ptree &pt)
    : tets_(tets), nods_(nods), pt_(pt) { }

int StvkReducedSolver::Init() {
    h_ = pt_.get<double>("elastic.time_step");
    alpha_ = pt_.get<double>("elastic.alpha");
    beta_ = pt_.get<double>("elastic.beta");
    disp_ = zeros<double>(3, nods_.size(2));
    fext_ = zeros<double>(3, nods_.size(2));
    grav_ = zeros<double>(3, nods_.size(2));

    double rho = pt_.get<double>("elastic.density");
    MassMatrix calculator(tets_, nods_, rho);
    bool lumped = true;
    calculator.Compute(M_, lumped);
    return 0;
}

int StvkReducedSolver::AddElasticEnergy(const double w) {
    const double E = pt_.get<double>("elastic.Young_modulus");
    const double v = pt_.get<double>("elastic.Poisson_ratio");
    const double lambda = E * v / ((1.0 + v) * (1.0 - 2.0 * v));
    const double miu = E / (2.0 * (1.0 + v));
    pe_.reset(BuildElasticEnergy(tets_, nods_, lambda, miu, "stvk"));
    if ( pe_.get() == nullptr )
        return __LINE__;
    return 0;
}

void StvkReducedSolver::SetPinnedVerts(const vector<size_t> &idx, const matrixd_t &uc) {
    for (size_t i = 0; i < idx.size(); ++i) {
        fixed_.insert(3 * idx[i] + 0);
        fixed_.insert(3 * idx[i] + 1);
        fixed_.insert(3 * idx[i] + 2);
    }
}

/// rebuild modal basis???
void StvkReducedSolver::FreePinnedVerts() {
    fixed_.clear();
}

void StvkReducedSolver::SetGravity(const double w) {
#pragma omp parallel for
    for (size_t i = 0; i < grav_.size(2); ++i)
        grav_(1, i) = -9.81 * w * M_.coeff(3 * i, 3 * i);
    fext_ += grav_;
}

void StvkReducedSolver::ClearGravity() {
    fext_ -= grav_;
}

void StvkReducedSolver::SetExternalForce(const size_t idx, const double *force) {
    fext_(colon(), idx) = itr_matrix<const double*>(3, 1, force);
}

void StvkReducedSolver::ClearExternalForce() {
    fext_ = zeros<double>(fext_.size(1), fext_.size(2));
}

int StvkReducedSolver::Prepare() {
    BuildModalBasis(fixed_);
    z_.resize(nbrBasis_, 1);
    dzdt_.setZero(nbrBasis_, 1);
    return 0;
}

int StvkReducedSolver::Advance() {
    // subspace integration
    // reduced equations of motion:
    // \ddot q + U^TD(Uq, U\dot q) + U^TR(Uq) = U^Tf_{ext}
    // apply Rayleigh damping
    // \ddot q + U^T(\alpha M + \beta K(Uq)U\dot q + U^TR(Uq) = U^Tf_{ext}

}

StvkReducedSolver::matrixd_t& StvkReducedSolver::get_disp() {
    Map<VectorXd>(&disp_[0], disp_.size()) = U_ * z_;
    return disp_;
}

int StvkReducedSolver::SolveModalDeriv() {

}

int StvkReducedSolver::BuildModalBasis(const unordered_set<size_t> &fix) {
    Eigen::SparseMatrix<double> K;
    matrixd_t u0 = zeros<double>(3, nods_.size(2));
    if ( pe_.get() )
        pe_->Hes(&u0[0], &K);
    else
        return __LINE__;

    nbrBasis_ = pt_.get<size_t>("elastic.basis_number");
    basis_builder_.reset(new ModalAnalyzer(K, M_, nbrBasis_, fix));
    basis_builder_->Compute();
    U_ = basis_builder_->get_modes();
    lambda_ = basis_builder_->get_freqs();

    // extend basis through modal derivative
    //.......

    return 0;
}

}}
