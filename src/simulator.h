#ifndef __ELASTIC_FULL_SIMULATOR_H__
#define __ELASTIC_FULL_SIMULATOR_H__

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <unordered_set>

namespace cj { namespace elastic {

class Energy;
class Constraint;
class ModalAnalyzer;
class ReducedToRS;

// displacement based stvk model

class StVKSimulator {
public :
    StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                  const zjucad::matrix::matrix<double> &nods,
                  boost::property_tree::ptree &pt);
    void SetFixedPoints(const std::vector<size_t> &idx,
                        const zjucad::matrix::matrix<double> &uc);
    void ClearFixedPoints();
    void SetExternalForce(const size_t idx, const double *force);
    void ClearExternalForce();
    int Advance();
    zjucad::matrix::matrix<double>& disp();

private :
    int AssembleLHS(Eigen::SparseMatrix<double> &A);
    int AssembleRHS(Eigen::VectorXd &rhs);

    ///< geometry
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;

    ///< energy and constraint
    std::shared_ptr<Energy>     pe_;
    std::shared_ptr<Constraint> pc_;

    ///< physics
    boost::property_tree::ptree &pt_;
    double h_, alpha_, beta_;
    Eigen::VectorXd x_;         ///< store velocity and lagragian multipliers
    zjucad::matrix::matrix<double> disp_;
    zjucad::matrix::matrix<double> fext_;
    Eigen::SparseMatrix<double> M_;
};

class ReducedSolver {
public :
    ReducedSolver(const zjucad::matrix::matrix<size_t> &tets,
                  const zjucad::matrix::matrix<double> &nods,
                  boost::property_tree::ptree &pt);
    int Init();
    int AddElasticEnergy(const double w);
    int SetPinnedVertices(const std::vector<size_t> &idx,
                          const zjucad::matrix::matrix<double> &uc);
    int SetExternalForce(const size_t idx, const double *force);
    int ClearExternalForce();

    /// @brief prepare for modal basis and discrete gradient operator
    /// @return 0 if succeed
    int Prepare();
    int Advance();
    zjucad::matrix::matrix<double>& get_disp(); // directly or warp
    size_t GetSubspaceDim() const { return nbrBasis_; }
public :
    int BuildModalBasis(const std::unordered_set<size_t> &fix);
    void VisualizeVibrationModes();
    // RS warping
    int ComputeRSCoords(const zjucad::matrix::matrix<double> &u);
    int ComputeRSCoords(const Eigen::VectorXd &z);
    int RSWarping();

public :
    ///< geometry
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;

    ///< for warping
    zjucad::matrix::matrix<double> vols_;   // #tets x 1
    zjucad::matrix::matrix<double> tetRS_;  // 9 x #tets
    zjucad::matrix::matrix<double> G_;      // 9 x #tets
    std::shared_ptr<Energy>     pe_warp_;
    std::shared_ptr<Constraint> pc_warp_;
    Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver_;
    Eigen::SparseMatrix<double> LHS_;
    std::shared_ptr<ReducedToRS> reducedtoRS;

    ///< energies and constraints
    std::shared_ptr<Energy>     pe_;
    std::shared_ptr<Constraint> pc_;

    ///< physics configuration
    boost::property_tree::ptree &pt_;
    double h_, alpha_, beta_;
    zjucad::matrix::matrix<double> disp_;
    zjucad::matrix::matrix<double> fext_;
    Eigen::SparseMatrix<double> M_;

    ///< Reduced base
    std::unordered_set<size_t> fixed_;
    std::shared_ptr<ModalAnalyzer> basis_builder_;
    size_t nbrBasis_ = 0;
    Eigen::VectorXd z_, dotz_;
    Eigen::MatrixXd U_;
    Eigen::VectorXd lambda_;

};

}}
#endif
