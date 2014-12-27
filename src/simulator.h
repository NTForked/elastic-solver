#ifndef __ELASTIC_FULL_SIMULATOR_H__
#define __ELASTIC_FULL_SIMULATOR_H__

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unordered_set>

namespace cj { namespace elastic {

class Energy;
class Constraint;
class ModalAnalyzer;

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
    int BuildModalBasis(const std::unordered_set<size_t> &fix);
    void VisualizeVibrationModes();
public :
    ///< geometry
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;

    ///< energies and constraints
    std::shared_ptr<Energy> pe_;
    std::shared_ptr<Constraint> pc_;

    ///< physics configuration
    boost::property_tree::ptree &pt_;
    double h_, alpha_, beta_;
    zjucad::matrix::matrix<double> disp_;
    zjucad::matrix::matrix<double> fext_;
    Eigen::SparseMatrix<double> M_;

    ///< Reduced base
    std::shared_ptr<ModalAnalyzer> basis_builder_;
    Eigen::MatrixXd U_;
    Eigen::VectorXd lambda_;

};

}}
#endif
