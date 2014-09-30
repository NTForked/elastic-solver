#ifndef __CJ_STVK_MODEL_H__
#define __CJ_STVK_MODEL_H__

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

class Energy;
class Constraint;

// displacement based stvk model

class StVKSimulator {
public :
    StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                  const zjucad::matrix::matrix<double> &nods,
                  boost::property_tree::ptree &pt);
    void SetFixedPoints(const std::vector<size_t> &idx,
                        const zjucad::matrix::matrix<double> &uc);
    void SetExternalForce(const size_t idx, const double *force);
    void ClearExternalForce();
    int Forward();
    zjucad::matrix::matrix<double>& disp();

private :
    int AssembleLHS(Eigen::SparseMatrix<double> &A);
    int AssembleRHS(Eigen::VectorXd &rhs);

    // geometry model
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;

    // energy and constraint
    std::shared_ptr<Energy>     pe_;
    std::shared_ptr<Constraint> pc_;

    // physics model
    boost::property_tree::ptree pt_;
    double h_, alpha_, beta_;
    Eigen::VectorXd x_;         // store velocity and lagragian multipliers
    zjucad::matrix::matrix<double> disp_;
    zjucad::matrix::matrix<double> fext_;
    Eigen::SparseMatrix<double> M_;
};

}}
#endif
