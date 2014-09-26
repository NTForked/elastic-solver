#ifndef __CJ_STVK_MODEL_H__
#define __CJ_STVK_MODEL_H__

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

class StVKEnergy;
class PositionCons;

// displacement based stvk model

class StVKSimulator  {
public :
    StVKSimulator(const zjucad::matrix::matrix<size_t> &tets,
                                 const zjucad::matrix::matrix<double> &nods,
                                 boost::property_tree::ptree &pt);

    int SetFixedPoint(std::vector<size_t> &idx);
    int SetExternalForce(const size_t idx, const double *force);
    int ClearExternalForce();
    int Forward();
    zjucad::matrix::matrix<double>& disp();

private :
    int AssembleLHS(Eigen::SparseMatrix<double> &A);
    int AssembleRHS(Eigen::VectorXd &rhs);

    // geometry model
    const zjucad::matrix::matrix<size_t> &tets_;
    const zjucad::matrix::matrix<double> &nods_;

    // energy and constraint
    std::unique_ptr<StVKEnergy> pe_;
    std::unique_ptr<PositionCons> pc_;

    // physics model
    const double h_, alpha_, beta_;
    Eigen::VectorXd x_;                                    // store velocity and lagragian multipliers
    zjucad::matrix::matrix<double> disp_;
    zjucad::matrix::matrix<double> fext_;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> M_;
};

#endif
