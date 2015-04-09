#ifndef ANIMATION_EDITOR_H
#define ANIMATION_EDITOR_H

#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace cj { namespace elastic {

class Energy;

/// refer to 'Interactive Editing of Deformable Simulations'
class AnimationEditor {
public:
    typedef zjucad::matrix::matrix<size_t> mati_t;
    typedef zjucad::matrix::matrix<double> matd_t;
    AnimationEditor();
    // IO
    int LoadAnimationSequence(const char *prefix);
    int SaveEditedSequence(const char *prefix);

    // config and init
    int LoadParameters(const boost::property_tree::ptree &pt);
    int Init() {
        BuildLinearElasticEnergy();
        CalcMassMatrix();
        CalcStiffnessMatrix();
    }

private:
    void BuildLinearElasticEnergy();
    void CalcMassMatrix();
    void CalcStiffnessMatrix();

    void BuildModalBasis();

    std::shared_ptr<Energy> e_;
    mati_t tets_;
    std::vector<matd_t> frame_;

    Eigen::SparseMatrix<double> K_, M_;
    Eigen::MatrixXd U_; // modal basis

    Eigen::VectorXd u_; // high dimensional displacement vector
    Eigen::VectorXd z_; // reduced coordinate
};

}}
#endif
