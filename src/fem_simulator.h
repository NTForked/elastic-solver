#ifndef __ELASTIC_FEM_SIMULATOR_H__
#define __ELASTIC_FEM_SIMULATOR_H__

namespace cj { namespace elastic {

///
/// @brief The FEMSimulator class, Not completed
///

class FEMSimulator {
public :
    FEMSimulator();
    virtual ~FEMSimulator();

    // energy
    void AddEnergy();
    void DeleteEnergy();
    void ClearEnergy();

    // constraint
    void SetFixedPoints();
    void RemoveFixedPoints();
    void ClearAllFixedPoints();

    // constraints
    void AddConstraints();
    void DeleteConstraints();

    // force
    void SetExternalForces();
    void ClearExternalForces();

    // workflow
    virtual void BuildEnvironment();
    virtual void Update();
    virtual void Forward();

    // IO
    virtual void DumpFrame();


protected :
    // Mesh mesh_;
    // Energy energy_;
    // Geometry geom_;
    // LinearSolver solver_;

    // other parameters:
    //   lame coffecients,
    //   time step;
    //   density;
    //   ........
};

class ReducedFEMSimulator : FEMSimulator {
public :
    ReducedFEMSimulator();
    virtual ~ReducedFEMSimulator();

    // basis
    void BuildU();
};

}}
#endif
