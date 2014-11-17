#include <iostream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/full_simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;
using boost::property_tree::ptree;

#define DRAW(name, frm, nods, tris)                                         \
do {                                                                        \
    stringstream ss;                                                        \
    ss << "./" << name << "_" << frm << ".vtk";                             \
    ofstream os(ss.str());                                                  \
    tet2vtk(os, nods.begin(), nods.size(2), tris.begin(), tris.size(2));    \
} while(0);

int main(int argc, char *argv[])
{
    if ( argc != 4 ) {
        cerr << "# Usage: test config.json input.vtk option\n";
        return __LINE__;
    }

    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(argv[2], &nods, &tets);

    ofstream os("./rest_sofa.vtk");
    tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));

    unique_ptr<StVKSimulator> sim;

    /// CONFIG PHYSICS CONSTRAINTS
    const int sofa_option = atoi(argv[3]);
    if ( sofa_option == 0 ) {
        /// init the simulator
        nods /= 3;
        sim.reset(new StVKSimulator(tets, nods, pt));

        /// set fixed points
        vector<size_t> cons{0, 1, 2, 3, 20, 11, 9, 8, 10};
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sim->SetFixedPoints(cons, uc);

        /// set force
        const double f[3] = {0, 0, -7500};
        vector<size_t> idx{21, 1023, 1024,
                           1037, 1051, 1055,
                           1043, 1033, 1032};
        for (size_t i = 0; i < idx.size(); ++i)
            sim->SetExternalForce(idx[i], f);
    } else if ( sofa_option == 1 ) {
        /// init the simulator
        sim.reset(new StVKSimulator(tets, nods, pt));

        /// set fixed points
        vector<size_t> cons;
        for (size_t i = 0; i < 38; ++i)
            cons.push_back(i);
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sim->SetFixedPoints(cons, uc);

        /// set force
        const double f[3] = {0, 0, -8500};
        vector<size_t> idx{857, 860, 229, 867, 862, 225};
        for (size_t i = 0; i < idx.size(); ++i)
            sim->SetExternalForce(idx[i], f);
        const double ff[3] = {0, 0, -9000};
        sim->SetExternalForce(864, ff);
    }

    /// SIMULATE
    matrix<double> curr_nods = nods;
    for (size_t i = 0; i < 30; ++i) {
        cout << "[INFO] step " << i << "\n";
        DRAW("sofa", i, curr_nods, tets);
        sim->Forward();
        curr_nods = nods + sim->disp();
    }

    cout << "[INFO] done\n";
    return 0;
}
