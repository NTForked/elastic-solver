#include <iostream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/stvk_simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;

int main(int argc, char *argv[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("stvk.model").c_str(), &nods, &tets);

    unique_ptr<StVKSimulator> sim(new StVKSimulator(tets, nods, pt));
    {
        // set fixed points
        vector<size_t> cons_nodes;
        for (int id = 0; id <= 20; ++id)
            cons_nodes.push_back(id);
        for (int id = 42; id <= 65; ++id)
            cons_nodes.push_back(id);
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sim->SetFixedPoints(cons_nodes, uc);
    }
    {
        // set external force
        for (size_t id = 66; id <= 89; ++id) {
            double f[3] = {0, 0, -6000};
            sim->SetExternalForce(id, f);
        }
    }
    // simulate
    matrix<double> curr_nods = nods;
    for (size_t frm = 0; frm < 200; ++frm) {
        cerr << "[INFO] this is " << frm << " frame.\n";
        {
            stringstream ss;
            ss << pt.get<string>("stvk.output_path") << "elastic_" << frm << ".vtk";
            ofstream os(ss.str());
            tet2vtk(os, &curr_nods[0], curr_nods.size(2), &tets[0], tets.size(2));
        }
        if ( frm == 30 )
            sim->ClearExternalForce();
        sim->Forward();
        curr_nods = nods + sim->disp();
    }
    cerr << "[INFO] done.\n";
    return 0;
}
