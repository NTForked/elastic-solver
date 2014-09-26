#include <iostream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/stvk_simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::json_parse::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("stvk.model").c_str(), &nods, &tets);

    unique_ptr<StVKSimulator> sim(new StVK());

    {
        // set fixed points
        vector<size_t> con_nods;
        con_nods.push_back();
        sim->SetFixedPoint();
    }
    {
        // set external force
        for (size_t id = 0; id < xxx; ++id) {
            double f[3] = { };
            sim->SetExternalForce(id, f);
        }
    }

    matrix<double> curr_nods;
    for (size_t frm = 0; frm < 100; ++frm) {
        cerr << "[INFO] this is " << i << " frame.\n";
        sim->Forward();
        curr_nods = nods + sim->disp();
        {
            stringstream ss;
            ss << pt.get<string>("stvk.out_path") << "elastic_" << i << ".vtk";
            ofstream os(ss.str());
            tet2vtk(os, &curr_nods[0], curr_nods.size(2), &tets_[0], tets_.size(2));
        }
    }

    return 0;
}
