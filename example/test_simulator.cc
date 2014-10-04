#include <iostream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/full_simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;
using boost::property_tree::ptree;

void OutputMesh(const matrix<size_t> &tets,
                const matrix<double> &nods,
                const size_t         frm,
                const ptree          &pt);

int main(int argc, char *argv[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("elastic.model").c_str(), &nods, &tets);

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

    // simulate
    matrix<double> curr_nods = nods;
    for (size_t frm = 0; frm < 300; ++frm) {
        cerr << "[INFO] this is " << frm << " frame.\n";
        OutputMesh(tets, curr_nods, frm, pt);
        // twist the model
        if ( frm < 100)
        {
            const double force = 2500;
            matrix<double> dir = curr_nods(colon(), 21) - curr_nods(colon(), 41);
            dir /= norm(dir);
            const double f[3] = {force*dir(0, 0), force*dir(1, 0), force*dir(2, 0)};
            sim->SetExternalForce(21, f);
            for (int id = 22; id <= 41; ++id) {
                matrix<double> dir = curr_nods(colon(), id) - curr_nods(colon(), id - 1);
                dir /= norm(dir);
                const double f[3] = {force*dir(0, 0), force*dir(1, 0), force*dir(2, 0)};
                sim->SetExternalForce(id, f);
            }
        }
        if ( frm == 100 ) {
            sim->ClearExternalForce();
            sim->ClearFixedPoints();
        }
        sim->Forward();
        curr_nods = nods + sim->disp();
    }
    cerr << "[INFO] done.\n";
    return 0;
}

void OutputMesh(const matrix<size_t> &tets,
                const matrix<double> &nods,
                const size_t         frm,
                const ptree          &pt) {
    stringstream ss;
    ss << pt.get<string>("elastic.output_path") << "elastic_" << frm << ".vtk";
    ofstream os(ss.str());
    tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
}
