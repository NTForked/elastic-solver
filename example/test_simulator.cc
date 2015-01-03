#include <iostream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/simulator.h"
#include "src/vtk.h"

#define _COMPRESS_MODEL_

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;
using boost::property_tree::ptree;

void OutputMesh(const matrix<size_t> &tets,
                const matrix<double> &nods,
                const size_t         frm,
                const ptree          &pt);

bool kCOMPRESS = false;
bool kTWIST = false;

int main(int argc, char *argv[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> TETS;
    matrix<double> NODS;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("elastic.model").c_str(), &NODS, &TETS);

//    matrix<size_t> tets(4, TETS.size(2) * 2);
//    matrix<double> nods(3, NODS.size(2) * 2);
//    tets(colon(), colon(0, TETS.size(2) - 1)) = TETS;
//    tets(colon(), colon(TETS.size(2), 2 * TETS.size(2) - 1)) = TETS + NODS.size(2) * ones<size_t>(4, TETS.size(2));

//    nods(colon(), colon(0, NODS.size(2) - 1)) = NODS;
//    nods(colon(), colon(NODS.size(2), 2 * NODS.size(2) - 1)) = NODS + ones<double>(3, NODS.size(2));
    matrix<size_t> tets = TETS;
    matrix<double> nods = NODS;

    unique_ptr<StVKSimulator> sim(new StVKSimulator(tets, nods, pt));
    {
        // set fixed points
        vector<size_t> cons_nodes;
        for (int id = 0; id <= 20; ++id)
            cons_nodes.push_back(id);
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sim->SetFixedPoints(cons_nodes, uc);
    }

    kCOMPRESS = pt.get<bool>("elastic.compress");
    kTWIST = pt.get<bool>("elastic.twist");

    if ( kCOMPRESS ) {
        // set compression force
        const double f[3] = {0, 0, -5000};
        for (size_t idx = 66; idx <= 89; ++idx) {
            sim->SetExternalForce(idx, f);
        }
    }

    // simulate
    matrix<double> curr_nods = nods;
    for (size_t frm = 0; frm < 300; ++frm) {
        cerr << "[INFO] this is " << frm << " frame.\n";
        OutputMesh(tets, curr_nods, frm, pt);

        if ( kTWIST ) { // twist the model
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
        }

        if ( kCOMPRESS ) {
            if ( frm == 50 )
                sim->ClearExternalForce();
        }

        sim->Advance();
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
