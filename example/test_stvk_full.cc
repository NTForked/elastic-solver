#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>
#include "src/simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace cj::elastic;

void draw(const string &directory, const size_t frm,
          const matrix<size_t> &tets,
          const matrix<double> &nods)
{
    stringstream ss;
    ss << directory << "/frame_" << frm << ".vtk";
    ofstream os(ss.str());
    tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
}

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Usage: " << argv[0] << "model.obj\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./stvk_full");
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json("../../config.json", pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(argv[1], &nods, &tets);

    shared_ptr<StVKSimulator> sol(new StVKSimulator(tets, nods, pt));
    {
        vector<size_t> pin;
        for (size_t id = 0; id <= 20; ++id)
            pin.push_back(id);
        for (size_t id = 42; id <= 65; ++id)
            pin.push_back(id);
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sol->SetFixedPoints(pin, uc);
    }
    sol->SetGravity(3.0);

    matrix<double> curr = nods;
    for (size_t frm = 0; frm < 300; ++frm) {
        cout << "[INFO] frame " << frm << "\n";
        draw("./stvk_full", frm, tets, curr);
        sol->Advance();
        curr = nods + sol->disp();
    }
    cout << "[INFO] done\n";
    return 0;
}
