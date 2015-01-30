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
        cerr << "# Usage: ./prog config.json\n";
        return __LINE__;
    }
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("elastic.model").c_str(), &nods, &tets);

    boost::filesystem::create_directory("./reduced_stvk");

    shared_ptr<StvkReducedSolver> sol(new StvkReducedSolver(tets, nods, pt));
    sol->Init();
    sol->AddElasticEnergy(1.0);
    {
        vector<size_t> pin;
        for (size_t id = 0; id <= 20; ++id)
            pin.push_back(id);
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sol->SetPinnedVerts(pin, uc);
    }
    sol->SetGravity(10.0);
    sol->Prepare();

    bool seeBasis = true;
    if ( seeBasis ) {
        for (size_t i = 0; i < sol->U_.cols(); ++i) {
            VectorXd u = sol->U_.col(i);
            itr_matrix<const double*> dx(3, nods.size(2), u.data());
            matrix<double> basis = nods + dx;
            stringstream ss;
            ss << "./reduced_stvk/basis_" << i << ".vtk";
            ofstream os(ss.str());
            tet2vtk(os, &basis[0], basis.size(2), &tets[0], tets.size(2));
        }
    }

    matrix<double> curr = nods;
    for (size_t frm = 0; frm < 300; ++frm) {
        cout << "[INFO] frame " << frm << "\n";
        draw("./reduced_stvk", frm, tets, curr);
        sol->Advance();
        curr = nods + sol->get_disp();
    }
    cout << "[INFO] done\n";
    return 0;
}
