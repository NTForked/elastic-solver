#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include "src/full_simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;
using boost::property_tree::ptree;
using jtf::mesh::face2tet_adjacent;

#define DRAW_TET(name, frm, nods, tris)                                     \
do {                                                                        \
    stringstream ss;                                                        \
    ss << "./" << name << "_" << frm << ".vtk";                             \
    ofstream os(ss.str());                                                  \
    tet2vtk(os, nods.begin(), nods.size(2), tris.begin(), tris.size(2));    \
} while(0);

#define DRAW_TRI(name, frm, nods, tris)                                     \
do {                                                                        \
    stringstream ss;                                                        \
    ss << "./" << name << "_" << frm << ".obj";                             \
    jtf::mesh::save_obj(ss.str().c_str(), tris, nods);                      \
} while(0);

int remove_extra_node(matrix<size_t>     &mesh,
                      matrix<double>     &nods,
                      matrix<size_t>     *l2g = nullptr)
{
    set<size_t> used_node_idx(mesh.begin(), mesh.end());
    if (used_node_idx.size() == nods.size(2))
        return 0;
    matrix<size_t> used_node_mat(used_node_idx.size(), 1);
    std::copy(used_node_idx.begin(), used_node_idx.end(), used_node_mat.begin());

    map<size_t,size_t> p2p;
    matrix<double> new_node(3, used_node_mat.size());

    for (size_t pi = 0; pi < used_node_mat.size(); ++pi) {
        new_node(colon(),pi) = nods(colon(), used_node_mat[pi]);
        p2p[used_node_mat[pi]] = pi;
    }
    if ( l2g )
        *l2g = used_node_mat;
    for (size_t pi = 0; pi < mesh.size(); ++pi)
        mesh[pi] = p2p[mesh[pi]];
    nods = new_node;
    return 0;
}

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
    const int sofa_option = atoi(argv[3]);

    ofstream os("./rest_sofa.vtk");
    tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));

    unique_ptr<StVKSimulator> sim(new StVKSimulator(tets, nods, pt));

    /// CONFIG PHYSICS CONSTRAINTS
    if ( sofa_option == 0 ) {
        /// set fixed points
        vector<size_t> cons{0, 1, 2, 3, 20, 11, 9, 8, 10};
        matrix<double> uc = zeros<double>(3, nods.size(2));
        sim->SetFixedPoints(cons, uc);

        /// set force
        const double f[3] = {0, 0, -3000};
        vector<size_t> idx{21, 1023, 1024,
                           1037, 1051, 1055,
                           1043, 1033, 1032};
        for (size_t i = 0; i < idx.size(); ++i)
            sim->SetExternalForce(idx[i], f);
    } else if ( sofa_option == 1 ) {
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

    /// EXTRACT THE SURFACE
    matrix<size_t> face;
    shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
    jtf::mesh::get_outside_face(*f2t, face);
    jtf::mesh::reorder_face(face);

    /// SIMULATE
    matrix<double> curr_nods = nods;
    for (size_t i = 0; i < 30; ++i) {
        cout << "[INFO] step " << i << "\n";

        DRAW_TET("sofa", i, curr_nods, tets);
        matrix<size_t> tris = face;
        matrix<double> verts = curr_nods;
        remove_extra_node(tris, verts);
        DRAW_TRI("sofa_tri", i, verts, tris);

        sim->Forward();
        curr_nods = nods + sim->disp();
    }

    cout << "[INFO] done\n";
    return 0;
}
