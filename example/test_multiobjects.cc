#include <iostream>
#include <boost/property_tree/json_parser.hpp>
#include <jtflib/mesh/io.h>

#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using boost::property_tree::ptree;

int main(int argc, char *argv[])
{
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("elastic.model").c_str(), &nods, &tets);

    matrix<size_t> TETS(4, tets.size(2) * 2);
    matrix<double> NODS(3, nods.size(2) * 2);

    TETS(colon(), colon(0, tets.size(2) - 1)) = tets;
    TETS(colon(), colon(tets.size(2), 2 * tets.size(2) - 1)) = tets + nods.size(2) * ones<size_t>(4, tets.size(2));

    NODS(colon(), colon(0, nods.size(2) - 1)) = nods;
    NODS(colon(), colon(nods.size(2), 2 * nods.size(2) - 1)) = nods + 0.8 * ones<double>(3, nods.size(2));

    stringstream ss;
    ss << pt.get<string>("elastic.output_path") << "mutibody.vtk";
    ofstream os(ss.str());
    tet2vtk(os, &NODS[0], NODS.size(2), &TETS[0], TETS.size(2));

    cerr << "[INFO] done.\n";
    return 0;
}
