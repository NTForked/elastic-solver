#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>

#include <jtflib/mesh/io.h>
#include "src/simulator.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Usage: ./program config.json\n";
        return __LINE__;
    }
    boost::property_tree::ptree pt;
    boost::property_tree::json_parser::read_json(argv[1], pt);

    matrix<size_t> tets;
    matrix<double> nods;
    jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("elastic.model").c_str(), &nods, &tets);

    shared_ptr<ReducedSolver> sol(new ReducedSolver(tets, nods, pt));
    sol->Init();
    sol->AddElasticEnergy(1.0);
    sol->BuildModalBasis();
    cout << sol->freq_.head(10) << "\n";
    cout << "[INFO] done\n";
    return 0;
}
