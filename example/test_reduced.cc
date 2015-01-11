#include <iostream>
#include <fstream>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>
#include "src/simulator.h"
#include "src/vtk.h"

#define ASSERT(x)                                               \
    do {                                                        \
      if (!(x)) {                                               \
        std::cerr << "[ERROR] Assertion failed at\n";           \
        std::cerr << __FILE__ << " " << __LINE__ << "\n";       \
        exit(0);                                                \
      }                                                         \
    } while(0);

using namespace std;
using namespace zjucad::matrix;
using namespace cj::elastic;
using namespace Eigen;

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

    boost::filesystem::create_directory("./modal_basis");
    {
        std::ofstream os("./modal_basis/model.vtk");
        tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
    }

    shared_ptr<LinearReducedSolver> sol(new LinearReducedSolver(tets, nods, pt));
    sol->Init();
    sol->AddElasticEnergy(1.0);
    std::unordered_set<size_t> fix;
    for (size_t id = 0; id <= 20; ++id) {
        fix.insert(id * 3 + 0);
        fix.insert(id * 3 + 1);
        fix.insert(id * 3 + 2);
    }
    sol->BuildModalBasis(fix);
    cout << sol->lambda_.head(10) << "\n";

    for (size_t i = 0; i < 10; ++i) {
        VectorXd u = sol->U_.col(i);
        ASSERT(nods.size(2) == u.size() / 3);
        itr_matrix<const double*> dx(3, nods.size(2), u.data());
        matrix<double> basis = nods + dx;
        stringstream ss;
        ss << "./modal_basis/basis_" << i << ".vtk";
        ofstream os(ss.str());
        tet2vtk(os, &basis[0], basis.size(2), &tets[0], tets.size(2));
    }
    cout << "[INFO] done\n";
    return 0;
}
