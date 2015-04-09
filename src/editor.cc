#include "editor.h"

#include <cstdio>
#include <jtflib/mesh/io.h>

#include "math_function.h"
#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace cj { namespace elastic {

int AnimationEditor::LoadAnimationSequence(const char *prefix) {
    char file[256];
    mati_t tets;
    matd_t nods;
    size_t frm;
    for (frm = 0; ; ++frm) {
        sprintf(file, "%s_%zu.vtk", prefix, frm);
        if ( !jtf::mesh::tet_mesh_read_from_vtk(file, &nods, &tets) ) {
            frame_.push_back(nods);
        } else {
            break;
        }
    }
    if ( frm == 0 ) {
        cerr << "# INFO: null input\n";
        return __LINE__;
    }
    return 0;
}

}}
