#ifndef __ELASTIC_ELASTIC_ENERGY_H__
#define __ELASTIC_ELASTIC_ENERGY_H__

#include <zjucad/matrix/matrix.h>

namespace cj { namespace elastic {

class Energy;

Energy* BuildElasticEnergy(const zjucad::matrix::matrix<size_t> &tets,
                           const zjucad::matrix::matrix<double> &nods,
                           const double lambda,
                           const double miu,
                           const std::string type,
                           const double w = 1);
}}

#endif
