#ifndef WARPING_ENERGY_H
#define WARPING_ENERGY_H

#include "math_function.h"

extern "C" {

void axb_energy_(double *val,
                 const double *X,
                 const double *A,
                 const double *B,
                 const double *w);

void axb_energy_jac_(double *gra,
                     const double *X,
                     const double *A,
                     const double *B,
                     const double *w);

void axb_energy_hes_(double *hes,
                     const double *X,
                     const double *A,
                     const double *B,
                     const double *w);
}

namespace cj { namespace elastic {

class WarpingEnergy : public Energy {

};

}}
#endif
