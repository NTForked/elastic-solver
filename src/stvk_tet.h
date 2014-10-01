#ifndef __CJ_STVK_TET_H__
#define __CJ_STVK_TET_H__

extern "C" {

void stvk_tet_(double        *val,
               const double  *x,
               const double  *Dm,
               const double  *volume,
               const double  *lambda,
               const double  *miu);

void stvk_tet_jac_(double        *jac,
                   const double  *x,
                   const double  *Dm,
                   const double  *volume,
                   const double  *lambda,
                   const double  *miu);

void stvk_tet_hes_(double        *hes,
                   const double  *x,
                   const double  *Dm,
                   const double  *volume,
                   const double  *lambda,
                   const double  *miu);

}

#endif
