#ifndef ELASTIC_CONSTITUTIVE_H
#define ELASTIC_CONSTITUTIVE_H

extern "C" {

void linear_elastic_tet_(double        *val,
                         const double  *x,
                         const double  *Dm,
                         const double  *volume,
                         const double  *lambda,
                         const double  *miu);

void linear_elastic_tet_jac_(double        *jac,
                             const double  *x,
                             const double  *Dm,
                             const double  *volume,
                             const double  *lambda,
                             const double  *miu);

void linear_elastic_tet_hes_(double        *hes,
                             const double  *x,
                             const double  *Dm,
                             const double  *volume,
                             const double  *lambda,
                             const double  *miu);



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


void neohookean_tet_(double        *hes,
                     const double  *x,
                     const double  *Dm,
                     const double  *volume,
                     const double  *lambda,
                     const double  *miu);

void neohookean_tet_jac_(double        *hes,
                         const double  *x,
                         const double  *Dm,
                         const double  *volume,
                         const double  *lambda,
                         const double  *miu);

void neohookean_tet_hes_(double        *hes,
                         const double  *x,
                         const double  *Dm,
                         const double  *volume,
                         const double  *lambda,
                         const double  *miu);

void fung_tet_(double        *val,
               const double  *x,
               const double  *Dm,
               const double  *volume,
               const double  *lambda,
               const double  *miu,
               const double  *lambda1,
               const double  *miu1,
               const double  *c);

void fung_tet_jac_(double        *jac,
                   const double  *x,
                   const double  *Dm,
                   const double  *volume,
                   const double  *lambda,
                   const double  *miu,
                   const double  *lambda1,
                   const double  *miu1,
                   const double  *c);

void fung_tet_hes_(double        *hes,
                   const double  *x,
                   const double  *Dm,
                   const double  *volume,
                   const double  *lambda,
                   const double  *miu,
                   const double  *lambda1,
                   const double  *miu1,
                   const double  *c);

}

void appro_corotational_tet_(double *val,
                             const double *x,
                             const double *Dm,
                             const double *volume,
                             const double *lambda,
                             const double *miu);

void appro_corotational_tet_jac_(double *jac,
                                 const double *x,
                                 const double *Dm,
                                 const double *volume,
                                 const double *lambda,
                                 const double *miu);

void appro_corotational_tet_hes_(double *hes,
                                 const double *x,
                                 const double *Dm,
                                 const double *volume,
                                 const double *lambda,
                                 const double *miu);

#endif
