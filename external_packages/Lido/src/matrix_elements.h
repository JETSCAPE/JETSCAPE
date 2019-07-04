#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>
#include <vector>
#include "lorentz.h"

////////// 1 <--> 2 //////////////////
// diffusion induced process
double LGV_q2qg(const double * x_, void *params_); 
double LGV_qg2q(const double * x_, void *params_); 
 
double LGV_g2gg(const double * x_, void *params_); 
double LGV_gg2g(const double * x_, void *params_);  

double LGV_g2qqbar(const double * x_, void *params_); 

////////// 2 <--> 2 //////////////////
double M2_gq2gq(const double t, void * params);
double dX_gq2gq_dt(const double t, void * params);

double M2_gg2gg(const double t, void * params);
double dX_gg2gg_dt(const double t, void * params);

double M2_Qq2Qq(const double t, void * params);
double dX_Qq2Qq_dt(const double t, void * params);

double M2_Qg2Qg(const double t, void * params);
double dX_Qg2Qg_dt(const double t, void * params);
double dX_Qg2Qg_dt_full(const double t, void * params);

// *** this one includes s, t, u channels and interferences
double M2_Qg2Qg_full(const double t, void * params);

////////// 2 <--> 3 //////////////////
double M2_Qq2Qqg(const double * x_, void * params_);
double M2_Qg2Qgg(const double * x_, void * params_);

double M2_Qqg2Qq(const double * x_,  void * params_);
double M2_Qgg2Qg(const double * x_, void * params_);

double M2_gq2gqg(const double * x_, void * params_);
double M2_gg2ggg(const double * x_, void * params_);

double M2_gqg2gq(const double * x_,  void * params_);
double M2_ggg2gg(const double * x_, void * params_);

double M2_gq2qqqbar(const double * x_, void * params_);
double M2_gg2qgqbar(const double * x_, void * params_);
#endif
