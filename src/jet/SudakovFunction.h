#ifndef SUDAKOVFUNCTION_H
#define SUDAKOVFUNCTION_H

class Matter;

namespace SudakovFunction{

double sudakov_Pgg(double g0, double g1, double loc_c, double E, Matter& matter_obj);
double sud_val_GG(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj);
double sud_z_GG(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj);
double sudakov_Pqq(double q0, double q1, double loc_c, double E, Matter& matter_obj) ;
double sudakov_Pqq_w_M_vac_only(double M, double q0, double q1,double loc_c, double E, Matter& matter_obj);
double sud_val_QQ(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj);
double sud_val_QQ_w_M_vac_only(double M, double h0, double h1,double h2, double loc_d, double E1, Matter& matter_obj);
double sud_z_QQ(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj);
double sud_z_QQ_w_M_vac_only(double M, double cg, double cg1,double loc_e, double l_fac, double E2, Matter& matter_obj);
double sudakov_Pqp(double g0, double g1, double loc_c, double E, Matter& matter_obj);
double sud_val_QP(double h0, double h1, double h2, double loc_d,double E1, Matter& matter_obj);
double sud_z_QP(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj);
double sudakov_Pqg(double g0, double g1, double loc_c, double E, Matter& matter_obj);
double sudakov_Pqg_w_M(double M, double g0, double g1, double loc_c,double E, Matter& matter_obj) ;
double sud_val_QG(double h0, double h1, double h2, double loc_d, double E1, Matter& matter_obj);
double sud_val_QG_w_M(double M, double h0, double h1, double h2,double loc_d, double E1, Matter& matter_obj) ;
double sud_z_QG(double cg, double cg1, double loc_e, double l_fac,double E2, Matter& matter_obj);
double sud_z_QG_w_M(double M, double cg, double cg1, double loc_e,double l_fac, double E2, Matter& matter_obj);
}

#endif