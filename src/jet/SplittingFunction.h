#ifndef SPLITTINGFUNCTION_H
#define SPLITTINGFUNCTION_H

class Matter;

namespace SplittingFunction{
	double P_z_gg_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2, Matter& matter_obj);
	double P_z_qq_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2, Matter& matter_obj);
	double P_z_qq_int_w_M_vac_only(double M, double cg, double cg1,double loc_e, double cg3, double l_fac,double E2,Matter& matter_obj);
	double P_z_qp_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2,Matter& matter_obj);
	double P_z_qg_int(double cg, double cg1, double loc_e, double cg3,double l_fac, double E2,Matter& matter_obj);
	double P_z_qg_int_w_M(double M, double cg, double cg1, double loc_e,double cg3, double l_fac, double E2,Matter& matter_obj); 
}

#endif