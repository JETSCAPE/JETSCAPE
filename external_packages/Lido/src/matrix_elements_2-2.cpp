#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"

///	g+g --> g+g
double M2_gg2gg(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = 0.;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp);

	double result = 72.*M_PI*M_PI*At*At*(-Q2s*Q2u/Q2t/Q2t);

	if (result < 0.) return 0.;
	return result;
}
double dX_gg2gg_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = 0.;
	return M2_gg2gg(t, params)/c16pi/(std::pow(s-M2, 2));
}
///	g+q --> g+q
double M2_gq2gq(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = 0.;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp);

	double result = At*At*64.*M_PI*M_PI/9.*(Q2s*Q2s+Q2u*Q2u)*( 9./4./Q2t/Q2t);
	if (result < 0.) return 0.;
	return result;
}
double dX_gq2gq_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = 0.;
	return M2_gq2gq(t, params)/c16pi/(std::pow(s-M2, 2));
}


/// Q+q --> Q+q
double M2_Qq2Qq(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);

	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*t)/Q2t/Q2t;
	if (result < 0.) return 0.;
	else return result;
}

double dX_Qq2Qq_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2, 2);
}

double M2_Qg2Qg(const double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);

	double result = c16pi2*2.*At*At * Q2s*(-Q2u)/Q2t/Q2t;

	if (result < 0.) return 0.;
	else return result;
}

double dX_Qg2Qg_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2);
}

///	full Q+g --> Q+g |s+u+t channel|^2
double M2_Qg2Qg_full(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = std::min(t, -cut*mt2), Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp),
		   Au = alpha_s(Q2u, Temp),
		   As = alpha_s(Q2s, Temp);
	double Q2s_reg = Q2s + mt2;
	double Q2u_reg = Q2u>0?(Q2u + mt2):(Q2u-mt2);

	double result = 0.0;
	// t*t
	result += 2.*At*At * Q2s*(-Q2u)/Q2t/Q2t;
	// s*s
	result += c4d9*As*As *
			( Q2s*(-Q2u) + 2.*M2*(Q2s + 2.*M2) ) / std::pow(Q2s_reg, 2);
	// u*u
	result += c4d9*Au*Au *
			( Q2s*(-Q2u) + 2.*M2*(Q2u + 2.*M2) ) / std::pow(Q2u_reg, 2);
	// s*u
	result += c1d9*As*Au * M2*(4.*M2 - t) / Q2s_reg / (-Q2u_reg);
	// t*s
	result += At*As * ( Q2s*(-Q2u) + M2*(Q2s - Q2u) ) / Q2t / Q2s_reg;
    // t*u
	result += -At*Au * ( Q2s*(-Q2u) - M2*(Q2s - Q2u) ) / Q2t / (-Q2u_reg);
	if (result < 0.) return 0.;
	return result*c16pi2;
}

double dX_Qg2Qg_dt_full(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg_full(t, params)/c16pi/std::pow(s-M2, 2);
}



