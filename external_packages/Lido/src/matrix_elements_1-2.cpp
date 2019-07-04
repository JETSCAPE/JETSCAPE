#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"

/////////////////////// SPLITTING ////////////////////////

// Diffusion-induced gluon radiation process q -> q + g
double LGV_q2qg(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double x = x_[0]; // k/p
	double y = x_[1]; // kT/k0

    double k0 = x*pabs;
    double kT = y*k0;
	double kT2 = kT*kT;
    
	double mg2 = t_channel_mD2->get_mD2(T)/2.;
	double x0 = k0/E;
	// No dead cone
	double Jacobian = 2*k0*kT;
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * P_q2qg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian;
    return dR_dxdy;
}

// Diffusion-induced gluon radiation process g -> g + g
double LGV_g2gg(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double x = x_[0]; // k/p
	double y = x_[1]; // kT/k0

    double k0 = x*pabs;
    double kT = y*k0;
	double kT2 = kT*kT;
    
	double mg2 = t_channel_mD2->get_mD2(T)/2.;
	double x0 = k0/E;
	if (x0 > 0.5) return 0.0; // aviod final state double counting
	// No dead cone
	double Jacobian = 2*k0*kT;
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * P_g2gg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian;
    return dR_dxdy;
}

// Diffusion-induced gluon splitting process g -> q + qbar
double LGV_g2qqbar(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double x = x_[0]; // k/p
	double y = x_[1]; // kT/k0

    double k0 = x*pabs;
    double kT = y*k0;
	double kT2 = kT*kT;
    
	double mg2 = t_channel_mD2->get_mD2(T)/2.;
	double x0 = k0/E;
	// No dead cone
	double Jacobian = 2*k0*kT;
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * P_g2qq(x0) * CF / 4.
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian;
    return dR_dxdy;
}

////////////////// ABSORPTION ////////////////////////


// Diffusion-induced gluon absorption process q + g -> q
double LGV_qg2q(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double xp = x_[0];// tanh(x)
	double x = std::atanh(xp); // kz/pabs
	double y = x_[1];// kT/k

	double kz = x*pabs;
	double k0 = std::abs(kz)/std::sqrt(1.001-y*y);
	double kT = y*k0;

	double kT2 = kT*kT;
	double x0 = k0/(k0+E); // k0/(E0)

	double mg2 = t_channel_mD2->get_mD2(T)/2.;

	// No dead cone
	double Jacobian = 2*(k0+E)*kT/(1-xp*xp);
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * x0 * P_q2qg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian * std::exp(-k0/T);
    return dR_dxdy;
}

// Diffusion-induced gluon absorption process g + g -> g
double LGV_gg2g(const double * x_, void * params_){
    double *params = static_cast<double*>(params_);
    double E = params[0];
    double T = params[1];
    double M = params[2];
	double pabs = std::sqrt(E*E-M*M);
    
	double xp = x_[0];// tanh(x)
	double x = std::atanh(xp); // kz/pabs
	double y = x_[1];// kT/k

	double kz = x*pabs;
	double k0 = std::abs(kz)/std::sqrt(1.001-y*y);
	double kT = y*k0;

	double kT2 = kT*kT;
	double x0 = k0/(k0+E); // k0/(E0)

	double mg2 = t_channel_mD2->get_mD2(T)/2.;

	// No dead cone
	double Jacobian = 2*(k0+E)*kT/(1-xp*xp);
    double dR_dxdy = alpha_s(kT2, T)/(2.*M_PI) * x0 * P_g2gg(x0)
                     * 1./std::pow(kT2+mg2, 2)
                     * Jacobian * std::exp(-k0/T);
    return dR_dxdy;
}
