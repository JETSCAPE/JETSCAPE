#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include <boost/math/tools/roots.hpp>
#include "simpleLogger.h"

double renormalization_scale = 2.0; // default

//=============running coupling=================================================
double alpha_s(double Q2, double T){
	double screen_scale2 = std::pow(renormalization_scale*M_PI*T, 2);
	double mu2;
    if (Q2 < 0.){
		mu2 = std::max(-Q2, screen_scale2);
		if (mu2 <= mu2_left) return alpha0;
		else return alpha0 / std::log(mu2/Lambda2);
	}
    else{
		mu2 = std::max(Q2, screen_scale2);
		if (mu2 <= mu2_right) return alpha0;
		else return alpha0 * ( .5 - std::atan(std::log(mu2/Lambda2)/M_PI) / M_PI);
	}
}

/// 					   time duration from last emission
///	LPM factor, x = ------------------------------------------
///							formation time of the process
/// They must be evalutate in the same frame)
double f_LPM(double x){
	return 1. - cos(x);
}

///             Debye mass pointer and class
Debye_mass * t_channel_mD2 = NULL;
Debye_mass::Debye_mass(const unsigned int _type):
TL(0.1), TH(1.0), NT(100), dT((TH-TL)/(NT-1.)),
type(_type), mD2(new double[NT])
{
	if (type==0) {
		LOG_INFO << "# leading order Debye mass";
		// type==0 use self-consistent Debye mass
		for (size_t i=0; i<NT; i++){
			double T = TL+dT*i;
			mD2[i] = pf_g*alpha_s(0., T)*T*T;
		}
	}
	if (type==1) {
		LOG_INFO << "# self-consistent Debye mass";
		// type==1 use self-consistent Debye mass
		for (size_t i=0; i<NT; i++){
			double T = TL+dT*i;
			size_t maxiter=100;
			boost::math::tools::eps_tolerance<double> tol{
			 (std::numeric_limits<double>::digits * 3) / 4};
			try{
				auto result = boost::math::tools::toms748_solve(
					[&T](double x) {return pf_g*alpha_s(-x, T)*T*T - x;},
					0.01, 20., tol, maxiter);
				mD2[i] = .5*(result.first + result.second);
			}
			catch (const std::domain_error&) {
				throw std::domain_error{
				"unable to calculate mD2"};
			}
		}
	}
}

double Debye_mass::get_mD2(double T){
	if (T<TL) T=TL;
	if (T>=TH-dT) T=TH-dT;
	double x = (T-TL)/dT;
	size_t index = std::floor(x);
	double r = x-index;
	return (1.-r)*mD2[index] + r*mD2[index+1];
}

void initialize_mD_and_scale(const unsigned int type, const double scale){
	t_channel_mD2 = new Debye_mass(type);
	renormalization_scale = scale;
	LOG_INFO << "Scale = " << renormalization_scale;
}

/// Q+q --> Q+q
double M2_Qq2Qq(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);
	double At = alpha_s(Q2t, Temp);
	double Q2t_reg = Q2t - mt2;
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*Q2t)/std::pow(Q2t_reg, 2);
	if (result < 0.) return 0.;
	else return result;
}

/// Q+q --> Q+q for radiation process, rightnow, this is the same as Q+q-->Q+q,
/// exept that we don't restrict (-t) > mD2 here, instead, we require later
/// that the emitted gluon has an energy larger than mD, excludes (-t) < epsilon
double M2_Qq2Qq_rad(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);
	double mt2 = t_channel_mD2->get_mD2(Temp);
	double Q2t_reg = Q2t - mt2;
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*Q2t)/std::pow(Q2t_reg, 2);
	if (result < 0.) return 0.;
	else return result;
}

double dX_Qq2Qq_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2, 2);
}

///	Q+g --> Q+g
double M2_Qg2Qg(const double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// Deybe mass (t-channel)
	double mt2 = t_channel_mD2->get_mD2(Temp);

	// define coupling constant for each channel
	double At = alpha_s(Q2t, Temp),
		   Au = alpha_s(Q2u, Temp),
		   As = alpha_s(Q2s, Temp);
	double Q2t_reg = Q2t - mt2;
	double Q2s_reg = Q2s + mt2;
	double Q2u_reg = Q2u>0?(Q2u + mt2):(Q2u-mt2);
	double result = 0.0;
	// t*t
	result += 2.*At*At * Q2s*(-Q2u)/std::pow(Q2t_reg, 2);
	// s*s
	result += c4d9*As*As *
			( Q2s*(-Q2u) + 2.*M2*(Q2s + 2.*M2) ) / std::pow(Q2s_reg, 2);
	// u*u
	result += c4d9*Au*Au *
			( Q2s*(-Q2u) + 2.*M2*(Q2u + 2.*M2) ) / std::pow(Q2u_reg, 2);
	// s*u
	result += c1d9*As*Au * M2*(4.*M2 - Q2t) / Q2s_reg / (-Q2u_reg);
	// t*s
	result += At*As * ( Q2s*(-Q2u) + M2*(Q2s - Q2u) ) / Q2t_reg / Q2s_reg;
    // t*u
	result += -At*Au * ( Q2s*(-Q2u) - M2*(Q2s - Q2u) ) / Q2t_reg / (-Q2u_reg);
	if (result < 0.) return 0.;
	return result*c16pi2;
}

/// Q+g --> Q+g for radiation process, rightnow, this only inclues t-channel
/// We don't restrict (-t) > mD2 here, instead, we require later that the
/// emitted gluon has an energy larger than mD, which exclude (-t) < epsilon
double M2_Qg2Qg_rad(const double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], Temp = p[1], M2 = p[2]*p[2];
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	double At = alpha_s(Q2t, Temp);
	double mt2 = t_channel_mD2->get_mD2(Temp);
	double Q2t_reg = Q2t - mt2;
	double result = c16pi2*2.*At*At * Q2s*(-Q2u)/std::pow(Q2t_reg, 2);
	if (result < 0.) return 0.;
	else return result;
}

double dX_Qg2Qg_dt(const double t, void * params){
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2);
}


/// Q + q --> Q + q + g
double M2_Qq2Qqg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double dt = params[3];
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.)
		return 0.;
	double ymax = std::acosh(Qmax/kt);
	double y = y_norm*ymax;
	// construct k^mu
	fourvec kmu{kt*std::cosh(y), kt, 0, kt*std::sinh(y)};
	double s34 = s - 2.*sqrts*kmu.t();
	double sqrts34 = std::sqrt(s34);
	double Q34 = (s34-M2)/2./sqrts34, E34 = (s34+M2)/2./sqrts34;
	// construct p3, p4 in (2) frame
	double cosphi34 = std::cos(phi34), sinphi34 = std::sin(phi34),
		   sintheta34 = std::sqrt(1.-std::pow(costheta34,2));
	fourvec p3mu{E34, Q34*sintheta34*cosphi34, 
				Q34*sintheta34*sinphi34, Q34*costheta34};
	fourvec p4mu{Q34, -Q34*sintheta34*cosphi34, 
				-Q34*sintheta34*sinphi34, -Q34*costheta34};
	// boost p3, p4 back to 3-body CoM Frame
	double V0 = sqrts - kmu.t();
	double v34[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
	p3mu = p3mu.boost_back(v34[0], v34[1], v34[2]);
	p4mu = p4mu.boost_back(v34[0], v34[1], v34[2]);
	
	// q-perp-vec, q = p2-p4, qperp = -p4perp
	double qx = -p4mu.x(), qy = -p4mu.y();
	double qt2 = qx*qx + qy*qy;
	
	double mD2 = t_channel_mD2->get_mD2(T);
	double kt2 = kt*kt;
	double x = (kmu.t()+kmu.z())/sqrts, xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1.-xbar;
	double basic_denominator = kt2 + x*x*M2 + one_minus_xbar*mD2/2.;
	double tauk = 2.*kmu.t()*one_minus_xbar/basic_denominator;
	// here u is the ratio of the mean-free-path over the formation length
	// mean-free-path \sim mean-free-time*v_HQ,
	// v_HQ = p/E = (s - M^2)/(s + M^2)
	// formation length = tau_k*v_k = tau_k
	double u = dt/tauk*(s-M2)/(s+M2);

	// 2->2
	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	double M2_elastic = M2_Qq2Qq_rad(t, params);

	// 1->2
	double alpha_rad = alpha_s(kt2, T);
	double iD1 = 1./basic_denominator,
	       iD2 = 1./(basic_denominator - 2.*qx*kmu.x()  + qt2);
	double Pg = alpha_rad*std::pow(one_minus_xbar, 2)
				*f_LPM(u)
				*(kt2*std::pow(iD1-iD2, 2.) + qt2*std::pow(iD2,2) + 2.*kmu.x()*qx*iD2*(iD1-iD2));
	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}

/// Q + g --> Q + g + g
double M2_Qg2Qgg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double dt = params[3];
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.)
		return 0.;
	double ymax = std::acosh(Qmax/kt);
	double y = y_norm*ymax;
	// construct k^mu
	fourvec kmu{kt*std::cosh(y), kt, 0, kt*std::sinh(y)};
	double s34 = s - 2.*sqrts*kmu.t();
	double sqrts34 = std::sqrt(s34);
	double Q34 = (s34-M2)/2./sqrts34, E34 = (s34+M2)/2./sqrts34;
	// construct p3, p4 in (2) frame
	double cosphi34 = std::cos(phi34), sinphi34 = std::sin(phi34),
		   sintheta34 = std::sqrt(1.-std::pow(costheta34,2));
	fourvec p3mu{E34, Q34*sintheta34*cosphi34, 
				Q34*sintheta34*sinphi34, Q34*costheta34};
	fourvec p4mu{Q34, -Q34*sintheta34*cosphi34, 
				-Q34*sintheta34*sinphi34, -Q34*costheta34};
	// boost p3, p4 back to 3-body CoM Frame
	double V0 = sqrts - kmu.t();
	double v34[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
	p3mu = p3mu.boost_back(v34[0], v34[1], v34[2]);
	p4mu = p4mu.boost_back(v34[0], v34[1], v34[2]);
	
	// q-perp-vec, q = p2-p4, qperp = -p4perp
	double qx = -p4mu.x(), qy = -p4mu.y();
	double qt2 = qx*qx + qy*qy;
	
	double mD2 = t_channel_mD2->get_mD2(T);
	double kt2 = kt*kt;
	double x = (kmu.t()+kmu.z())/sqrts, xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1.-xbar;
	double basic_denominator = kt2 + x*x*M2 + one_minus_xbar*mD2/2.;
	double tauk = 2.*kmu.t()*one_minus_xbar/basic_denominator;
	// here u is the ratio of the mean-free-path over the formation length
	// mean-free-path \sim mean-free-time*v_HQ,
	// v_HQ = p/E = (s - M^2)/(s + M^2)
	// formation length = tau_k*v_k = tau_k
	double u = dt/tauk*(s-M2)/(s+M2);

	// 2->2
	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	double M2_elastic = M2_Qg2Qg_rad(t, params);

	// 1->2
	double alpha_rad = alpha_s(kt2, T);
	double iD1 = 1./basic_denominator,
	       iD2 = 1./(basic_denominator - 2.*qx*kmu.x()  + qt2);
	double Pg = alpha_rad*std::pow(one_minus_xbar, 2)
				*f_LPM(u)
				*(kt2*std::pow(iD1-iD2, 2.) + qt2*std::pow(iD2,2) + 2.*kmu.x()*qx*iD2*(iD1-iD2));
	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}


//=============Basic for 3->2===========================================
double Ker_Qqg2Qq(double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta24 = x_[0], phi24 = x_[1];
	if (costheta24<=-1. || costheta24>=1. || phi24 <=0. || phi24 >=2.*M_PI) return 0.;
	double sintheta24 = std::sqrt(1. - costheta24*costheta24), sinphi24 = std::sin(phi24), cosphi24 = std::cos(phi24);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12k, T, M, ...
	double E2 = params[3];
	double E4 = params[4];
	double TwoE2E4 = 2.*E2*E4;
	double k = params[5];
	/// The energy of the in coming gluon must have energy larger than mD in
	/// CoM frame
		double cos2 = params[6];
		double sin2 = -std::sqrt(1. - cos2*cos2);
	  	double cosk = params[7];
	  	double sink = std::sqrt(1. - cosk*cosk);
		double x2M2 = params[8];
		double coeff_mg2 = params[9];
		// 2->2
		double t = TwoE2E4 * (costheta24 - 1);
		double M2_elastic = M2_Qq2Qq_rad(t, params);

		// 1->2
		double p2x = E2*sin2, p2z = E2*cos2;
		double kx = k*sink, kz = k*cosk;
		double p4x = E4*(cos2*sintheta24*cosphi24 + sin2*costheta24),
		     p4y = E4*sintheta24*sinphi24,
		     p4z = E4*(-sin2*sintheta24*cosphi24 + cos2*costheta24);
		double qx = p2x - p4x, qy = - p4y, qz = p2z - p4z;
		double q_project = (qx*p4x + qy*p4y + qz*p4z)/E4,
		     k_project = (kx*p4x + kz*p4z)/E4;
		double qt2 = qx*qx+qy*qy+qz*qz - std::pow(q_project,2);
		double kt2 = kx*kx + kz*kz - std::pow(k_project, 2);
		double two_kt_dot_qt = kx*qx + kz*qz - k_project*q_project;
		double iD1 = 1./(kt2 + x2M2 + coeff_mg2);
		double iD2 = 1./(kt2 + qt2 + two_kt_dot_qt + x2M2 + coeff_mg2);
		double Pg = kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2)
	              - two_kt_dot_qt*(iD1-iD2)*iD2;

		// 2->3 = 2->2 * 1->2
		return M2_elastic*Pg/16.;
}

double Ker_Qgg2Qg(double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta24 = x_[0], phi24 = x_[1];
	if (costheta24<=-1. || costheta24>=1. || phi24 <=0. || phi24 >=2.*M_PI) return 0.;
	double sintheta24 = std::sqrt(1. - costheta24*costheta24), sinphi24 = std::sin(phi24), cosphi24 = std::cos(phi24);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12k, T, M, ...
	double E2 = params[3];
	double E4 = params[4];
	double TwoE2E4 = 2.*E2*E4;
	double k = params[5];
	/// The energy of the in coming gluon must have energy larger than mD in
	/// CoM frame
		double cos2 = params[6];
		double sin2 = -std::sqrt(1. - cos2*cos2);
		double cosk = params[7];
		double sink = std::sqrt(1. - cosk*cosk);
		double x2M2 = params[8];
		double coeff_mg2 = params[9];
		// 2->2
		double t = TwoE2E4 * (costheta24 - 1);
		double M2_elastic = M2_Qg2Qg_rad(t, params);

		// 1->2
		double p2x = E2*sin2, p2z = E2*cos2;
		double kx = k*sink, kz = k*cosk;
		double p4x = E4*(cos2*sintheta24*cosphi24 + sin2*costheta24),
		   p4y = E4*sintheta24*sinphi24,
		   p4z = E4*(-sin2*sintheta24*cosphi24 + cos2*costheta24);
		double qx = p2x - p4x, qy = - p4y, qz = p2z - p4z;
		double q_project = (qx*p4x + qy*p4y + qz*p4z)/E4,
		   k_project = (kx*p4x + kz*p4z)/E4;
		double qt2 = qx*qx+qy*qy+qz*qz - std::pow(q_project,2);
		double kt2 = kx*kx + kz*kz - std::pow(k_project, 2);
		double two_kt_dot_qt = kx*qx + kz*qz - k_project*q_project;
		double iD1 = 1./(kt2 + x2M2 + coeff_mg2);
		double iD2 = 1./(kt2 + qt2 + two_kt_dot_qt + x2M2 + coeff_mg2);
		double Pg = kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2)
				- two_kt_dot_qt*(iD1-iD2)*iD2;

		// 2->3 = 2->2 * 1->2
		return M2_elastic*Pg/16.;
}
