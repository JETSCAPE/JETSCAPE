#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include <boost/math/tools/roots.hpp>
#include "simpleLogger.h"

double renormalization_scale; // default

//=============running coupling=================================================
double alpha_s(double Q2, double T){
	double screen_scale2 = std::pow(renormalization_scale*M_PI*T, 2);
	double mu2;
    if (Q2 < 0.) mu2 = std::max(-Q2, screen_scale2);
		else mu2 = std::max(Q2, screen_scale2);

		if (mu2 <= mu2_left) return alpha0;
		else return alpha0 / std::log(mu2/Lambda2);
}

/// 					   time duration from last emission
///	LPM factor, x = ------------------------------------------
///							formation time of the process
/// They must be evalutate in the same frame)
double f_LPM(double x){
	return 2*(1. - std::cos(x));
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
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*Q2t)/Q2t_reg/(Q2t-Lambda2);
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
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*Q2t)/Q2t_reg/(Q2t-Lambda2);
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
	result += 2.*At*At * Q2s*(-Q2u)/Q2t_reg/(Q2t-Lambda2);
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
	double result = c16pi2*2.*At*At * Q2s*(-Q2u)/Q2t_reg/(Q2t-Lambda2);
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
	double kt_qt2 = kt2 - 2.*qx*kmu.x() + qt2;
	double x = (kmu.t()+kmu.z())/sqrts,
	       xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1.-xbar;

	double iD1 = 1./(kt2 + x*x*M2 + (1.-xbar)*mD2/2.),
	       iD2 = 1./(kt_qt2 + x*x*M2 + (1.-xbar)*mD2/2.);
	double tauk = 2.*one_minus_xbar*kmu.t()*iD1;
	// here u is the ratio of the mean-free-path over the formation length
	// mean-free-path \sim mean-free-time*v_HQ,
	// v_HQ = p/E = (s - M^2)/(s + M^2)
	// formation length = tau_k*v_k = tau_k
	double u = dt/tauk*(s-M2)/(s+M2);

	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	double M2_elastic = M2_Qq2Qq_rad(t, params);

	double alpha_rad = alpha_s(kt2, T);
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
    double kt_qt2 = kt2 - 2.*qx*kmu.x() + qt2;
    double x = (kmu.t()+kmu.z())/sqrts,
           xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
    double one_minus_xbar = 1.-xbar;
    double iD1 = 1./(kt2 + x*x*M2 + (1.-xbar)*mD2/2.),
           iD2 = 1./(kt_qt2 + x*x*M2 + (1.-xbar)*mD2/2.);
    double tauk = 2.*one_minus_xbar*kmu.t()*iD1;
	double u = dt/tauk*(s-M2)/(s+M2);

        double t = -2.*Qmax*(p4mu.t()+p4mu.z());
        double M2_elastic = M2_Qg2Qg_rad(t, params);

        double alpha_rad = alpha_s(kt2, T);
        double Pg = alpha_rad*std::pow(one_minus_xbar, 2)
      *f_LPM(u)
      *(kt2*std::pow(iD1-iD2, 2.) + qt2*std::pow(iD2,2) + 2.*kmu.x()*qx*iD2*(iD1-iD2));

	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}


//=============Basic for 3->2===========================================
// Not CoM frame of  1,2,k, but CoM frame of 1 and 2
// sampled final states within 3+4 CoM frame
double prefix_3to2(double s, double s12, double s1k, double dt, double M, double T){
	double M2 = M*M;
	double sqrts = std::sqrt(s);
	double sqrts12 = std::sqrt(s12);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double kt2 = k*k*(1. - costhetak*costhetak);
	double kz = k*costhetak;
	double x = (k+kz)/(sqrts12+k+kz), xbar = (k+std::abs(kz))/(sqrts12+k+std::abs(kz));
	double x2M2 = x*x*M2;
	double mD2 = t_channel_mD2->get_mD2(T);
	double tauk = 2.*(1.-xbar)*k/(kt2 + x2M2 + (1.-xbar)*mD2/2.);
	return f_LPM(dt/tauk*(s12-M2)/(s12+M2));
}

double Ker_Qqg2Qq(const double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta34 = x_[0], phi34 = x_[1];

	if (costheta34<=-1. || costheta34>=1. || phi34 <=0. || phi34 >=2.*M_PI) return 0.;
	double sintheta34 = std::sqrt(1. - costheta34*costheta34),
		sinphi34 = std::sin(phi34), cosphi34 = std::cos(phi34);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12, T, M, k, costhetak
	double s = params[0];
	double T = params[1];
	double M = params[2];
	double M2 = M*M;
	double s12 = params[3]*(s-M2) + M2; // 0 < params[3] = xinel = (s12-M2)/(s-M2) < 1
	double s1k = (params[4]*(1-s12/s)*(1-M2/s12)+M2/s12)*s; // 0 < params[4] = yinel = (s1k/s-M2/s12)/(1-s12/s)/(1-M2/s12) < 1
	double sqrts12 = std::sqrt(s12);
	double sqrts = std::sqrt(s);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double E3 = (s+M2)/2./sqrts, p3 = (s-M2)/2./sqrts;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double sinthetak = std::sqrt(1. - costhetak*costhetak);
	double kt = k*sinthetak;
	double kt2 = kt*kt;
	double kz = k*costhetak;
	double x = (k+kz)/(sqrts12+k+kz), xbar = (k+std::abs(kz))/(sqrts12+k+std::abs(kz));
  double mD2 = t_channel_mD2->get_mD2(T);

	// get final state
	fourvec Ptot{sqrts12+k, kt, 0., kz};
	double vcom[3] = {Ptot.x()/Ptot.t(), Ptot.y()/Ptot.t(),Ptot.z()/Ptot.t()};
	// final state in 34-com frame
	fourvec p3mu{E3, p3*sintheta34*cosphi34, p3*sintheta34*sinphi34, p3*costheta34};
	fourvec p4mu{p3, -p3*sintheta34*cosphi34, -p3*sintheta34*sinphi34, -p3*costheta34};

	// boost final state back to 12-com frame
	p3mu = p3mu.boost_back(vcom[0], vcom[1], vcom[2]);
	p4mu = p4mu.boost_back(vcom[0], vcom[1], vcom[2]);

	// 2->2 part
	fourvec qmu{p1-p4mu.t(), -p4mu.x(), -p4mu.y(), -p1-p4mu.z()};
	double t = -2.*p1*(p4mu.t()+p4mu.z());

	double qt2 = std::pow(qmu.x(),2) + std::pow(qmu.y(),2);

	double new_params[3] = {s, T, M};
	double M2_elastic = M2_Qq2Qq_rad(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + x2M2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + x2M2 + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	//double anti_prefix = kt2 + x2M2 + (1.-xbar)*mD2/2.;
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

double Ker_Qgg2Qg(const double * x_, void * params_){
	// unpack variables costheta42 = x_[0]
	double costheta34 = x_[0], phi34 = x_[1];

	if (costheta34<=-1. || costheta34>=1. || phi34 <=0. || phi34 >=2.*M_PI) return 0.;
	double sintheta34 = std::sqrt(1. - costheta34*costheta34),
		sinphi34 = std::sin(phi34), cosphi34 = std::cos(phi34);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12, T, M, k, costhetak
	double s = params[0];
	double T = params[1];
	double M = params[2];
	double M2 = M*M;
	double s12 = params[3]*(s-M2) + M2; // 0 < params[3] = xinel = (s12-M2)/(s-M2) < 1
	double s1k = (params[4]*(1-s12/s)*(1-M2/s12)+M2/s12)*s; // 0 < params[4] = yinel = (s1k/s-M2/s12)/(1-s12/s)/(1-M2/s12) < 1
	double sqrts12 = std::sqrt(s12);
	double sqrts = std::sqrt(s);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double E3 = (s+M2)/2./sqrts, p3 = (s-M2)/2./sqrts;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double sinthetak = std::sqrt(1. - costhetak*costhetak);
	double kt = k*sinthetak;
	double kt2 = kt*kt;
	double kz = k*costhetak;
	double x = (k+kz)/(sqrts12+k+kz), xbar = (k+std::abs(kz))/(sqrts12+k+std::abs(kz));
    double mD2 = t_channel_mD2->get_mD2(T);

	// get final state
	fourvec Ptot{sqrts12+k, kt, 0., kz};
	double vcom[3] = {Ptot.x()/Ptot.t(), Ptot.y()/Ptot.t(),Ptot.z()/Ptot.t()};
	// final state in 34-com frame
	fourvec p3mu{E3, p3*sintheta34*cosphi34, p3*sintheta34*sinphi34, p3*costheta34};
	fourvec p4mu{p3, -p3*sintheta34*cosphi34, -p3*sintheta34*sinphi34, -p3*costheta34};

	// boost final state back to 12-com frame
	p3mu = p3mu.boost_back(vcom[0], vcom[1], vcom[2]);
	p4mu = p4mu.boost_back(vcom[0], vcom[1], vcom[2]);

	// 2->2 part
	fourvec qmu{p1-p4mu.t(), -p4mu.x(), -p4mu.y(), -p1-p4mu.z()};
	double t = -2.*p1*(p4mu.t()+p4mu.z());

	double qt2 = std::pow(qmu.x(),2) + std::pow(qmu.y(),2);

	double new_params[3] = {s, T, M};
	double M2_elastic = M2_Qg2Qg_rad(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + x2M2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + x2M2 + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	//double anti_prefix = kt2+x2M2 + (1.-xbar)*mD2/2.;
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}
