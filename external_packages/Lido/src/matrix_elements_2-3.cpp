#include <cmath>
#include <iostream>
#include "predefine.h"
#include "matrix_elements.h"
#include "lorentz.h"
#include "simpleLogger.h"



////////////////////// Hard Quark //////////////////////////
/// Q + q --> Q + q + g
double M2_Qq2Qqg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>=.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	
	if(y < 0){
		qx = -p3mu.x();
		qy = -p3mu.y();
		t = -2.*Qmax*(p3mu.t()-p3mu.z());
	}
	
	double mg2 = t_channel_mD2->get_mD2(T) / 2.;
    double xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1. - xbar;
	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - qx;
	MA[1] = kmu.y() - qy;
	MB[0] = kmu.x() - xbar*qx;
	MB[1] = kmu.y() - xbar*qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1.-xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1.-xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1.-xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		

    double M2_elastic = M2_Qq2Qq(t, params);
    double Pg = alpha_s(kt2, T) * one_minus_xbar
			* (1. + std::pow(one_minus_xbar, 2)) / 2.
 			* (CF*A2 + CF*B2 - (2.*CF-CA)*AB)/CA;

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
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>=.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
	double t = -2.*Qmax*(p4mu.t()+p4mu.z());
	
	if(y < 0){
		qx = -p3mu.x();
		qy = -p3mu.y();
		t = -2.*Qmax*(p3mu.t()-p3mu.z());
	}
	
	double mg2 = t_channel_mD2->get_mD2(T) / 2.;
    double xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	if (y<0 && xbar>0.5) return 0.;
	double one_minus_xbar = 1. - xbar;
	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - qx;
	MA[1] = kmu.y() - qy;
	MB[0] = kmu.x() - xbar*qx;
	MB[1] = kmu.y() - xbar*qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1.-xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1.-xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1.-xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		

    double M2_elastic = M2_Qg2Qg(t, params);
	double Pg = 1.;
	if (y>0){
		Pg = alpha_s(kt2, T) * one_minus_xbar
			* (1. + std::pow(one_minus_xbar, 2)) / 2.
 			* (CF*A2 + CF*B2 - (2.*CF-CA)*AB)/CA;
	} else {
		Pg = alpha_s(kt2, T) * xbar * one_minus_xbar
			* (1. + std::pow(xbar, 4) + std::pow(one_minus_xbar, 4)) / 2.
 			* (A2 + B2 - AB);
	}

	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}


//=============Basic for 3->2===========================================
// Not CoM frame of  1,2,k, but CoM frame of 1 and 2
// sampled final states within 3+4 CoM frame
double M2_Qqg2Qq(const double * x_, void * params_){
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
	double M2_elastic = M2_Qq2Qq(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

double M2_Qgg2Qg(const double * x_, void * params_){
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
	double M2_elastic = M2_Qg2Qg(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

////////////////// Hard Gluon ////////////////////////////////
/// g + q --> g + q + g
double M2_gq2gqg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
	double t = -2.*Qmax*(p4mu.t()+p4mu.z());

	double mg2 = t_channel_mD2->get_mD2(T) / 2.;
	double xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
	double one_minus_xbar = 1.-xbar;

	if (y>0 && xbar>0.5) return 0.0;
	if(y < 0){
		qx = -p3mu.x();
		qy = -p3mu.y();
		t = -2.*Qmax*(p3mu.t()-p3mu.z());
	}
	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - xbar*qx;
	MA[1] = kmu.y() - xbar*qy;
	MB[0] = kmu.x() - qx;
	MB[1] = kmu.y() - qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1.-xbar+xbar*xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1.-xbar+xbar*xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1.-xbar+xbar*xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		
    double M2_elastic = M2_gq2gq(t, params);
	
	double Pg = 1.;
	if (y > 0.){
		Pg = alpha_s(kt2, T)
				* (1. + std::pow(xbar, 4) + std::pow(one_minus_xbar, 4)) / 2.
	 			* (A2 + B2 - AB);
	}else{
		Pg = alpha_s(kt2, T) * one_minus_xbar
				* (1. + std::pow(one_minus_xbar, 2)) / 2.
	 			* (CF*A2 + CF*B2 - (2.*CF-CA)*AB)/CA;
	}
	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}

/// g + g --> g + g + g
double M2_gg2ggg(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
    double t = -2.*Qmax*(p4mu.t()+p4mu.z());

	double mg2 = t_channel_mD2->get_mD2(T) / 2.;
    double xbar = (kmu.t()+std::abs(kmu.z()))/sqrts;
    double one_minus_xbar = 1.-xbar;

	if (y>0 && xbar>0.5) return 0.0;
	if (y<0 && xbar>0.5) return 0.0;

	if(y < 0){
		qx = -p3mu.x();
		qy = -p3mu.y();
		t = -2.*Qmax*(p3mu.t()-p3mu.z());
	}

	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - xbar*qx;
	MA[1] = kmu.y() - xbar*qy;
	MB[0] = kmu.x() - qx;
	MB[1] = kmu.y() - qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1.-xbar+xbar*xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1.-xbar+xbar*xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1.-xbar+xbar*xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		

    double M2_elastic = M2_gg2gg(t, params);
    double Pg = alpha_s(kt2, T)
			* (1. + std::pow(xbar, 4) + std::pow(one_minus_xbar, 4)) / 2.
 			* (A2 + B2 - AB );
	
	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}


//=============Basic for 3->2===========================================
// Not CoM frame of  1,2,k, but CoM frame of 1 and 2
// sampled final states within 3+4 CoM frame
double M2_gqg2gq(const double * x_, void * params_){
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
	double M2_elastic = M2_gq2gq(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

double M2_ggg2gg(const double * x_, void * params_){
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
	double M2_elastic = M2_gg2gg(t, params);
	double x2M2 = x*x*M2;

	// 1->2
	double iD1 = 1./(kt2 + (1.-xbar)*mD2/2.);
	double iD2 = 1./(kt2 + qt2 + 2.*kt*qmu.x() + (1.-xbar)*mD2/2.);
	double Pg = 48.*M_PI*alpha_s(kt2, T)*std::pow(1.-xbar, 2)*(
			kt2*std::pow(iD1-iD2, 2) + qt2*std::pow(iD2, 2) - 2.*kt*qmu.x()*(iD1-iD2)*iD2
		);

	double detail_balance_factor = 1./16.;
	double Jacobian = 1./std::pow(2*M_PI, 2)/8.*(1. - M2/s);
	// 2->3 = 2->2 * 1->2
	return M2_elastic * Pg * Jacobian * detail_balance_factor;
}

/////////// Gluon to q qbar /////////////////

/// g + g --> q + g + qbar
double M2_gg2qgqbar(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
    double t = -2.*Qmax*(p4mu.t()+p4mu.z());

	double mg2 = t_channel_mD2->get_mD2(T)/2.;
    double xbar = 2*kmu.t()/sqrts;
    double one_minus_xbar = 1.-xbar;

	if(y < 0){
		return 0.;
	}

	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - xbar*qx;
	MA[1] = kmu.y() - xbar*qy;
	MB[0] = kmu.x() - qx;
	MB[1] = kmu.y() - qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		

    double M2_elastic = M2_gg2gg(t, params);
    double Pg = alpha_s(kt2, T) 
			* one_minus_xbar * xbar
			* (std::pow(xbar,2) + std::pow(one_minus_xbar, 2))/ 2.
 			* ( CF*A2 + CF*B2 - (2*CF-CA)*AB ) / 2.;

	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}

/// g + q --> q + q + qbar
double M2_gq2qqqbar(const double * x_, void * params_){
	// unpack variables, parameters and check integration range
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T = params[1];
	double M2 = params[2]*params[2];
	double Qmax = (s-M2)/2./sqrts;
	double log_1_ktT = x_[0], y_norm = x_[1], // in 3-body com frame---(1)
		   costheta34 = x_[2], phi34 =x_[3]; // in p3-p4 com frame---(2)
		   // (1) and (2) share the same z-direction, and are related by a boost
	// check bounds
	double kt = T*(std::exp(log_1_ktT)-1.);
	if (std::abs(costheta34)>=1.||phi34<=0.||phi34>=2.*M_PI||kt>=Qmax||kt<=0.||std::abs(y_norm)>.999)
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
	double kt2 = kmu.x()*kmu.x() + kmu.y()*kmu.y();
    double t = -2.*Qmax*(p4mu.t()+p4mu.z());

	
	double mg2 = t_channel_mD2->get_mD2(T)/2.;
    double xbar = 2*kmu.t()/sqrts;
    double one_minus_xbar = 1.-xbar;
	if (y<0) return 0.;
   	double MA[2], MB[2], MC[2], A[2], B[2];
	MA[0] = kmu.x() - xbar*qx;
	MA[1] = kmu.y() - xbar*qy;
	MB[0] = kmu.x() - qx;
	MB[1] = kmu.y() - qy;
	MC[0] = kmu.x();
	MC[1] = kmu.y();
	double DA = MA[0]*MA[0] + MA[1]*MA[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;
	double DB = MB[0]*MB[0] + MB[1]*MB[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;
	double DC = MC[0]*MC[0] + MC[1]*MC[1] + (1-CA/CF*xbar+CA/CF*xbar*xbar)*mg2;

	A[0] = MA[0]/DA-MB[0]/DB;
	A[1] = MA[1]/DA-MB[1]/DB;
	B[0] = MA[0]/DA-MC[0]/DC;
	B[1] = MA[1]/DA-MC[1]/DC;
	
	double A2 = std::pow(A[0],2) + std::pow(A[1],2);
	double B2 = std::pow(B[0],2) + std::pow(B[1],2);
	double AB = A[0]*B[0] + A[1]*B[1];
		

    double M2_elastic = M2_gq2gq(t, params);
    double Pg = alpha_s(kt2, T) 
			* one_minus_xbar * xbar
			* (std::pow(xbar,2) + std::pow(one_minus_xbar, 2))/ 2.
 			* ( CF*A2 + CF*B2 - (2*CF-CA)*AB ) / 2.;

	// Jacobian
	double J = (1.0 - M2/s34) * M_PI/8./std::pow(2*M_PI,5) * kt * (kt + T) * ymax;
	// 2->3 = 2->2 * 1->2
	return c48pi*M2_elastic*Pg*J;
}

