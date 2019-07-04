#include "predefine.h"
#include <cmath>
#include <boost/math/tools/roots.hpp>
#include "simpleLogger.h"


char LO[] = "Leading-order";
char GB[] = "Gunion-Bertsch";
char GBHT[] = "Gunion-Bertsch-higher-twist";

bool type1_warned = false;
bool type2_warned = false; 
bool type3_warned = false; 

const double c4d9 = 4./9.;
const double c1d9 = 1./9.;
const double c16pi = 16.*M_PI;
const double c48pi = 48.*M_PI;
const double c16pi2 = 16.*M_PI*M_PI;
const double c64d9pi2 = 64./9.*M_PI*M_PI;
const double c256pi4 = 256.*std::pow(M_PI, 4);
const double fmc_to_GeV_m1 = 5.026;

// number of color=3 (3*3-1 = 8 gluons), number of flavor=3, (u,d,s quark)
const int Nc = 3, nf = 3;
const double CF = 4./3.;
const double CA = 3.;
const double CF_over_CA = CF/CA;
const double TR = 0.5;

// the prefractor for gluon debye mass with Boltzmann statistics
// mD^2 = 8\pi*(Nc+nf)*alpha_s*T^2 ~ 15*alpha_s*T^2
// Note that if using quantum statistics, this will be
// mD^2 = 4\pi/3*(Nc+nf/2)*alpha_s*T^2 ~ 18*alpha_s*T^2
const double pf_g = 8./M_PI*(Nc + nf); // prefractor for gluon self energy^2

// For QCD coupling constant
// alpha_s = alpha_0 = 4\pi/(11Nc/3-2Nf/3)/log(Q^2/LambdaQCD^2)
const double alpha0 = 4.*M_PI/(11./3.*Nc - 2./3.*nf); // alpha_s(Q2 = e*Lambda2)
const double alpha_max = 1.0;
const double Lambda = 0.2; // [GeV] Lambda QCD = 0.2 GeV
const double Lambda2 = Lambda*Lambda; // [GeV^2] Lambda QCD squared
const double mu2_NP = Lambda2*std::exp(alpha0/alpha_max); // minimum cut on Q2, where alpha = alpha_0
const double Tc = 0.154;
double cut;
double scale;
double afix;
double Rvac;
Debye_mass * t_channel_mD2;
const double LPM_prefactor = 0.75; // to match analytic calculation, 0.75 by default

qhat_params_struct qhat_params;

void initialize_mD_and_scale(int _mD_type, double _scale, double _afix, double _cut, double _Rvac){
        cut = _cut;
        scale = _scale;
        afix = _afix;
	Rvac = _Rvac;
        t_channel_mD2 = new Debye_mass(_mD_type);
}

void initialize_transport_coeff(double _K, double _a, double _b, double _p, double _q, double _gamma){
	qhat_params.K = _K;
	qhat_params.a = _a;
	qhat_params.b = _b;
	qhat_params.p = _p;
	qhat_params.q = _q;
	qhat_params.gamma = _gamma;
}

void echo(void){
	std::cout << "mu\tafix\tK\ta\tb\tp\tq\tgamma\tQcut\tRvac\n";
	std::cout <<scale<<"\t"<<afix<<"\t"<<qhat_params.K<<"\t"<<qhat_params.a
               <<"\t"<<qhat_params.b<<"\t"<<qhat_params.p<<"\t"<<qhat_params.q
               <<"\t"<<qhat_params.gamma<<"\t"<<cut<<"\t"<<Rvac<<"\n";
}

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
                        size_t maxiter=200;
                        boost::math::tools::eps_tolerance<double> tol{
                         (std::numeric_limits<double>::digits * 3) / 4};
                        try{
                                auto result = boost::math::tools::toms748_solve(
                                        [&T](double x) {return pf_g*alpha_s(-x, T)*T*T - x;},
                                        0.001, 100., tol, maxiter);
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


// splitting function
double P_q2qg(double x){
	return CF*(1. + std::pow(1.-x, 2))/x;
}
double P_q2gq(double x){
	return CF*(1. + x*x)/(1.-x);
}
double P_g2gg(double x){
	return CA*(1.+std::pow(x, 4)+std::pow(1.-x,4))/x/(1.-x);
}
double P_g2qq(double x){
	return nf*(x*x + std::pow(1.-x, 2));
}


//=============running coupling=================================================
double alpha_s(double Q2, double T){
        if (afix > 0) {
                return afix;
        }
        else{
                double screen_scale2 = std::pow(scale*M_PI*T, 2);
                double mu2;
                if (Q2 < 0.) mu2 = std::max(-Q2, screen_scale2);
                else mu2 = std::max(Q2, screen_scale2);

                if (mu2 <= mu2_NP) return 1.0;
                else return alpha0 / std::log(mu2/Lambda2);
        }
}

bool is_file_exist(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

