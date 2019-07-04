#include "Langevin.h"
#include "random.h"
#include "predefine.h"
#include "matrix_elements.h"
double const tiny = 1e-10;

double delta_qhat(int pid, double E, double M, double T){
	double CR = (pid==21) ? CA : CF;
	if (pid == 21) E = std::sqrt(E*E + t_channel_mD2->get_mD2(T)/2.);
	double delta_qhat = CR/CF * qhat_params.K * std::pow(T, 3)
		/(1. + std::pow(qhat_params.a*(T+tiny)/Tc, qhat_params.p))
		/(1. + std::pow(qhat_params.b*(E+tiny)/(T+tiny), qhat_params.q));
	return delta_qhat;
}

double qhat_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double alphas_at_T = alpha_s(0, T); 
        double mD2 = t_channel_mD2->get_mD2(T);
        double Q2cut = cut*mD2;
        double thermal2 = std::pow(scale*M_PI*T, 2);
        double logs;
        if (Q2cut > thermal2 && afix < 0.){
            double log0 = std::log(1.+thermal2/mD2);
            double log1 = std::log(thermal2/Lambda2);
            double log2 = std::log(Q2cut/Lambda2);
            logs = log0 + log1*(1. - log1/log2);
        }
        else{
            logs = std::log(1.+Q2cut/mD2);
        }
	return alphas_at_T * CR * T * mD2 * ( logs);
}

double qhat_L_small_angle_LOpQCD(int pid, double E, double M, double T){
        double CR = (pid==21) ? CA : CF;
        double alphas_at_T = alpha_s(0, T);
        double mD2 = t_channel_mD2->get_mD2(T);
        double Q2cut = cut*mD2;
		double minf2 = .5*mD2;
        double thermal2 = std::pow(scale*M_PI*T, 2);
        double logs;
        if (Q2cut > thermal2 && afix < 0.){
            double log0 = std::log(1.+thermal2/minf2);
            double log1 = std::log(thermal2/Lambda2);
            double log2 = std::log(Q2cut/Lambda2);
            logs = log0 + log1*(1. - log1/log2);
        }
        else{
            logs = std::log(1.+Q2cut/minf2);
        }
        return alphas_at_T * CR * T * minf2 * (logs);
}


double qhat(int pid, double E, double M, double T){
	return  qhat_small_angle_LOpQCD(pid, E, M, T) 
	      + delta_qhat(pid, E, M, T);
}


double qhat_L(int pid, double E, double M, double T){
	double m0;
	double minf = std::sqrt(t_channel_mD2->get_mD2(T)/2.);
    m0 = std::max(minf, M);
    return  qhat_L_small_angle_LOpQCD(pid, E, M, T)
              + delta_qhat(pid, E, m0, T)/2. * std::pow(E/m0, qhat_params.gamma);                       
}

double dqhat_L_dp2(int pid, double E, double M, double T){
	double p2 = E*E - M*M + tiny;
	double dp2 = p2*.05;
	double Eprime = std::sqrt(E*E + dp2);
	return (qhat_L(pid, Eprime, M, T) - qhat_L(pid, E, M, T) ) /dp2;
}

void Ito_update(int pid, double dt_lab, double M, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
	// Boost dt to medium frame
	double dt = dt_lab*pIn_cell.t()/pIn.t();
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M+1e-9);

	double kt = qhat(pid, E0, M, T)/2.;
	double kl = qhat_L(pid, E0, M, T);
	double dkl_dp2 = dqhat_L_dp2(pid, E0, M, T);
	double drag = kl/(2.*E0*T) - (kl - kt)/std::pow(p0, 2) - dkl_dp2;
		   
	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

        pOut.a[1] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[2] = Ct * Srandom::white_noise(Srandom::gen);
        pOut.a[3] = p0 * (1. - drag * dt) + Cl * Srandom::white_noise(Srandom::gen);
        pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
		  	+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

