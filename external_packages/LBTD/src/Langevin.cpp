#include "Langevin.h"
#include "random.h"

double A=0., B=0.;
double const tiny = 1e-10;

double kperp(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

double kpara(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

double dkpara_dp2(double E, double M, double T){
	return std::pow(T,3)*B*(-1.)/(2.*E*E*E*T);
}


void initialize_transport_coeff(double _A, double _B){
	A = _A; B = _B;
	std::cout << "A = " << A << ", B = " << B << std::endl;
};

void postpoint_update(	double dt, double M, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M + 1e-9);

	// step-1
	double kt = kperp(E0, M, T),
		   kl = kpara(E0, M, T);
	double drag = kl/(2.*E0*T) - std::pow((std::sqrt(kl)-std::sqrt(kt))/p0, 2);
		   
    double white_noise_holder[3];
    for (size_t i=0; i<3; ++i) 
		white_noise_holder[i] = Srandom::white_noise(Srandom::gen);

	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

    pOut.a[1] = Ct * white_noise_holder[0];
    pOut.a[2] = Ct * white_noise_holder[1];
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
					+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// step-2
	kt = kperp(pOut.t(), M, T);
	kl = kpara(pOut.t(), M, T);

	Ct = std::sqrt(kt*dt);
	Cl = std::sqrt(kl*dt);

    pOut.a[1] = Ct * white_noise_holder[0];
    pOut.a[2] = Ct * white_noise_holder[1];
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
					+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

void Ito_update(	double dt, double M, double T, std::vector<double> v, 
						const fourvec & pIn, fourvec & pOut){
	// Boost pIn to medium frame
	auto pIn_cell = pIn.boost_to(v[0], v[1], v[2]);
	// imaging rotating to a frame where pIn lies on z-axis
	double E0 = pIn_cell.t();
	double p0 = std::sqrt(E0*E0 - M*M+1e-9);

	double kt = kperp(E0, M, T),
		   kl = kpara(E0, M, T),
		   dkl_dp2 = dkpara_dp2(E0, M, T);
	double drag = kl/(2.*E0*T) - (kl - kt)/std::pow(p0, 2) - dkl_dp2;
		   
    double white_noise_holder[3];
    for (size_t i=0; i<3; ++i) 
		white_noise_holder[i] = Srandom::white_noise(Srandom::gen);

	double Ct = std::sqrt(kt*dt);
	double Cl = std::sqrt(kl*dt);

    pOut.a[1] = Ct * white_noise_holder[0];
    pOut.a[2] = Ct * white_noise_holder[1];
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) 
					+ std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

	// rotate back to the original frame
	pOut = pOut.rotate_back(pIn_cell);
	// boost back to lab frame
	pOut = pOut.boost_back(v[0], v[1], v[2]);
}

