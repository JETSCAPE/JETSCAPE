#include "approx_functions.h"

// Xsection
scalar approx_X22(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	return scalar{1.0/T/T};
}

scalar approx_dX22_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	return scalar{1.0/std::pow(T, 2)};
}

scalar approx_X23(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{1./(1.+1./std::pow(3*delta_t*T,2))
					/std::pow(T, 2)*sqrts*sqrts};
}

scalar approx_dX23_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{1./(1.+1./std::pow(3*delta_t*T,2))
					/std::pow(T, 4)*sqrts*sqrts};
}
// Rate
scalar approx_R22(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	return scalar{T};
}

scalar approx_dR22_max(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	return scalar{T};
}

scalar approx_R23(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{T};
}

scalar approx_dR23_max(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{1./(1.+1./std::pow(delta_t*T,2))};
}

