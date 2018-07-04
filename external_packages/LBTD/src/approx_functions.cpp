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
	double a = std::pow(delta_t*T,2);
	return scalar{a/(1.+a)/std::pow(T,3)};
}

scalar approx_dX23_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[2];
	double a = std::pow(delta_t*T,2);
	return scalar{a/(2.+a)/std::pow(T,4)};
}

scalar approx_X32(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[4];
	double a = std::pow(delta_t*T,2);
	return scalar{a/(1+a)/std::pow(T, 2)/sqrts};
}

scalar approx_dX32_max(std::vector<double> params){
	double sqrts = params[0];
	double T = params[1];
	double delta_t = params[4];
	double a = delta_t*T;
	return scalar{a/(.2+a)/std::pow(T, 2)};
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
	double a = std::pow(delta_t*T, 2);
	return scalar{a/(a+5.)*T};
}

scalar approx_dR23_max(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	double delta_t = params[2];
	double a = std::pow(delta_t*T, 2);
	return scalar{a/(a+3.)*T};
}

scalar approx_R32(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	double delta_t = params[2];
	return scalar{delta_t*std::pow(T,3)/E};
}

scalar approx_dR32_max(std::vector<double> params){
	double E = params[0];
	double T = params[1];
	double delta_t = params[2];
	double a = std::pow(delta_t*T, 2);
	return scalar{a/(.25*a+1)/E};
}
