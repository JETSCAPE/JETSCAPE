#include "Xsection.h"
#include "matrix_elements.h"
#include "integrator.h"
#include "minimizer.h"
#include "sampler.h"
#include "predefine.h"
#include <fstream>
#include "random.h"
#include "approx_functions.h"

template<>
Xsection<2, double(*)(const double, void*)>::
	Xsection(std::string Name, std::string configfile,
			double(*f)(const double, void*)):
StochasticBase<2>(Name+"/xsection", configfile),
_f(f), fast_exp_(0., 15., 1000)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);

	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string model_name = strs[0];
	std::string process_name = strs[1];
	auto tree = config.get_child(model_name+"."+process_name);
	_mass = tree.get<double>("mass");

	// Set Approximate function for X and dX_max
	StochasticBase<2>::_ZeroMoment->SetApproximateFunction(approx_X22);
	StochasticBase<2>::_FunctionMax->SetApproximateFunction(approx_dX22_max);
}

template<>
Xsection<3, double(*)(const double *, void*)>::
	Xsection(std::string Name, std::string configfile,
			double(*f)(const double*, void*)):
StochasticBase<3>(Name+"/xsection", configfile),
_f(f), fast_exp_(0., 15., 1000)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);

	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string model_name = strs[0];
	std::string process_name = strs[1];
	auto tree = config.get_child(model_name+"."+process_name);
	_mass = tree.get<double>("mass");

	// Set Approximate function for X and dX_max
	StochasticBase<3>::_ZeroMoment->SetApproximateFunction(approx_X23);
	StochasticBase<3>::_FunctionMax->SetApproximateFunction(approx_dX23_max);
}

template<>
Xsection<4, double(*)(const double *, void*)>::
	Xsection(std::string Name, std::string configfile,
			double(*f)(const double*, void*)):
StochasticBase<4>(Name+"/xsection", configfile),
_f(f), fast_exp_(0., 15., 1000)
{
	// read configfile
	boost::property_tree::ptree config;
	std::ifstream input(configfile);
	read_xml(input, config);

	std::vector<std::string> strs;
	boost::split(strs, Name, boost::is_any_of("/"));
	std::string model_name = strs[0];
	std::string process_name = strs[1];
	auto tree = config.get_child(model_name+"."+process_name);
	_mass = tree.get<double>("mass");

	// Set Approximate function for X and dX_max
	//StochasticBase<4>::_ZeroMoment->SetApproximateFunction(approx_X32);
	//StochasticBase<4>::_FunctionMax->SetApproximateFunction(approx_dX32_max);
}

/*****************************************************************/
/*************************sample dX/dPS **************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
void Xsection<2, double(*)(const double, void*)>::
		sample(std::vector<double> parameters,
				std::vector< fourvec > & FS){
	double sqrts = std::max(parameters[0], 1.4), temp = parameters[1];
	double s = std::pow(sqrts,2);
    // transform w = -log(1-t/Temp^2)
    // since nearly everything happens when t ~ T^2
	auto dXdw = [s, temp, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return this->_f(t, params)*Jacobian;
	};
	double tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0.;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp);
	double w = sample_1d(dXdw, {wmin, wmax}, StochasticBase<2>::GetFmax(parameters).s);
	double t = temp*temp*(1.-fast_exp_(-w));
	// sample phi
	double phi = Srandom::dist_phi(Srandom::gen);
	double cosphi = std::cos(phi), sinphi = std::sin(phi);
	double E = (s+_mass*_mass)/2./sqrts; // EQ = (s+M^2)/2sqrts
	double p = (s-_mass*_mass)/2./sqrts; // pQ = (s-M^2)/2sqrts
	double costheta = 1. + t/(2.*p*p); // deflection angle
	double sintheta = std::sqrt(1.-costheta*costheta);
	FS.clear();
	FS.resize(2); // HQ + light parton
	FS[0] = {E, p*sintheta*cosphi, p*sintheta*sinphi, p*costheta};
	FS[1] = {p, -p*sintheta*cosphi, -p*sintheta*sinphi, -p*costheta};
}

/*------------------Implementation for 2 -> 3--------------------*/
template<>
void Xsection<3, double(*)(const double*, void*)>::
	sample(std::vector<double> parameters,
			std::vector< fourvec > & FS){
	double sqrts = std::max(parameters[0], 1.4), temp = parameters[1],
		   delta_t = parameters[2];
	double s = sqrts*sqrts;
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	// x0 = log(1+kt/T), x1 = y/ymax, x2 = -log(1+(1-costheta34)/(T/Qmax)^2), c3 = phi34
	double x2min = -std::log(1.+2./std::pow(temp/Qmax, 2)),
		   x2max = 0.;
	auto dXdPS = [s, temp, delta_t, Qmax, this](double * PS){
		double x2 = PS[2];
		double x[4] = {PS[0], PS[1],
				1.0 - (std::exp(-x2)-1.)*std::pow(temp/Qmax, 2), PS[3]};
		double Jacobian = std::exp(-x2)*std::pow(temp/Qmax, 2);
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(x, params)/2./(s-_mass*_mass)*Jacobian;
	};

	double xmin[4] = {0., -1.,  x2min, 0.};
	double xmax[4] = {umax, 1., x2max, 2.*M_PI};
	double fmax = StochasticBase<3>::GetFmax(parameters).s;
	bool status = true;
	auto res = sample_nd(dXdPS, 4, {{xmin[0], xmax[0]}, {xmin[1], xmax[1]},
									{xmin[2], xmax[2]}, {xmin[3], xmax[3]}},
									fmax, status);
	// deconvolve parameter
	double kt = temp*(std::exp(res[0])-1.);
	double yk = res[1]*std::acosh(Qmax/kt);
	double costheta34 = 1.0 - (std::exp(-res[2])-1.)*std::pow(temp/Qmax, 2);
	double phi34 = res[3];
	// sample phi
	double phi = Srandom::dist_phi(Srandom::gen);
	// reconstruct momentums
	double M2 = _mass*_mass;
	fourvec kmu{kt*std::cosh(yk), kt*std::cos(phi),
			    kt*std::sin(phi), kt*std::sinh(yk)}; // k
	double s34 = s - 2.*sqrts*kmu.t();
	double sqrts34 = std::sqrt(s34);
	double Q34 = (s34-M2)/2./sqrts34, E34 = (s34+M2)/2./sqrts34;
	// construct p3, p4 in (2) frame
	double cosphi34 = std::cos(phi34), sinphi34 = std::sin(phi34),
		   sintheta34 = std::sqrt(1.-std::pow(costheta34,2));
	fourvec p3mu{E34, Q34*sintheta34*cosphi34, Q34*sintheta34*sinphi34, Q34*costheta34},
		p4mu{Q34, -Q34*sintheta34*cosphi34, -Q34*sintheta34*sinphi34, -Q34*costheta34};

	// boost p3, p4 back to 3-body CoM Frame
	double V0 = sqrts - kmu.t();
	double v34[3] = {-kmu.x()/V0, -kmu.y()/V0, -kmu.z()/V0};
	p3mu = p3mu.boost_back(v34[0], v34[1], v34[2]); // p3
	p4mu = p4mu.boost_back(v34[0], v34[1], v34[2]); // p4

	//FS.resize(3); // HQ + light parton
	FS.clear();
	FS.push_back(p3mu);
	FS.push_back(p4mu);
	FS.push_back(kmu);
}
/*------------------Implementation for 3 -> 2--------------------*/
template<>
void Xsection<4, double(*)(const double*, void*)>::
	sample(std::vector<double> parameters,
			std::vector< fourvec > & FS){
	double sqrts = parameters[0], temp = parameters[1],
		   xinel = parameters[2], yinel = parameters[3];
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, xinel, yinel, this](const double * PS){
		double M = this->_mass;
		double params[5] = {s, temp, M, xinel, yinel};
		return this->_f(PS, params);
	};
	double fmax = std::exp(StochasticBase<4>::GetFmax(parameters).s);
	//LOG_INFO << "dX(sqrts, T, x, y, dt) " << sqrts << " " << temp << " " << xinel << " " << yinel << " " << dt;
	bool status = true;
	auto res = sample_nd(dXdPS, 2, {{-1., 1.}, {0., 2.*M_PI}}, fmax, status);
	/*if (status == false){
		FS.resize(2);
		double s12 = xinel*(s-_mass*_mass) + _mass*_mass;
		double sqrts12 = std::sqrt(s12);
		fourvec a{(s12+_mass*_mass)/2./sqrts12, 0, 0, (s12-_mass*_mass)/2./sqrts12};
		fourvec b{(s12-_mass*_mass)/2./sqrts12, 0, 0, -(s12-_mass*_mass)/2./sqrts12};
		FS[0] = a;
		FS[1] = b;
		return;
	}*/
	double costheta34 = res[0], phi34 = res[1];
	double sintheta34 = std::sqrt(1. - costheta34*costheta34),
			sinphi34 = std::sin(phi34), cosphi34 = std::cos(phi34);
	double M2 = _mass*_mass;
	double s12 = xinel*(s-M2) + M2; // 0 < params[3] = xinel = (s12-M2)/(s-M2) < 1
	double s1k = (yinel*(1-s12/s)*(1-M2/s12)+M2/s12)*s; // 0 < params[4] = yinel = (s1k/s-M2/s12)/(1-s12/s)/(1-M2/s12) < 1
	double sqrts12 = std::sqrt(s12);
	double E1 = (s12+M2)/2./sqrts12,  p1 = (s12-M2)/2./sqrts12;
	double E3 = (s+M2)/2./sqrts, p3 = (s-M2)/2./sqrts;
	double k = (s-s12)/2./sqrts12;
	double costhetak = (M2 + 2.*E1*k - s1k)/2./p1/k;
	double sinthetak = std::sqrt(1. - costhetak*costhetak);
	double kt = k*sinthetak, kz = k*costhetak;
	// get final state
	fourvec Ptot{sqrts12+k, kt, 0., kz};
	double vcom[3] = {Ptot.x()/Ptot.t(), Ptot.y()/Ptot.t(),Ptot.z()/Ptot.t()};
	// final state in 34-com frame
	fourvec p3mu{E3, p3*sintheta34*cosphi34, p3*sintheta34*sinphi34, p3*costheta34};
	fourvec p4mu{p3, -p3mu.x(), -p3mu.y(), -p3mu.z()};
	// boost final state back to 12-com frame
	p3mu = p3mu.boost_back(vcom[0], vcom[1], vcom[2]);
	p4mu = p4mu.boost_back(vcom[0], vcom[1], vcom[2]);
	FS.resize(2);
	FS[0] = p3mu;
	FS[1] = p4mu;
}
/*****************************************************************/
/*******************find max of dX/dPS ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<2, double(*)(const double, void*)>::
	find_max(std::vector<double> parameters){
    double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2)
    // since nearly everything happens when t ~ T^2
	auto minus_dXdw = [s, temp, this](const double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return -(this->_f(t, params)*Jacobian);
	};
	double tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp);
	double res = -minimize_1d(minus_dXdw, {wmin, wmax}, 1e-8, 100, 1000);
	return scalar{res*1.5};
}
/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<3, double(*)(const double*, void*)>::
	find_max(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1],
		   delta_t = parameters[2];

	double s = sqrts*sqrts;
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	// x0 = log(1+kt/T), x1 = yk, x2 = -log(1+(1-costheta34)/(T/Qmax)^2), c3 = phi34
	double x2min = -std::log(1.+2./std::pow(temp/Qmax, 2)),
		   x2max = 0.;
	double Lx2 = x2max - x2min;
	auto dXdPS = [s, temp, delta_t, Qmax, x2min, x2max, this](double * PS){
		double x2 = PS[2];
		double x[4] = {PS[0], PS[1],
					1.0 - (std::exp(-x2)-1.)*std::pow(temp/Qmax, 2), PS[3]};
		double Jacobian = std::exp(-x2)*std::pow(temp/Qmax, 2);
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(x, params)/2./(s-_mass*_mass)*Jacobian;
	};
	auto nega_dXdPS = [s, temp, delta_t, Qmax, x2min, x2max, this](double * PS){
		double x2 = PS[2];
		double x[4] = {PS[0], PS[1],
					1.0 - (std::exp(-x2)-1.)*std::pow(temp/Qmax, 2), PS[3]};
		double Jacobian = std::exp(-x2)*std::pow(temp/Qmax, 2);
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return -this->_f(x, params)/2./(s-_mass*_mass)*Jacobian;
	};
	// use MC_maximize to get into the vincinity ot the extrma
	auto startloc = MC_maximize(dXdPS, 4,
			{{umax*0.4,umax*0.6}, {-0.2, 0.2},
			 {x2min+Lx2/3., x2max-Lx2/3.}, {-0.5, 0.5}}, 400);
	// use the best result of MC_maximize and determine the step of the simplex minimization method
	std::vector<double> step = {umax/20., 0.1, (x2max-x2min)/20., 0.1};
	double L[4] = {0, -1, x2min, 0};
	double H[4] = {umax, 1, x2max, 2*M_PI};
	for(int i=0; i<4; i++){
		double dx = std::min(H[i]-startloc[i], startloc[i]-L[i])/2.;
		step[i] = std::min(dx, step[i]);
	}
	// find the more precise maximum by the simplex method
    double val = -minimize_nd(nega_dXdPS, 4, startloc, step,
										4000, Lx2*umax*4*M_PI/1e12);
	// save max*1.5 just to be safe
	return scalar{val*1.5};
}
/*------------------Implementation for 3 -> 2--------------------*/
template<>
scalar Xsection<4, double(*)(const double*, void*)>::
	find_max(std::vector<double> parameters){
		double sqrts = parameters[0], temp = parameters[1],
			   xinel = parameters[2], yinel = parameters[3];
		double s = sqrts*sqrts;
		auto dXdPS = [s, temp, xinel, yinel, this](const double * PS){
			double M = this->_mass;
			double params[5] = {s, temp, M, xinel, yinel};
			return this->_f(PS, params);
		};
		auto nega_dXdPS = [s, temp, xinel, yinel, this](const double * PS){
			double M = this->_mass;
			double params[5] = {s, temp, M, xinel, yinel};
			return -this->_f(PS, params);
		};
		// use MC_maximize to get into the vincinity ot the extrma
		auto startloc = MC_maximize(dXdPS, 2, {{-1., 1.}, {0.,2.*M_PI}}, 100);
		// use the best result of MC_maximize and determine the step of the simplex minimization method
		std::vector<double> step = {0.1, 0.1};
		double L[2] = {-1., 0.};
		double H[2] = {1., 2.*M_PI};
		for(int i=0; i<2; i++){
			double dx = std::min(H[i]-startloc[i], startloc[i]-L[i])/2.;
			step[i] = std::min(dx, step[i]);
		}
		// find the more precise maximum by the simplex method
		double val = std::log(-minimize_nd(nega_dXdPS, 2, startloc, step, 1000, 2*2*M_PI/1e8)*2.);
		// save max*1.5 just to be safe

		return scalar{val};
}
/*****************************************************************/
/*************************Integrate dX ***************************/
/*****************************************************************/
/*------------------Implementation for 2 -> 2--------------------*/
template<>
scalar Xsection<2, double(*)(const double, void*)>::
		calculate_scalar(std::vector<double> parameters){
	double s = std::pow(parameters[0],2), temp = parameters[1];
    // transform w = -log(1-t/Temp^2)
    // since nearly everything happens when t ~ T^2
	auto dXdw = [s, temp, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
		return this->_f(t, params)*Jacobian;
	};
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0.;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double res = quad_1d(dXdw, {wmin, wmax}, error);
	return scalar{res};
}
/*------------------Implementation for 2 -> 3--------------------*/
template<>
scalar Xsection<3, double(*)(const double*, void*)>::
				calculate_scalar(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1],
		   delta_t = parameters[2];
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, delta_t, this](const double * PS){
		double M = this->_mass;
		double params[4] = {s, temp, M, delta_t};
		return this->_f(PS, params)/2./(s-_mass*_mass);
	};
	double Qmax = (s-_mass*_mass)/2./sqrts;
	double umax = std::log(1.+Qmax/temp);
	double xmin[4] = {0., -1., -1., 0.};
	double xmax[4] = {umax, 1., 1., 2.*M_PI};
	double error;
	double res = vegas(dXdPS, 4, xmin, xmax, error);
	return scalar{res};
}
/*------------------Implementation for 3 -> 2--------------------*/
template<>
scalar Xsection<4, double(*)(const double*, void*)>::
				calculate_scalar(std::vector<double> parameters){

	// xiel = s12/s
	double sqrts = parameters[0], temp = parameters[1],
		   xinel = parameters[2], yinel = parameters[3];
	double s = sqrts*sqrts;
	auto dXdPS = [s, temp, xinel, yinel, this](const double * PS){
		double M = this->_mass;
		double params[5] = {s, temp, M, xinel, yinel};
		std::vector<double> res{this->_f(PS, params)};
		return res;
	};
	double xmin[2] = {-1., 0.};
	double xmax[2] = {1., 2.*M_PI};
	double error;
	auto res = quad_nd(dXdPS, 2, 1, xmin, xmax, error);
	double result = std::log(res[0]);

	return scalar{result};
}
/*****************************************************************/
/**************Integrate dX \Delta p^mu***************************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<size_t N, typename F>
fourvec Xsection<N, F>::calculate_fourvec(std::vector<double> parameters){
	return fourvec::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
fourvec Xsection<2, double(*)(const double, void*)>::
		calculate_fourvec(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz*X>
	auto dpz_dXdw = [s, temp, p0, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
        double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz*Jacobian;
	};
	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
	double dpz = quad_1d(dpz_dXdw, {wmin, wmax}, error);
	return fourvec{0., 0., 0., dpz};
}

/*****************************************************************/
/**************Integrate dX \Delta p^mu*p^nu *********************/
/*****************************************************************/
/*------------------Default Implementation-----------------------*/
template<size_t N, typename F>
tensor Xsection<N, F>::calculate_tensor(std::vector<double> parameters){
	return tensor::unity();
}
/*------------------Implementation for 2 -> 2--------------------*/
template<>
tensor Xsection<2, double(*)(const double, void*)>::
	calculate_tensor(std::vector<double> parameters){
	double sqrts = parameters[0], temp = parameters[1];
	double s = std::pow(sqrts,2);
	const double p0 = (s-_mass*_mass)/2./sqrts;
	// <dpz^2*X>
	auto dpzdpz_dXdw = [s, temp, p0, this](double w) {
        double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-fast_exp_(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
        double dpz = p0*(cos_theta13-1.);
		return this->_f(t, params)*dpz*dpz*Jacobian;
	};
	// <dpt^2*X>
	auto dptdpt_dXdw = [s, temp, p0, this](double w) {
		double T2 = temp*temp;
		double params[3] = {s, temp, this->_mass};
		double t = T2*(1.-std::exp(-w));
		double Jacobian = T2 - t;
        double cos_theta13 = 1. + t/(2*p0*p0);
		double dptdpt = p0*p0*(1.-cos_theta13*cos_theta13);
		return this->_f(t, params)*dptdpt*Jacobian;
	};

	double error, tmin = -std::pow(s-_mass*_mass, 2)/s, tmax=0;
    double wmin = -std::log(1.-tmin/temp/temp),
		   wmax = -std::log(1.-tmax/temp/temp+1e-9);
    double dpzdpz = quad_1d(dpzdpz_dXdw, {wmin, wmax}, error);
	double dptdpt = quad_1d(dptdpt_dXdw, {wmin, wmax}, error);
	return tensor{0., 	0., 		0., 		0.,
				  0., 	dptdpt/2., 	0., 		0.,
				  0., 	0., 		dptdpt/2.,	0.,
				  0., 	0., 		0., 		dpzdpz};
}
// instance:
template class Xsection<2, double(*)(const double, void*)>;
template class Xsection<3, double(*)(const double*, void*)>;
template class Xsection<4, double(*)(const double*, void*)>;
