#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>
#include <map>

struct particle{
	// mass, x, p, t, all in units of [GeV^a]
	int pid;
	double mass, weight;
	bool is_vac, is_virtual;
	
	double T0, mfp0, Tf; // production temperature, local mfp
	fourvec x0; // production location
	fourvec x; // current location
	fourvec p0; // production momentum
	fourvec p; // current momentum

	std::vector<particle> radlist;
	std::vector<double> vcell;

	void freestream(double dt){
		double a = dt/p.t();
		x.a[0] = x.t() + dt;
		x.a[1] = x.x() + p.x()*a;
		x.a[2] = x.y() + p.y()*a;
		x.a[3] = x.z() + p.z()*a;
	}
};

typedef Rate<LO, 2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<GB, 2, 2, double(*)(const double*, void*)> Rate23;
typedef Rate<GB, 2, 4, double(*)(const double*, void*)> Rate32;
typedef EffRate12<2, double(*)(const double*, void*)> Rate12;
typedef EffRate21<2, double(*)(const double*, void*)> Rate21;
typedef boost::variant<Rate22, Rate23, Rate32, Rate12, Rate21> Process;
extern std::map<int, std::vector<Process>> AllProcesses;

void init_process(Process& r, std::string mode, std::string table_path);
void initialize(std::string mode, std::string setting_path, std::string table_path, double mu, double afix, 
		double K, double a, double b, double p, double q, double gamma, double cut, double Rvac);

int update_particle_momentum_Lido(double dt, double temp, std::vector<double> v3cell, particle & pIn, std::vector<particle> & pOut_list);

double formation_time(fourvec p, fourvec k, double T, int split);
double calcualte_dt_from_dtau(fourvec x, fourvec p, double tau, double dtau);
void output(const std::vector<particle> plist, std::string fname);
void output_oscar(const std::vector<particle> plist, int abspid, std::string fname);
double mean_pT(const std::vector<particle> plist);
double mean_E(const std::vector<particle> plist);
#endif
