#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>
#include <map>

struct particle{
	int pid;
	bool freezeout;
	double mass;
	fourvec x;
	fourvec p;
	double t_rad, t_absorb;
	fourvec p0;
	std::vector<double> vcell;
	double Tf;
	void freestream(double dt){
		double a = dt/p.t();
		x.a[0] = x.t() + dt;
		x.a[1] = x.x() + p.x()*a;
		x.a[2] = x.y() + p.y()*a;
		x.a[3] = x.z() + p.z()*a;
	}
};

typedef Rate<2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<3, 3, double(*)(const double*, void*)> Rate23;
typedef Rate<3, 4, double(*)(const double*, void*)> Rate32;
typedef boost::variant<Rate22, Rate23, Rate32> Process;
extern std::map<int, std::vector<Process>> AllProcesses;
void initialize(std::string, std::string path, double mu);
int update_particle_momentum(double dt, double temp, std::vector<double> v3cell, int pid,
				double D_formation_t23, double D_formation_t32, fourvec incoming_p, std::vector<fourvec> & FS);

std::vector<double> probe_test(double E0, double T, double dt, int Nsteps,
				int Nparticles, std::string mode);
std::vector<std::vector<double>> rate_test(double E0, double T, double dt,
	int Nsteps, int Nparticles, std::string mode, double rescale);


#endif
