#ifndef WORKFLOW_H
#define WORKFLOW_H

#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>

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
typedef boost::variant<Rate22, Rate23> Process;
extern std::vector<Process> MyProcesses;

void initialize(std::string);
int update_particle_momentum(double dt, double temp, std::vector<double> v3cell, 
				double D_formation_t, fourvec incoming_p, std::vector<fourvec> & FS);

void probe_test(double E0, double T, double dt, int Nsteps, 
				int Nparticles, std::string mode);		
#endif
