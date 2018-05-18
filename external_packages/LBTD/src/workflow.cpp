#include "workflow.h"
#include "simpleLogger.h"
#include <fstream>
#include "random.h"
#include "matrix_elements.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include "logo.h"

std::map<int, std::vector<Process> > AllProcesses;

void init_process(Process& r, std::string mode){
     switch(r.which()){
                        case 0:
                                if (boost::get<Rate22>(r).IsActive())
                                        if(mode == "new"){
                                                boost::get<Rate22>(r).initX("table.h5");
                                                boost::get<Rate22>(r).init("table.h5");
                                        } else{
                                                boost::get<Rate22>(r).loadX("table.h5");
                                                boost::get<Rate22>(r).load("table.h5");
                                        }
                                else return;
                                break;
                        case 1:
                                if (boost::get<Rate23>(r).IsActive())
                                        if(mode == "new"){
                                                boost::get<Rate23>(r).initX("table.h5");
                                                boost::get<Rate23>(r).init("table.h5");
                                        } else{
                                                boost::get<Rate23>(r).loadX("table.h5");
                                                boost::get<Rate23>(r).load("table.h5");
                                        }
                                else return;
                                break;
						case 2:
								if (boost::get<Rate32>(r).IsActive())
										if(mode == "new"){
												boost::get<Rate32>(r).initX("table.h5");
												boost::get<Rate32>(r).init("table.h5");
										} else{
												boost::get<Rate32>(r).loadX("table.h5");
												boost::get<Rate32>(r).load("table.h5");
										}
								else return;
								break;
                        default:
                                exit(-1);
                                break;
                }
}


void initialize(std::string mode, std::string path, double mu){
	print_logo();
	boost::property_tree::ptree config;
    initialize_mD_and_scale(1, mu);

	AllProcesses[4] = std::vector<Process>();
	AllProcesses[4].push_back( Rate22("Boltzmann/cq2cq", path, dX_Qq2Qq_dt) );
	AllProcesses[4].push_back( Rate22("Boltzmann/cg2cg", path, dX_Qg2Qg_dt) );
	AllProcesses[4].push_back( Rate23("Boltzmann/cq2cqg", path, M2_Qq2Qqg) );
	AllProcesses[4].push_back( Rate23("Boltzmann/cg2cgg", path, M2_Qg2Qgg) );
	AllProcesses[4].push_back( Rate32("Boltzmann/cqg2cq", path, Ker_Qqg2Qq) );
	AllProcesses[4].push_back( Rate32("Boltzmann/cgg2cg", path, Ker_Qgg2Qg) );

    AllProcesses[5] = std::vector<Process>();
    AllProcesses[5].push_back( Rate22("Boltzmann/bq2bq", path, dX_Qq2Qq_dt) );
    AllProcesses[5].push_back( Rate22("Boltzmann/bg2bg", path, dX_Qg2Qg_dt) );
    AllProcesses[5].push_back( Rate23("Boltzmann/bq2bqg", path, M2_Qq2Qqg) );
    AllProcesses[5].push_back( Rate23("Boltzmann/bg2bgg", path, M2_Qg2Qgg) );
	AllProcesses[5].push_back( Rate32("Boltzmann/bqg2bq", path, Ker_Qqg2Qq) );
	AllProcesses[5].push_back( Rate32("Boltzmann/bgg2bg", path, Ker_Qgg2Qg) );

	BOOST_FOREACH(Process& r, AllProcesses[4]) init_process(r, mode);
	BOOST_FOREACH(Process& r, AllProcesses[5]) init_process(r, mode);
}

int update_particle_momentum(double dt, double temp, std::vector<double> v3cell,
			int pid, double D_formation_t23, double D_formation_t32, fourvec incoming_p, std::vector<fourvec> & FS){
	int absid = std::abs(pid);
	auto p_cell = incoming_p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double D_formation_t23_cell = D_formation_t23 / incoming_p.t() * p_cell.t();
	double D_formation_t32_cell = D_formation_t32 / incoming_p.t() * p_cell.t();
	double dt_cell = dt / incoming_p.t() * p_cell.t();
	double E_cell = p_cell.t();
        std::vector<double> P_channels(AllProcesses[absid].size());
	double P_total = 0.;
	int channel = 0;
	double dR;
	BOOST_FOREACH(Process& r, AllProcesses[absid]){
		switch(r.which()){
			case 0:
				if (boost::get<Rate22>(r).IsActive())
					dR = boost::get<Rate22>(r).GetZeroM(
												{E_cell, temp}).s * dt_cell;
				else dR = 0.0;
				P_channels[channel] = P_total + dR;
				break;
			case 1:
				if (boost::get<Rate23>(r).IsActive())
					dR = boost::get<Rate23>(r).GetZeroM(
									{E_cell, temp, D_formation_t23_cell}).s * dt_cell;
				else dR = 0.0;
				P_channels[channel] = P_total + dR;
				break;
			case 2:
				if (boost::get<Rate32>(r).IsActive())
					dR = boost::get<Rate32>(r).GetZeroM(
									{E_cell, temp, D_formation_t32_cell}).s * dt_cell;
				else dR = 0.0;
				P_channels[channel] = P_total + dR;
				break;
			default:
				exit(-1);
				break;
		}
		P_total += dR;
		channel ++;
	}
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15) LOG_WARNING << "P_total = " << P_total << " may be too large";
	if ( Srandom::init_dis(Srandom::gen) > P_total) return -1;
	else{
		double p = Srandom::init_dis(Srandom::gen);
		for(int i=0; i<P_channels.size(); ++i){
			if (P_channels[i] > p) {
				channel = i;
				break;
			}
		}
	}
	// Do scattering
	switch(AllProcesses[absid][channel].which()){
		case 0:
			boost::get<Rate22>(AllProcesses[absid][channel]).sample({E_cell, temp}, FS);
			break;
		case 1:
			boost::get<Rate23>(AllProcesses[absid][channel]).sample(
											{E_cell, temp, D_formation_t23_cell}, FS);
			break;
		case 2:
			boost::get<Rate32>(AllProcesses[absid][channel]).sample(
											{E_cell, temp, D_formation_t32_cell}, FS);
			break;
		default:
			LOG_FATAL << "Channel = " << channel << " not exists";
			exit(-1);
			break;
	}
	// rotate it back and boost it back
	for(auto & pmu : FS) {
		pmu = pmu.rotate_back(p_cell);
		pmu = pmu.boost_back(v3cell[0], v3cell[1], v3cell[2]);
	}
	return channel;
}


std::vector<double> probe_test(double E0, double T, double dt=0.05, int Nsteps=100, int Nparticles=10000, std::string mode="old"){
	double fmc_to_GeV_m1 = 5.026;
	initialize(mode, "./settings.xml", 1.0);
	double M = 1.3;
	std::vector<double> dE;
	std::vector<particle> plist(Nparticles);
  fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)};
	for (auto & p : plist) {
		p.pid = 4;
		p.x = fourvec{0,0,0,0};
		p.p = p0;
		p.t_rad = 0.;
		p.t_absorb = 0.;
	}
	double time = 0.;
  double sum = 0.;
	for (int it=0; it<Nsteps; ++it){
		LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		for (auto & p : plist) sum += E0-p.p.t();
		dE.push_back(sum/Nparticles);

		time += dt;
		for (auto & p : plist){
      p.p = p0; // reset energy for probe test
			std::vector<fourvec> FS;
			int channel = update_particle_momentum(dt*fmc_to_GeV_m1, T,
				{0.0, 0.0, 0.0}, p.pid, (p.x.t()-p.t_rad)*fmc_to_GeV_m1, (p.x.t()-p.t_absorb)*fmc_to_GeV_m1, p.p, FS);

			p.freestream(dt);
			if (channel>=0) {
				p.p = FS[0];
				if (channel == 2 || channel ==3) p.t_rad = time;
				if (channel == 4 || channel ==5) p.t_absorb = time;
			}
		}
	}
	return dE;
}

std::vector<std::vector<double>> rate_test(double E0, double T, double dt=0.05, int Nsteps=100, int Nparticles=10000, std::string mode="old", double rescale=1.0){
  double fmc_to_GeV_m1 = 5.026;
	initialize(mode, "./settings.xml", 1.0);
	double M = 1.3;
	std::vector<std::vector<double>> Rate;
  Rate.resize(Nsteps);
	std::vector<particle> plist(Nparticles);
  fourvec p0{E0, 0, 0, std::sqrt(E0*E0-M*M)};
	for (auto & p : plist) {
		p.pid = 4;
		p.x = fourvec{0,0,0,0};
		p.p = p0;
		p.t_rad = 0.;
		p.t_absorb = 0.;
	}
	double time = 0.;
  double sum = 0.;
	for (int it=0; it<Nsteps; ++it){
    Rate[it].resize(6);
    for (int c=0; c<6; c++){
        Rate[it][c] = 0.;
    }
		LOG_INFO << it << " steps, " << "time = " << time << " [fm/c]";
		time += dt;
		for (auto & p : plist){
      p.p = p0; // reset energy for probe test
			std::vector<fourvec> FS;
			int channel = update_particle_momentum(dt*fmc_to_GeV_m1, T,
				{0.0, 0.0, 0.0}, p.pid, (p.x.t()-p.t_rad)*fmc_to_GeV_m1*rescale, (p.x.t()-p.t_absorb)*fmc_to_GeV_m1*rescale, p.p, FS);

			p.freestream(dt);
			if (channel>=0) {
				p.p = FS[0];
				if (channel == 2 || channel ==3) p.t_rad = time;
				if (channel == 4 || channel ==5) p.t_absorb = time;
        Rate[it][channel] += 1.;
			}
		}
    for (int c=0; c<6; c++){
        Rate[it][c] /= Nparticles;
    }
	}
	return Rate;
}
