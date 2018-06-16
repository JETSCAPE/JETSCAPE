// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "LBTD.h"
#include "JetScapeLogger.h"
//#include "JetScapeXML.h"
#include <string>
#include <random>
#include <chrono>
#include <math.h>
#include "simpleLogger.h"
#include <fstream>
#include "Srandom.h"
//#include "predefine.h"
#include "matrix_elements.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/variant/get.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>


unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> ZeroOneDistribution(0.0, 1.0);

//#include "tinyxml2.h"
#include <iostream>
#include <fstream>

#include "FluidDynamics.h"

const double fmc_to_GeV_m1 = 5.026;

using namespace Jetscape;

using std::ofstream;
using std::ifstream;
using std::ostream;
using std::ios;


LBTD::LBTD()
{
  SetId("LBTD");
  //VERBOSE(8);
}

LBTD::~LBTD()
{
  VERBOSE(8);
}

void LBTD::Init()
{
  INFO<<"Initialize LBTD ...";
  std::string mode="old";
  std::string path="./settings.xml";
  boost::property_tree::ptree config;
  double mu=1.0;
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

void LBTD::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  
  // particle info
  int Id, newId;
  FourVector pin;            // 4 vector of incoming parton
  FourVector pinrest;        // 4 vector of incoming parton in rest frame of fluid cell
  FourVector pout;           
  FourVector poutrest;       
      
  FourVector xin;            // 4 vector for incoming position
  FourVector xout;           // 4 vector for outgoing position (for next time step!)
  double eta;                // pseudo-rapidity

  // flow info
  double vx, vy, vz;         // 3 components of flow velocity
  double T;                  // Temperature of fluid cell

  //INFO<<"---------------------pin size: "<< pIn.size();

  for (int i=0;i<pIn.size();i++)
   {
     Id = pIn[i].pid();

     if (abs(Id) == 4||abs(Id) == 5)
       {
	  INFO << "--------------------------------particle id: " << Id << " channel " <<pIn[i].hq_channel() << " mother id: "<< pIn[i].hq_mother_id();

	  pin = FourVector ( pIn[i].px(), pIn[i].py(), pIn[i].pz(), pIn[i].e());
	  xin = FourVector (pIn[i].x_in().x(),pIn[i].x_in().y(),pIn[i].x_in().z(), Time);

	  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
	  GetHydroCellSignal(Time, xin.x(), xin.y(), xin.z(), check_fluid_info_ptr);
          //VERBOSE(7)<<"Inputs: "<<Time<<" "<<xin.x()<<" "<< xin.y()<<" "<< xin.z();
	  //VERBOSE(0)<<"Temperature (Signal) = "
	  //<<check_fluid_info_ptr->temperature<<" "<<check_fluid_info_ptr->vx<<" "<<check_fluid_info_ptr->vy<<" "<<check_fluid_info_ptr->vz;

	  vx = check_fluid_info_ptr->vx;
	  vy = check_fluid_info_ptr->vy;
	  vz = check_fluid_info_ptr->vz;
           T = check_fluid_info_ptr->temperature;

          std::vector<fourvec> FS;
	  //INFO <<"Input for update_particle_momentum: T: "<<T<<" vx: "<<vx<<" vy: "<<vy<<" vz: "<<vz<<" pin: "<< pin.t() << " " << pin.x() << " " << pin.y() << " " << pin.z()<<" time: "<<Time-pIn[i].form_time()<<" "<<Time-pIn[i].absorp_time();
	  int channel = update_particle_momentum(deltaT*fmc_to_GeV_m1, T, {vx, vy, vz}, Id, (Time-pIn[i].form_time())*fmc_to_GeV_m1, (Time-pIn[i].absorp_time())*fmc_to_GeV_m1, fourvec{pin.t(),pin.x(),pin.y(),pin.z()}, FS);

          FourVector p1out;
          FourVector p2out;
          int qId=random_choose_Id();
          if (channel>=0) 
          {
	    pout.Set(FS[0].x(),FS[0].y(),FS[0].z(), FS[0].t());  
	    p1out.Set(FS[1].x(),FS[1].y(),FS[1].z(), FS[1].t());  
	  }
	  else{
	    pout.Set(pin.x(), pin.y(), pin.z(), pin.t());
	  }

          //free streaming
          double a = deltaT/pout.t();
          //INFO<<"a: "<<a;
          xout.Set(xin.x() + pout.x()*a, xin.y() + pout.y()*a, xin.z() + pout.z()*a, xin.t() + deltaT);
          //INFO<<"x out: "<<xout.t()<<" "<<xout.x()<<" "<<xout.y()<<" "<<xout.z();

          //add outgoing heavy particle
          pOut.push_back(Parton(0, Id, 0, pout, xout));
          pOut[pOut.size()-1].set_form_time(pIn[i].form_time());
          pOut[pOut.size()-1].set_absorp_time(pIn[i].absorp_time());
          pOut[pOut.size()-1].set_hq_channel(pIn[i].hq_channel());
          pOut[pOut.size()-1].set_hq_mother_id(pIn[i].hq_mother_id());
          

          //add high energy light partons
          
	  switch(channel)
	  {
	    case 0: {
		      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,qId,0,xout,p1out));
                      }
                      break;
                    }
            case 1: {  
		      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,21,0,xout,p1out));
                      }
                      break;
                    }
            case 2: {
                      pOut[pOut.size()-1].set_form_time(Time);
	              p2out.Set(FS[2].x(),FS[2].y(),FS[2].z(), FS[2].t());
                      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,qId,0,xout,p1out));
                      }
                      if(p2out.t()>T)
                      {
                        pOut.push_back(Parton(0,21,0,xout,p2out));
                      }
                      break;
                    }
            case 3: {
                      pOut[pOut.size()-1].set_form_time(Time);
	              p2out.Set(FS[2].x(),FS[2].y(),FS[2].z(), FS[2].t());
                      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,21,0,xout,p1out));
                      }
                      if(p2out.t()>T)
                      {
                        pOut.push_back(Parton(0,21,0,xout,p2out));
                      }
                      break;
                    }
            case 4: {
                      pOut[pOut.size()-1].set_absorp_time(Time);
                      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,qId,0,xout,p1out));
                      }
                      break;
                    }
            case 5: {
                      pOut[pOut.size()-1].set_absorp_time(Time);
                      if(p1out.t()>T)
                      {
                        pOut.push_back(Parton(0,21,0,xout,p1out));
                      }
                      break;
                    }
            default: 
                      //INFO<<"No valid channel!";
                      break;
	            
	    }
            

       }
     //not a heavy quark
     else
       {
	 pOut.push_back(pIn[i]);
       }
   }
}

int LBTD::update_particle_momentum(double dt, double temp, std::vector<double> v3cell,
			int pid, double D_formation_t23, double D_formation_t32, fourvec incoming_p, std::vector<fourvec> & FS)
{
	int absid = std::abs(pid);
	auto p_cell = incoming_p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	//INFO<<"p_cell: "<< p_cell;
	//INFO<<"incoming_p: "<<incoming_p;
	double D_formation_t23_cell = D_formation_t23 / incoming_p.t() * p_cell.t();
	double D_formation_t32_cell = D_formation_t32 / incoming_p.t() * p_cell.t();
	double dt_cell = dt / incoming_p.t() * p_cell.t();
	double E_cell = p_cell.t();
	//INFO<<"E_cell: "<<E_cell<<"dt_cell: "<<dt_cell;
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
		//INFO<<"dR: "<<dR;
		P_total += dR;
		channel ++;
	}
	for(auto& item : P_channels) {item /= P_total;}
	if (P_total > 0.15) LOG_WARNING << "P_total = " << P_total << " may be too large";
	//LOG_WARNING << "P_total = " << P_total;
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

void LBTD::init_process(Process& r, std::string mode){
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

int LBTD::random_choose_Id()
{
  double x=ZeroOneDistribution(generator);
  if(x<1./3.) 
  {
    return 1;
  } 
  else if(x<2./3.) 
  {
    return 2;
  }
  else 
  {
    return 3;
  }
}
