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
  //set mode to old??
  std::string mode="old";
  boost::property_tree::ptree config;
  std::ifstream input("settings.xml");
  read_xml(input, config);
  double mu = config.get_child("Boltzmann.QCD").get<double>("mu");
  initialize_mD_and_scale(0, mu);

  MyProcesses.clear();
  MyProcesses.push_back( Rate22("Boltzmann/Qq2Qq", "settings.xml", dX_Qq2Qq_dt) );
  MyProcesses.push_back( Rate22("Boltzmann/Qg2Qg", "settings.xml", dX_Qg2Qg_dt) );
	
	MyProcesses.push_back( Rate23("Boltzmann/Qq2Qqg", "settings.xml", M2_Qq2Qqg) );
	MyProcesses.push_back( Rate23("Boltzmann/Qg2Qgg", "settings.xml", M2_Qg2Qgg) );

	BOOST_FOREACH(Process& r, MyProcesses){
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
				else continue;
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
				else continue;
				break;
			default:
				exit(-1);
				break;
		}
	} 
 
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
     INFO << "--------------------------------particle id: "<<Id;
     if (abs(Id) == 4)
       {
	  //INFO << "--------------------------------particle id: "<<Id;

	  pin = FourVector ( pIn[i].px(), pIn[i].py(), pIn[i].pz(), pIn[i].e());
	  xin = FourVector (pIn[i].x_in().x(),pIn[i].x_in().y(),pIn[i].x_in().z(), Time);

	  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
	  GetHydroCellSignal(Time, xin.x(), xin.y(), xin.z(), check_fluid_info_ptr);
	  VERBOSE(0)<<"Temperature from Brick (Signal) = "
		    <<check_fluid_info_ptr->temperature;

	  vx = check_fluid_info_ptr->vx;
	  vy = check_fluid_info_ptr->vy;
	  vz = check_fluid_info_ptr->vz;
           T = check_fluid_info_ptr->temperature;

          std::vector<fourvec> FS;
	  INFO << pin.t() << " " << pin.x() << " " << pin.y() << " " << pin.z();
	  int channel = update_particle_momentum(deltaT*fmc_to_GeV_m1, T, 
				{vx, vy, vz}, (Time-pIn[i].form_time())*fmc_to_GeV_m1, fourvec{pin.t(),pin.x(),pin.y(),pin.z()}, FS);     //needs to be modified?

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
          xout.Set(xin.t() + deltaT,xin.x() + pout.x()*a,xin.y() + pout.y()*a,xin.z() + pout.z()*a);

          //add outgoing heavy particle
          pOut.push_back(Parton(0, Id, 0, pout, xout));

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
            default: 
                      //INFO<<"No valid channel!";
                      break;
	            
	    }

       }
   }
}

int LBTD::update_particle_momentum(double dt, double temp, std::vector<double> v3cell, 
				   double D_formation_t, fourvec incoming_p, std::vector<fourvec> & FS)
{
	auto p_cell = incoming_p.boost_to(v3cell[0], v3cell[1], v3cell[2]);
	double D_formation_t_cell = D_formation_t / incoming_p.t() * p_cell.t();
	double dt_cell = dt / incoming_p.t() * p_cell.t();
	double E_cell = p_cell.t();
    std::vector<double> P_channels(MyProcesses.size());
	double P_total = 0.;
	int channel = 0;
	double dR;
	BOOST_FOREACH(Process& r, MyProcesses){
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
									{E_cell, temp, D_formation_t_cell}).s * dt_cell;
				else dR = 0.0;
				P_channels[channel] = P_total + dR;
				break;
			default:
				exit(-1);
				break;
		}
		P_total += dR;
		//LOG_INFO << dR;
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
	switch(MyProcesses[channel].which()){ 
		case 0:
			boost::get<Rate22>(MyProcesses[channel]).sample({E_cell, temp}, FS);
			break;
		case 1:
			boost::get<Rate23>(MyProcesses[channel]).sample(
											{E_cell, temp, D_formation_t_cell}, FS);
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
