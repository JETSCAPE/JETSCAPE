// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "Lido.h"
#include "JetScapeLogger.h"
//#include "JetScapeXML.h"
#include <string>
#include <random>
#include <chrono>

#define MAGENTA "\033[35m"

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> ZeroOneDistribution(0.0, 1.0);

//double fmc_to_GeV_m1 = 5.026;

//#include "tinyxml2.h"
#include <iostream>

#include "FluidDynamics.h"

using namespace Jetscape;

using std::ifstream;
using std::ios;
using std::ofstream;
using std::ostream;

Lido::Lido()
{
  SetId("Lido");
  //VERBOSE(8);
}

Lido::~Lido()
{
  VERBOSE(8);
}

void Lido::Init()
{
  JSINFO << "Initialize Lido ...";
  std::string mode = "old";
  std::string path = "settings.xml";
  std::string table_path = "table.h5";
  boost::property_tree::ptree config;
  double mu = 2.0;
  double const_alphas = 0.3;

  initialize(mode, path, table_path, mu, -const_alphas, 0., 0., 0., 0., 0., 0.,
             1., 0.75);

  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss )     throw std::runtime_error("Eloss not properly initialized in XML file ...");

  Q0=1.0;
  hydro_Tc = 0.16;
  //martini->FirstChildElement("hydro_Tc")->QueryDoubleText(&hydro_Tc);
  hydro_tStart = 0.6;

}

void Lido::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton> &pIn, vector<Parton> &pOut)
{

  //VERBOSESHOWER(0)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;

  // particle info
  int Id, newId;
  FourVector pin;     // 4 vector of incoming parton
  FourVector pinrest; // 4 vector of incoming parton in rest frame of fluid cell
  FourVector pout;
  FourVector poutrest;

  FourVector xin;  // 4 vector for incoming position
  FourVector xout; // 4 vector for outgoing position (for next time step!)
  double eta;      // pseudo-rapidity

  // flow info
  double vx, vy, vz; // 3 components of flow velocity
  double T;          // Temperature of fluid cell

  //INFO<<"---------------------pin size: "<< pIn.size();

  for (int i = 0; i < pIn.size(); i++)
  {
    
    Id = pIn[i].pid();

    //VERBOSE(8) << "--------------------------------particle id: " << Id << " channel " << pIn[i].user_info<HQInfoBase>().hq_channel() << " mother id: " << pIn[i].user_info<HQInfoBase>().hq_mother_id();

    VERBOSE(0)<<"Lido: mass: "<<pIn[i].restmass()<<" time: "<<Time;    
    pin = FourVector(pIn[i].px(), pIn[i].py(), pIn[i].pz(), pIn[i].e());
    xin = FourVector(pIn[i].x_in().x(), pIn[i].x_in().y(), pIn[i].x_in().z(), Time);

    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

    double tt = xin.t();
    double xx = xin.x() + (Time-tt)*pin.x()/pin.t();
    double yy = xin.y() + (Time-tt)*pin.y()/pin.t();
    double zz = xin.z() + (Time-tt)*pin.z()/pin.t();

    GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    VERBOSE(0)<<"Inputs: "<<Time<<" "<<xx<<" "<< yy<<" "<< zz;
    VERBOSE(0) << "Temperature (Signal) = "
               << check_fluid_info_ptr->temperature << " " << check_fluid_info_ptr->vx << " " << check_fluid_info_ptr->vy << " " << check_fluid_info_ptr->vz;

    vx = check_fluid_info_ptr->vx;
    vy = check_fluid_info_ptr->vy;
    vz = check_fluid_info_ptr->vz;
    T = check_fluid_info_ptr->temperature;

    if (pIn[i].t() > Q0*Q0 + 0.00001 || Time <=hydro_tStart || T < hydro_Tc) continue;
    TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

    std::vector<fourvec> FS;
    std::vector<fourvec> FSG;
    particle Lido_particle;
    Lido_particle.mass = pIn[i].restmass();                                                                                        // mass
    Lido_particle.pid = Id;                                                                                                        // charm quark
    Lido_particle.x = fourvec{xin.t() * fmc_to_GeV_m1, xin.x() * fmc_to_GeV_m1, xin.y() * fmc_to_GeV_m1, xin.z() * fmc_to_GeV_m1}; //position
    Lido_particle.x0 = Lido_particle.x;
    Lido_particle.p = fourvec{pin.t(), pin.x(), pin.y(), pin.z()};
    Lido_particle.p0 = Lido_particle.p;
    Lido_particle.weight = 1.;
    if (pIn[i].has_user_info<LidoVirtualList>())
    {
      Lido_particle.radlist = pIn[i].user_info<LidoVirtualList>().radlist();
    }
    else
    {
      vector<particle> rad;
      Lido_particle.radlist = rad;
    }

    std::vector<particle> Lido_pOut;
    update_particle_momentum_Lido(deltaT * fmc_to_GeV_m1, T, {vx, vy, vz}, Lido_particle, Lido_pOut);
    //INFO <<"Time: "<<Time << " Temp: "<<T<<" vx: "<<vx<<" vy: "<<vy<<" vz: "<<vz<<" pin: "<< pin.t() << " " << pin.x() << " " << pin.y() << " " << pin.z()<<" mass: "<< sqrt(pin.t()*pin.t()-pin.x()*pin.x()-pin.y()*pin.y()-pin.z()*pin.z());
    //int channel = update_particle_momentum(deltaT*fmc_to_GeV_m1, T, {vx, vy, vz}, Id, (Time-pIn[i].form_time())*fmc_to_GeV_m1, (Time-pIn[i].absorp_time())*fmc_to_GeV_m1, fourvec{pin.t(),pin.x(),pin.y(),pin.z()}, FS);

    for (int i = 0; i < Lido_pOut.size(); i++)
    {
      int pid = Lido_pOut[i].pid;
      FourVector x = FourVector(Lido_pOut[i].x.x(), Lido_pOut[i].x.y(), Lido_pOut[i].x.z(), Lido_pOut[i].x.t());
      FourVector p = FourVector(Lido_pOut[i].p.x(), Lido_pOut[i].p.y(), Lido_pOut[i].p.z(), Lido_pOut[i].p.t());
      pOut.push_back(Parton(0, pid, 0, p, x));
      pOut.back().set_user_info(new LidoVirtualList(Lido_pOut[i].radlist));
    }
  }
}
