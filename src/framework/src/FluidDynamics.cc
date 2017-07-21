// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "FluidDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;
 
#define MAGENTA "\033[35m"

namespace Jetscape {

FluidDynamics::FluidDynamics()
{
  VERBOSE(8);
  eta=-99.99;
  SetId("FluidDynamics");
}

FluidDynamics::~FluidDynamics()
{
  VERBOSE(8);
  disconnect_all();
}

void FluidDynamics::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize FluidDynamics : "<<GetId()<< " ...";
 
  fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" );

  if (!fd) {
     WARN << "Not a valid JetScape XML Hydro section file or no XML file loaded!";
	 exit(-1);
  }
  
  VERBOSE(8);
  
  InitTask();

  initialize_hydro(parameter_list);

  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
      WARN << "No initialization module!";
  } else {
      INFO << "length of entropy density vector=" << ini->entropy_density_distribution_.size();
  }
  
  JetScapeTask::InitTasks();
}

void FluidDynamics::Exec()
{
  INFO<<"Run Hydro : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

  evolve_hydro();
  
  JetScapeTask::ExecuteTasks();
}

void FluidDynamics::UpdateEnergyDeposit(int t, double edop)
{
  //sigslot::lock_block<multi_threaded_local> lock(this);
  DEBUG<<MAGENTA<<"Jet Signal received : "<<t<<" "<<edop;
}

void FluidDynamics::GetEnergyDensity(int t,double &edensity)
{
  //sigslot::lock_block<multi_threaded_local> lock(this);
  edensity=0.5;
  DEBUG<<"Edensity to Jet = "<<edensity<<" at t="<<t;
}

} // end namespace Jetscape
