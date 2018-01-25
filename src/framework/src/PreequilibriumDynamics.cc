// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------

#include "PreequilibriumDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;

#define MAGENTA "\033[35m"

namespace Jetscape {

PreequilibriumDynamics::PreequilibriumDynamics()
{
  VERBOSE(8);
  SetId("PreequilibriumDynamics");
}

PreequilibriumDynamics::~PreequilibriumDynamics()
{
  VERBOSE(8);
  disconnect_all();
}

void PreequilibriumDynamics::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize PreequilibriumDynamics : "<<GetId()<< " ...";

  fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Preequilibrium" );

  if (!fd) {
     WARN << "Not a valid JetScape XML Preequilibrium Dynamics section file or no XML file loaded!";
	 exit(-1);
  }

  VERBOSE(8);

  //this is grabbing the initial entropy density ?
  ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
  if (!ini) {
      WARN << "No initialization module, try: auto trento = make_shared<TrentoInitial>(); jetscape->Add(trento);";
  }

  initialize_preequilibrium(parameter_list);

  InitTask();

  JetScapeTask::InitTasks();
}

void PreequilibriumDynamics::Exec()
{
  INFO <<"Run Preequilibrium : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

  if (ini) {
      INFO << "length of energy density vector=" << ini->energy_density_distribution_.size();
  }

  evolve_preequilibrium();

  JetScapeTask::ExecuteTasks();
}

//Freestreaming = no interactions, so the jets should not be modified by the medium in this stage
//if we exchange freestreaming for EKT or something including interactions, this needs to be revisited
/*
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
*/

} // end namespace Jetscape
