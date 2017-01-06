// Framework test (dummy) FluidDynamics class implementation (to be changed with real implemenation)
#include "FluidDynamics.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

using namespace std;

//#define CYAN    "\033[36m" 
#define MAGENTA "\033[35m"

FluidDynamics::FluidDynamics()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Constructor called."<<endl;
  VERBOSE(8);
  eta=-99.99;
  SetId("FluidDynamics");
}

FluidDynamics::~FluidDynamics()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Destructor called."<<endl;
  VERBOSE(8);
  disconnect_all();
}

void FluidDynamics::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize FluidDynamics ...";
  
  tinyxml2::XMLElement *fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" );
  tinyxml2::XMLElement *hyd1=fd->FirstChildElement("Hydro1");

  if (hyd1)
    {
      string s = hyd1->FirstChildElement( "name" )->GetText();

      DEBUG << s << " to be initilizied ...";
      
      hyd1->FirstChildElement("eta")->QueryDoubleText(&eta);

      DEBUG << s << " with eta = "<<eta;      
    }

  //JetScapeSignalManager::Instance()->SetHydroPointer(shared_from_this());
  
  JetScapeTask::InitTasks();
}

void FluidDynamics::Exec()
{
  INFO<<"Run Hydro ...";
  //INFO<<"Found "<<GetNumberOfTasks()<<" Hydro Tasks/Modules Execute them ... ";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();
  // With current event number as static it is easy to define now for hydro event reading from file how often to be reused
  // can add aslo an other (static) variable for that in Fluiddynamics ...
  JetScapeTask::ExecuteTasks();
}

void FluidDynamics::UpdateEnergyDeposit(int t, double edop)
{
  sigslot::lock_block<multi_threaded_local> lock(this);
  DEBUG<<MAGENTA<<"Jet Signal received : "<<t<<" "<<edop;
}

void FluidDynamics::GetEnergyDensity(int t,double &edensity)
{
  sigslot::lock_block<multi_threaded_local> lock(this);
  edensity=0.5;
  DEBUG<<"Edensity to Jet = "<<edensity<<" at t="<<t;
}
