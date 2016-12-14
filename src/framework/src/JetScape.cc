// JetScape class implementation
#include "JetScape.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"

#include<iostream>

using namespace std;

JetScape::JetScape()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Constructor called."<<endl;
  n_events=1;
  VERBOSE(8);
}

JetScape::~JetScape()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Destructor called."<<endl;
  VERBOSE(8);
  //needed because of shared pointers and connection to instance
  // In general explore: unique_pointers rewrite!!!
  //JetScapeSignalManager::Instance()->Clear(); //not needed, use weak_ptr in JetScapeSignalManager class (=not owning)
}

void JetScape::Show()
{
  INFO_NICE;
  INFO_NICE<<" JetScape Event Generator";
  INFO_NICE<<" ------------------------";
  INFO_NICE;
}

void JetScape::Init()
{
  //Show();
  
  INFO<<"Intialize JetScape ...";
  
  JetScapeXML::Instance()->OpenXMLFile(GetXMLFileName());
  
  string log_debug = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "debug" )->GetText();
  string log_remark = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "remark" )->GetText();
  int m_vlevel=0;
  
  JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("vlevel")->QueryIntText(&m_vlevel);
  
  // Some general JetScape interface settings can also be read in from XML file ...
  // Order of task could also be verifed via XML or via sort of TaskList with order parameter
  if ((int) log_debug.find("off")>=0)
    JetScapeLogger::Instance()->SetDebug(false);
  
  if ((int) log_remark.find("on")>=0)
    JetScapeLogger::Instance()->SetRemark(true);

  DEBUG<<"JetScape Debug from XML = "<< log_debug;
  DEBUG<<"JetScape Remark from XML = "<< log_remark;

   if ((int) m_vlevel>0)
     {
       JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
       DEBUG<<"JetScape Verbose Level from XML = "<<m_vlevel;
     }

   SetPointers();
   
   // Has to be called explicitly since not really fully recursively (if ever needed)
   // So JetScape is "Task Manager" of all modules
   //JetScapeTask::Init();
   INFO<<"Found "<<GetNumberOfTasks()<<" Modules Initialize them ... ";
   JetScapeTask::InitTasks();
}

// kind of cluncky, maybe a better way ... !?
// Handle signal/slots in JetScape !?? (avoid passing pointers to sub tasks ...)
void JetScape::SetPointers()
{
  
   // to get hydro pointer for signals, use signal?
  INFO<<"Set Hydro Pointer for SignalManager to create Signal/Slots";
  
  //shared_ptr<FluidDynamics> m_hydro;
  //shared_ptr<JetEnergyLossManager> m_jloss_manager;
  
  for (auto it : GetTaskList())
    {
      if (dynamic_pointer_cast<FluidDynamics>(it))
	{
	  //m_hydro=dynamic_pointer_cast<FluidDynamics>(it);
	  JetScapeSignalManager::Instance()->SetHydroPointer(dynamic_pointer_cast<FluidDynamics>(it));
	}
    }

  for (auto it : GetTaskList())
	{
	  if (dynamic_pointer_cast<JetEnergyLossManager>(it))
	    {
	      //m_jloss_manager=dynamic_pointer_cast<JetEnergyLossManager>(it);	      
	      //dynamic_pointer_cast<JetEnergyLossManager>(it)->SetHydroPointer(m_hydro);
	      JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(dynamic_pointer_cast<JetEnergyLossManager>(it));
	    }
	}
  
  //JetScapeSignalManager::Instance()->SetHydroPointer(m_hydro);
  //JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(m_jloss_manager);
  // do it here and cluncky, but avoid issues with getting *this in shared pointers ...
  // add more, like IS, Hadronization ...
  
  }

void JetScape::Exec()
{
  INFO<<"Run JetScape ...";
  INFO<<BOLDBLACK<<"Number of Events = "<<GetNumberOfEvents(); 
  //SetCurrentEvent(0);
  
  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // So JetScape is "Task Manager" of all modules
  // Hmm, this change still not doing what I wanted (see old Init way ...)

  //SetPointers();
  
  for (int i=0;i<GetNumberOfEvents();i++)
    {
      INFO<<BOLDBLACK<<"Run Event # = "<<i;
      INFO<<"Found "<<GetNumberOfTasks()<<" Modules Execute them ... ";
      
      JetScapeTask::ExecuteTasks();
      IncrementCurrentEvent();
      
      //cout<<JetScapeSignalManager::Instance()->GetNumberOfJetSignals()<<endl;
      //JetScapeSignalManager::Instance()->PrintJetSignalMap();
      //JetScapeSignalManager::Instance()->PrintEdensitySignalMap();
    }
}
