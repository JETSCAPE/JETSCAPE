// Framework test (dummy) JetEnergyLoss class implementation (to be changed with real implemenation)
#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "ELossModulesTest.h"
#include "tinyxml2.h"
#include<iostream>

using namespace std;

JetEnergyLoss::JetEnergyLoss()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Constructor called."<<endl;
  qhat=-99.99;
  SetId("JetEnergyLoss");
  jetSignalConnected=false;
  edensitySignalConnected=false;
  VERBOSE(8);
}

JetEnergyLoss::~JetEnergyLoss()
{
  //Simple Debug replace --> logger
  //cout<<"JetScape : Default Destructor called."<<endl;  
  VERBOSE(8);
  disconnect_all();
}

//JetEnergyLoss::JetEnergyLoss(const JetEnergyLoss &j) //check !!!!
JetEnergyLoss::JetEnergyLoss(const JetEnergyLoss &j)
{
  //VERBOSE(9);
  qhat=j.GetQhat();
  SetActive(j.GetActive());
  SetId(j.GetId());
  SetJetSignalConnected(false);//j.GetJetSignalConnected());
  SetEdensitySignalConnected(false);
  //cout<<qhat<<" "<<j.GetQhat()<<endl;
  //VERBOSE(9) << GetTaskList().size();
  //cout<<typeid(j).name()<<endl;
  
  VERBOSE(9) << "To be copied : # Subtasks = "<<j.GetTaskList().size();
  for (auto it : j.GetTaskList())
    {
      //cout<<typeid(it).name()<<endl;
      // std::dynamic_pointer_cast<A>(bar)
      //auto st=make_shared<dynamic_pointer_cast<it>> (*it);
      //auto st=make_shared<JetScapeTask> (*it);
      //auto st=make_shared<JetScapeTask> (*it);
      //cout<<dynamic_pointer_cast<JetScapeTask>(it)->GetQhat()<<endl;

      // Not happy !!!!!!
      // There must be a smarter way to get this from Base (see idea of Clone ...)!?
      // or use: string id from JetScapeTask ... (still a bit cluncky ...)
      if (dynamic_pointer_cast<Matter>(it))
	{
	  auto st=make_shared<Matter> (*dynamic_pointer_cast<Matter>(it));	
	  //cout<<st.get()<<" "<<it.get()<<endl;
	  //st = dynamic_pointer_cast<JetScapeTask>(it);
	  Add(st);
	}
      else if (	dynamic_pointer_cast<Martini>(it))
	{
	  auto st=make_shared<Martini> (*dynamic_pointer_cast<Martini>(it));       
	  //out<<st.get()<<" "<<it.get()<<endl;
	  //st = dynamic_pointer_cast<JetScapeTask>(it);
	  Add(st);
	}
    }
  
}

void JetEnergyLoss::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize JetEnergyLoss ..."; 
  
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );  

  if (!eloss)
    {
      WARN << " : Not a valid JetScape Energy Loss XML section in file!";
      exit(-1);
    }
  
  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss modules found ...";
      exit(-1);
    }

  // Has to be called explicitly since not really fully recursively yet (if ever needed)
  // So JetEnergyLoss is "Task Manager" of all energy loss modules
  // JetScapeTask::Init();
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

  //SetActive(false); // needed to access tasks not recursively by default but individually ...
                    // also usefull to prevent hydro if multiple read ins of the same event ...
}

void JetEnergyLoss::Exec()
{
  INFO<<"Run JetEnergyLoss ...";
  
  //test signal
  //DEBUG<<"Emit Signal: jetSignal(2,2.3)";
  //jetSignal(2,2.3);

  // Has to be called explicitly since not really fully recursively yet (if ever needed)
  // So JetEnergyLoss is "Task Manager" of all energy loss modules
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";
  JetScapeTask::ExecuteTasks();
  //GetTaskAt(0)->Exec(); // better with vector
  //GetTaskAt(1)->Exec();
  //Allows more felixbility in when used ...
}
