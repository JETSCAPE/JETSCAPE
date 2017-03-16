// Framework test (dummy) JetEnergyLoss class implementation (to be changed with real implemenation)
#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include<iostream>

using namespace std;

JetEnergyLoss::JetEnergyLoss()
{
  qhat=-99.99;
  SetId("JetEnergyLoss");
  jetSignalConnected=false;
  edensitySignalConnected=false;
  AddJetSourceSignalConnected=false;
  GetTemperatureSignalConnected=false;
  GetHydroCellSignalConnected=false;
  
  VERBOSE(8);
}

JetEnergyLoss::~JetEnergyLoss()
{
  VERBOSE(8);
  disconnect_all();
}

JetEnergyLoss::JetEnergyLoss(const JetEnergyLoss &j)
{
  qhat=j.GetQhat();
  SetActive(j.GetActive());
  SetId(j.GetId());
  SetJetSignalConnected(false);
  SetEdensitySignalConnected(false);
  SetAddJetSourceSignalConnected(false);
  SetGetTemperatureSignalConnected(false);
  SetGetHydroCellSignalConnected(false);
  inP=j.inP;
  
  VERBOSE(8) << "To be copied : # Subtasks = "<<j.GetTaskList().size();
  for (auto it : j.GetTaskList())
    {
      // Working via CRTP JetEnergyLossModule Clone function !
      auto st=dynamic_pointer_cast<JetEnergyLoss>(it)->Clone();
      Add(st);
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

  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();
}

void JetEnergyLoss::Exec()
{
  INFO<<"Run JetEnergyLoss ...";
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";

  // hmm, this is probably an other exampe for signals rather adding recursively ...
  for (auto it : GetTaskList())
    { 
      dynamic_pointer_cast<JetEnergyLoss>(it)->AddShowerInitiatingParton(inP);
    }
  
  JetScapeTask::ExecuteTasks();
}

void JetEnergyLoss::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  w.lock()->Write("Energy loss Shower Initating Parton: "+GetId());
  w.lock()->Write(inP);
  
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLoss::PrintShowerInitiatingParton()
{
  DEBUG<<inP->pid();
}
