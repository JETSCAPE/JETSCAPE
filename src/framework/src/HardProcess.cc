// Framework test (dummy) FluidDynamics class implementation (to be changed with real implemenation)
#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;

//#define CYAN    "\033[36m" 
#define MAGENTA "\033[35m"

HardProcess::HardProcess()
{
  VERBOSE(8);
  SetId("HardProcess");
}

HardProcess::~HardProcess()
{
  VERBOSE(8);
  hp_list.clear();
  disconnect_all();
}

void HardProcess::Init()
{
  JetScapeModuleBase::Init();

  INFO<<"Intialize HardProcess : "<<GetId()<< " ...";
 
  fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" );

  if (!fd)
     {
         WARN << "Not a valid JetScape XML Hard section file or no XML file loaded!";
          exit(-1);
     }
  
  VERBOSE(8);
  
  InitTask();
  
  JetScapeTask::InitTasks();
}

void HardProcess::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

  //VERBOSE(8);
  // DEBUG: Test ...
  /*
  double pAssign[4], xLoc[4];
  for (int i=0;i<=3; i++) {
     xLoc[i] = 0.0;
   };
   
   pAssign[0] = 11.0; 
   pAssign[3] = 10.0;
   pAssign[1] = pAssign[2] = 1.0;

   for (int i=0;i<1;i++)
     {
       hp_list.push_back(Parton(1,21,0,0,pAssign,xLoc));
     }
   
   VERBOSE(8)<<hp_list.size();
   VERBOSE(8)<<GetPartonAt(0).pid();
  */
  
  JetScapeTask::ExecuteTasks();
}

void HardProcess::Clear()
{
  DEBUG<<"Clear Hard Process : "<<GetId()<< " ...";
  //cout<<hp_list.size()<<endl;
  hp_list.clear();
  VERBOSE(8)<<hp_list.size();
}

void HardProcess::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  VERBOSE(8)<<w.lock()->GetOutputFileName();
  w.lock()->Write("HardProcess Parton List: "+GetId());
  for (int i=0;i<hp_list.size();i++)
    w.lock()->Write(GetPartonAt(i));
}
