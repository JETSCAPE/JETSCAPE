// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "HardProcess.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include <string>

#include<iostream>

using namespace std;

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
  
  JetScapeTask::ExecuteTasks();
}

void HardProcess::Clear()
{
  DEBUG<<"Clear Hard Process : "<<GetId()<< " ...";

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
