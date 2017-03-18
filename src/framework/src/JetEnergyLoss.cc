// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include<iostream>
#include "JetScapeWriterAscii.h"

#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */

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
  SentInPartonsConnected=false;
  GetOutPartonsConnected=false;
  
  deltaT=0;
  maxT=0;
  
  inP=nullptr; pShower=nullptr;
  
  VERBOSE(8);
}

JetEnergyLoss::~JetEnergyLoss()
{
  VERBOSE(8);
  
  //pShower->clear();
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
  SetSentInPartonsConnected(false);
  SetGetOutPartonsConnected(false);
  
  //inP=j.inP;
  //pShower=j.pShower; // to be checked ...

  deltaT=j.deltaT;
  maxT=j.maxT;
  
  inP=nullptr; pShower=nullptr;
  
  VERBOSE(8) << "To be copied : # Subtasks = "<<j.GetTaskList().size();
  for (auto it : j.GetTaskList())
    {
      // Working via CRTP JetEnergyLossModule Clone function !
      auto st=dynamic_pointer_cast<JetEnergyLoss>(it)->Clone(); //shared ptr with clone !!????
      Add(st);
    }
  
}

void JetEnergyLoss::Clear()
{
  VERBOSESHOWER(8);
  if (pShower)
    pShower->clear();
  
  //inP=nullptr;pShower=nullptr; // kind of defeating the porpose of shared pointers somehow ...
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
  else
    {
      eloss->FirstChildElement("deltaT")->QueryDoubleText(&deltaT);
      eloss->FirstChildElement("maxT")->QueryDoubleText(&maxT);
      INFO<<"Eloss shower with deltaT = "<<deltaT<<" and maxT = "<<maxT;
    }
  
  if (GetNumberOfTasks()<1)
    {
      WARN << " : No valid Energy Loss modules found ...";
      exit(-1);
    }

  inP=nullptr;pShower=nullptr;
  
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();
}


void::JetEnergyLoss::DoShower()
{
  double tStart=0;
  double currentTime=0;

  VERBOSESHOWER(8)<<"Hard Parton from Initial Hard Process ...";
  VERBOSEPARTON(6,*GetShowerInitiatingParton());

  // consider pointers for speed up ... (!?)
  vector<Parton> pIn; vector<Parton> pOut;
  vector<Parton> pInTemp; vector<Parton> pOutTemp;
  vector<node> vStartVec; vector<node> vStartVecOut;
  
  pIn.push_back(*GetShowerInitiatingParton());

  // Add here the Hard Shower emitting parton ...
  vStart=pShower->new_vertex(make_shared<VertexBase>());
  vEnd=pShower->new_vertex(make_shared<VertexBase>());
  pShower->new_parton(vStart,vEnd,make_shared<Parton>(*GetShowerInitiatingParton()));

  // start then the recursive shower ...
  vStartVec.push_back(vEnd);
  
  do
    {
      VERBOSESHOWER(8)<<"Current time = "<<currentTime<<" with #Input "<<pIn.size();     
      currentTime += deltaT;
      
      for (int i=0;i<pIn.size();i++)
	{
	  //cout<<currentTime<<" "<<i<<" "<<pIn[i].pt()<<endl;
	  pInTemp.push_back(pIn[i]);
	  
	  SentInPartons(currentTime,pIn[i].pt(),pInTemp,pOutTemp);
	  //vStart=vStartVec[i+pIn.size()-1]; // could avoid the temps ... check ..
	  // works is always splitt, not with random ... fix below with temp vector ...!!!!
	  vStart=vStartVec[i];
	  //cout<<vStart<<endl;
	    
	  for (int k=0;k<pOutTemp.size();k++)
	    {
	      //cout<<vStartVec.size()-pIn.size()+(i)-k<<endl;
	      
	      vEnd=pShower->new_vertex(make_shared<VertexBase>(0,0,0,currentTime));
	      //cout<<vStart<<"-->"<<vEnd<<endl;
	      
	      pShower->new_parton(vStart,vEnd,make_shared<Parton>(pOutTemp[k]));
	      //cout<<pOutTemp[k];
	      
	      vStartVecOut.push_back(vEnd);
	      //cout<<vStartVec.size()<<endl;
	      
	      pOut.push_back(pOutTemp[k]);
	    }

	  //not working !? (check assignemt operator !??)
	  //pOut.insert(pOut.end(),pOutTemp.begin(),pOutTemp.end());
	 
	  pOutTemp.clear();
	  pInTemp.clear();
	}
      
      if (pOut.size()>0)
	{
	  pIn.clear();
	  pIn=pOut;
	  pOut.clear();
	  vStartVec.clear();
	  vStartVec=vStartVecOut;
	  vStartVecOut.clear();
	}
    }
  while (currentTime<maxT);

  pShower->PrintNodes();
  pShower->PrintEdges();
  // real ceck using mom. consveration at vertex ...!!!!
  //pShower->save(&(cout<<BOLDCYAN));
  
  pIn.clear();pOut.clear();pInTemp.clear();pOutTemp.clear();
}

void JetEnergyLoss::Exec()
{
  INFO<<"Run JetEnergyLoss ...";
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";

  if (GetShowerInitiatingParton())
     {
       pShower=make_shared<PartonShower>();
       DoShower();
     }
  else
    {WARN<<"NO Initial Hard Parton for Parton shower received ...";}  
  
  //JetScapeTask::ExecuteTasks(); // prevent Further modules to be execute, everything done by JetEnergyLoss ... (also set the no active flag ...!?)
}

void JetEnergyLoss::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  w.lock()->Write("Energy loss Shower Initating Parton: "+GetId());
  w.lock()->Write(inP);

  // check with gzip version later ...
  // Also allow standard output/not using GTL graph structure ....
  if (dynamic_pointer_cast<JetScapeWriterAscii> (w.lock()))
    {
      w.lock()->Write("Parton Shower in gml format from GTL:");
      pShower->save(dynamic_pointer_cast<JetScapeWriterAscii> (w.lock())->GetFileStream());
      w.lock()->Write("");
    }
  
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLoss::PrintShowerInitiatingParton()
{
  DEBUG<<inP->pid();
}
