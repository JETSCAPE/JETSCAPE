// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include<iostream>

#include <thread>        
//#include <mutex>         
//#include <condition_variable>
//#include <future>

#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "JetScapeWriterAscii.h"
#include "JetScapeWriterHepMC.h"

//#include "PartonShowerGenerator.h"

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

  // consider pointers for speed up ... 
  vector<Parton> pIn; vector<Parton> pOut;
  vector<Parton> pInTemp; vector<Parton> pOutTemp;
  vector<Parton> pInTempModule;
  
  vector<node> vStartVec; vector<node> vStartVecOut;
  vector<node> vStartVecTemp;
  
  //DEBUG this guy isn't linked to anything - put in test particle for now
  //pIn.push_back(*GetShowerInitiatingParton());

  double pAssign[4] = {10,14,2,20};
  double xLoc[4] = {2,3,4,5};
  Parton pTemp(1,21,0,pAssign,xLoc);
  INFO<<"id: "<<pTemp.pid()<<" x: "<< pTemp.x_in().x() << " " << pTemp.x_in().y() << " " << pTemp.x_in().z() << " " << pTemp.x_in().t();
  INFO<< " p: "<< pTemp.p_in().x() << " " << pTemp.p_in().y() << " " << pTemp.p_in().z() << " " << pTemp.p_in().t();
  pIn.push_back(pTemp);

  // Add here the Hard Shower emitting parton ...
  vStart=pShower->new_vertex(make_shared<VertexBase>());
  vEnd=pShower->new_vertex(make_shared<VertexBase>());
  pShower->new_parton(vStart,vEnd,make_shared<Parton>(*GetShowerInitiatingParton()));

  // start then the recursive shower ...
  vStartVec.push_back(vEnd);
  //vStartVecTemp.push_back(vEnd);
  
  // ISSUE: Probably not yet 100% wrt to time step evolution ...
  // Logic mistake to remove the original ones when no split occured !!??? Follow up!!!!
  REMARK<<"DoShower() Splitting including time evolution (allowing non-splits at later times) implemeted correctly (should be made nicer/pointers). To be checked!!!";

  // --------------------------------------------
  
  do
    {
      VERBOSESHOWER(7)<<"Current time = "<<currentTime<<" with #Input "<<pIn.size();     
      currentTime += deltaT;

      // --------------------------------------------
      
      for (int i=0;i<pIn.size();i++)
	{
	  //DEBUG:
	  //cout<<currentTime<<" pIn size = "<<pIn.size()<<" "<<i<<" "<<pIn[i].pt()<<endl;

	  pInTemp.push_back(pIn[i]);
	  pInTempModule.push_back(pIn[i]);
	  
	  SentInPartons(currentTime,pIn[i].pt(),pInTempModule,pOutTemp);

	  vStart=vStartVec[i];
	  vStartVecTemp.push_back(vStart);

	  //DEBUG:
	  //cout<<vStart<<endl;
	  // --------------------------------------------
	  for (int k=0;k<pOutTemp.size();k++)
	    {
	      vEnd=pShower->new_vertex(make_shared<VertexBase>(0,0,0,currentTime));	    	      
	      pShower->new_parton(vStart,vEnd,make_shared<Parton>(pOutTemp[k]));	     
	      	     
	      //DEBUG:
	      //cout<<vStart<<"-->"<<vEnd<<endl;
	      //cout<<pOutTemp[k];
	      //cout<<vStartVec.size()<<endl;
	      //cout<<pInTempModule.size()<<endl;
	      
	      vStartVecOut.push_back(vEnd);
	      pOut.push_back(pOutTemp[k]);

	      // --------------------------------------------
	      // Add new roots from ElossModules ...
	      // (maybe add for clarity a new vector in the signal!???)
	      // Otherwise keep track of input size (so far always 1
	      // and check if size > 1 and create additional root nodes to that vertex ...
	      // Simple Test here below:
	      // DEBUG:
	      //cout<<"In JetEnergyloss : "<<pInTempModule.size()<<endl;
	      
	      if (pInTempModule.size()>1)
		{
		  VERBOSESHOWER(7)<<pInTempModule.size()-1<<" new root node(s) to be added ...";
		  //cout<<pInTempModule.size()-1<<" new root node(s) to be added ..."<<endl;
		  
		  for (int l=1;l<pInTempModule.size();l++)
		    {
		      node vNewRootNode=pShower->new_vertex(make_shared<VertexBase>(0,0,0,currentTime-deltaT));
		      pShower->new_parton(vNewRootNode,vEnd,make_shared<Parton>(pInTempModule[l]));
		    }
		}
	      // --------------------------------------------
	      
	      if (k==0)
		{
		  pInTemp.pop_back();       		  
		  vStartVecTemp.pop_back();		 
		}	  
	    }
	  // --------------------------------------------
	 
	  pOutTemp.clear();
	  pInTempModule.clear();
	}

      // --------------------------------------------
      
      pIn.clear();
      
      pIn.insert(pIn.end(),pInTemp.begin(),pInTemp.end());
      pIn.insert(pIn.end(),pOut.begin(),pOut.end());
      
      pOut.clear();
      pInTemp.clear();
          
      vStartVec.clear();
      
      vStartVec.insert(vStartVec.end(),vStartVecTemp.begin(),vStartVecTemp.end());
      vStartVec.insert(vStartVec.end(),vStartVecOut.begin(),vStartVecOut.end());
           
      vStartVecOut.clear();
      vStartVecTemp.clear();
    }
  while (currentTime<maxT); //other criteria (how to include; TBD)

  // --------------------------------------------
  
  pIn.clear();pOut.clear();pInTemp.clear();pOutTemp.clear();
  vStartVec.clear();
}

void JetEnergyLoss::Exec()
{
  INFO<<"Run JetEnergyLoss ...";
  INFO<<"Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";
  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" | Run JetEnergyLoss ...";
  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" | Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";
    
  if (GetShowerInitiatingParton())
     {
       pShower=make_shared<PartonShower>();

       /*
       //Check Memory ...
       VERBOSE(8)<<"Use PartonShowerGenerator to do Parton shower stored in PartonShower Graph class";
       DEBUG<<"Use PartonShowerGenerator to do Parton shower stored in PartonShower Graph class";
       
       PartonShowerGenerator PSG;
       PSG.DoShower(*shared_from_this()); //needed otherwise all signal slots have to be recreated for shower module ....
       // Overall not the nicest logic though .... Just to make changing and expanding the shower code in the future ...
       // (basically, just now to remove the code out of the jet energy loss class ...) TBD
       // also not really nice, since now the energy loss part in the parton shower and not really visible in this class ...
       // Keep both codes so far ...
       */
       
       // Shower handled in this class ...
       DoShower();
        
       pShower->PrintNodes();
       pShower->PrintEdges();
     }
  else
    {WARN<<"NO Initial Hard Parton for Parton shower received ...";}  

  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" Finished!";
  //JetScapeTask::ExecuteTasks(); // prevent Further modules to be execute, everything done by JetEnergyLoss ... (also set the no active flag ...!?)
}

void JetEnergyLoss::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  INFO<<"In JetEnergyLoss::WriteTask";
  w.lock()->WriteComment("Energy loss Shower Initating Parton: "+GetId());
  w.lock()->Write(inP);

  // check with gzip version later ...
  // Also allow standard output/not using GTL graph structure ....

  //If you want HepMC output, pass the whole shower along...
  //if (dynamic_pointer_cast<JetScapeWriterHepMC> (w.lock())){
     dynamic_pointer_cast<JetScapeWriterHepMC>(w.lock())->JetScapeWriterHepMC::WriteEvent(pShower);
  //}

  if (dynamic_pointer_cast<JetScapeWriterAscii> (w.lock()))
    {
      /*
      w.lock()->Write("Parton Shower in gml format from GTL:");
      pShower->save(dynamic_pointer_cast<JetScapeWriterAscii> (w.lock())->GetFileStream());
      w.lock()->Write("");
      */
      /*
      // Test: Each shower in seperate gml file for testing of JetScape format ...
      string fGMLname="shower_"; fGMLname += to_string(GetCurrentEvent()); fGMLname+=".gml";
      //cout<<fGMLname<<" "<<GetCurrentEvent()<<endl;
      pShower->save((char*) fGMLname.c_str());
      */
    }

  //Own storage of graph structure, needs separate PartonShower reader ...
  if (pShower)
    {
      w.lock()->WriteComment("Parton Shower in JetScape format to be used later by GTL graph:");
      
      // write vertices ...
      PartonShower::node_iterator nIt,nEnd;
      
      for (nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd; ++nIt)
	{
	  w.lock()->WriteWhiteSpace("["+to_string(nIt->id())+"]");
	  w.lock()->Write(pShower->GetVertex(*nIt));
	}
      
      PartonShower::edge_iterator eIt,eEnd;
      
      for (eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd; ++eIt)
	{
	  w.lock()->WriteWhiteSpace("["+to_string(eIt->source().id())+"]-->["+to_string(eIt->target().id())+"]");
	  w.lock()->Write(pShower->GetParton(*eIt));
	}
    }
  else
    {
      w.lock()->WriteComment("No EnergyLoss Modules were run - No Parton Shower information stored");
    }
  
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLoss::PrintShowerInitiatingParton()
{
  //DEBUG<<inP->pid();
}
