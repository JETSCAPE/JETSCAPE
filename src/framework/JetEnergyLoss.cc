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
#include <iostream>
#include "tinyxml2.h"
#include "JetScapeSignalManager.h"
#include "JetScapeWriterAscii.h"
#include "HardProcess.h"

#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

//#include "PartonShowerGenerator.h"

#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */

using namespace std;

/**
   DELETE ME convenient output
*/

namespace Jetscape {
ostream &operator<<(ostream &ostr, const fjcore::PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
	 << " m = " << jet.m()
	 << " y = " << jet.rap()
	 << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }                                                      
  return ostr;
}


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

void JetEnergyLoss::DoShower()
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
  pIn.push_back(*GetShowerInitiatingParton());

  // Add here the Hard Shower emitting parton ...
  vStart=pShower->new_vertex(make_shared<Vertex>());
  vEnd=pShower->new_vertex(make_shared<Vertex>());
  // Add original parton later, after it had a chance to acquire virtuality
  // pShower->new_parton(vStart,vEnd,make_shared<Parton>(*GetShowerInitiatingParton()));

  // start then the recursive shower ...
  vStartVec.push_back(vEnd);
  //vStartVecTemp.push_back(vEnd);
  
  // --------------------------------------------

  // cerr << " ---------------------------------------------- " << endl;
  // cerr << "Start with " << *GetShowerInitiatingParton() << "  -> " << GetShowerInitiatingParton()->t() << endl;
  bool foundchangedorig=false;
  do
    {
      VERBOSESHOWER(7)<<"Current time = "<<currentTime<<" with #Input "<<pIn.size();     
      currentTime += deltaT;

      // --------------------------------------------
      
      for (int i=0;i<pIn.size();i++)
	{
	  // INFO << pIn.at(i).edgeid();
	  pInTempModule.push_back(pIn[i]);
	  
	  SentInPartons(deltaT,currentTime,pIn[i].pt(),pInTempModule,pOutTemp);
	  if ( !foundchangedorig ) {
	    // cerr  << " End with "<< pInTempModule.at(0) << "  -> " << pInTempModule.at(0).t() << endl;
	    // cerr << " ---------------------------------------------- " << endl;
	    pShower->new_parton(vStart,vEnd,make_shared<Parton>(pInTempModule.at(0)));
	    foundchangedorig=true;
	  }
	  
	  pInTemp.push_back(pInTempModule[0]);

	  vStart=vStartVec[i];
	  vStartVecTemp.push_back(vStart);

	  for (int k=0;k<pOutTemp.size();k++)
	    {
	      vEnd=pShower->new_vertex(make_shared<Vertex>(0,0,0,currentTime));	
	      int edgeid = pShower->new_parton(vStart,vEnd,make_shared<Parton>(pOutTemp[k]));
	      pOutTemp[k].set_shower( pShower );
	      pOutTemp[k].set_edgeid( edgeid );
		      
	      vStartVecOut.push_back(vEnd);
	      pOut.push_back(pOutTemp[k]);

	      Parton& particle = pOut.back();
	      // Parton& particle = pOut[iout];
	      // Parton& particle = pIn.at(i);
	      // if ( particle.parents().size() )
	      // 	cout << particle << "  " << particle.parents().at(0)<< endl;

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
		      node vNewRootNode=pShower->new_vertex(make_shared<Vertex>(0,0,0,currentTime-deltaT));
		      pShower->new_parton(vNewRootNode,vEnd,make_shared<Parton>(pInTempModule[l]));
		    }
		}
	      // --------------------------------------------
	      // 
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
  
  pIn.clear();pOut.clear();pInTemp.clear();pInTempModule.clear();pOutTemp.clear();
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
       JSDEBUG<<"Use PartonShowerGenerator to do Parton shower stored in PartonShower Graph class";
       
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

       weak_ptr<HardProcess> hproc = JetScapeSignalManager::Instance()->GetHardProcessPointer();

       for(unsigned int ipart=0; ipart<pShower->GetNumberOfPartons(); ipart++){
	 //   Uncomment to dump the whole parton shower into the parton container
	 //           hproc.lock()->AddParton(pShower->GetPartonAt(ipart));
       }
	
       shared_ptr<PartonPrinter> pPrinter = JetScapeSignalManager::Instance()->GetPartonPrinterPointer().lock();
       if ( pPrinter ){
	 pPrinter->GetFinalPartons2(pShower);
       }
    }
  else
    {WARN<<"NO Initial Hard Parton for Parton shower received ...";}  

  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" Finished!";
  //JetScapeTask::ExecuteTasks(); // prevent Further modules to be execute, everything done by JetEnergyLoss ... (also set the no active flag ...!?)
}

void JetEnergyLoss::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  VERBOSE(4)<<"In JetEnergyLoss::WriteTask";
  auto f = w.lock();
  if ( !f ) return;

  f->WriteComment("Energy loss Shower Initating Parton: "+GetId());
  f->Write(inP);

  // check with gzip version later ...
  // Also allow standard output/not using GTL graph structure ....

#ifdef USE_HEPMC
  //If you want HepMC output, pass the whole shower along...
  if (dynamic_pointer_cast<JetScapeWriterHepMC> (f)){
      VERBOSE(4) << " writing partons... found " << pShower->GetNumberOfPartons();
      (f)->Write(pShower);
  }
#endif

  //Own storage of graph structure, needs separate PartonShower reader ...
  if (pShower) {
    f->WriteComment("Parton Shower in JetScape format to be used later by GTL graph:");
    
    // write vertices
    PartonShower::node_iterator nIt,nEnd;
    
    for (nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd; ++nIt){ 
      f->WriteWhiteSpace("["+to_string(nIt->id())+"] V");
      f->Write(pShower->GetVertex(*nIt));
    }
    
    PartonShower::edge_iterator eIt,eEnd;      
    for (eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd; ++eIt) {
      f->WriteWhiteSpace("["+to_string(eIt->source().id())+"]=>["+to_string(eIt->target().id())+"] P");
      f->Write(pShower->GetParton(*eIt));
    }
  }
  else  {
    f->WriteComment("No EnergyLoss Modules were run - No Parton Shower information stored");
  }
  
  JetScapeTask::WriteTasks(w);
}

void JetEnergyLoss::PrintShowerInitiatingParton()
{
  //JSDEBUG<<inP->pid();
}



} // end namespace Jetscape
