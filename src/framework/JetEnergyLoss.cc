/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <iostream>
#include <thread>
//#include <mutex>
//#include <condition_variable>
//#include <future>

#include "JetEnergyLoss.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include "tinyxml2.h"
#include "JetScapeSignalManager.h"
#include "JetScapeWriterStream.h"
#include "HardProcess.h"
#include "JetScapeModuleMutex.h"
#include "LiquefierBase.h"
#include "MakeUniqueHelper.h"
#include "FluidDynamics.h"
#include <GTL/dfs.h>

#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif

#define BOLDCYAN "\033[1m\033[36m" /* Bold Cyan */

using namespace std;

namespace Jetscape {

JetEnergyLoss::JetEnergyLoss() {
  qhat = -99.99;
  SetId("JetEnergyLoss");
  jetSignalConnected = false;
  edensitySignalConnected = false;
  GetHydroCellSignalConnected = false;
  SentInPartonsConnected = false;

  deltaT = 0;
  maxT = 0;

  inP = nullptr;
  pShower = nullptr;

  VERBOSE(8);
}

JetEnergyLoss::~JetEnergyLoss() {
  VERBOSE(8);

  //pShower->clear();
  disconnect_all();
}

JetEnergyLoss::JetEnergyLoss(const JetEnergyLoss &j) {
  qhat = j.GetQhat();
  SetActive(j.GetActive());
  SetId(j.GetId());
  SetJetSignalConnected(false);
  SetEdensitySignalConnected(false);
  SetGetHydroCellSignalConnected(false);
  SetSentInPartonsConnected(false);
  AddModuleClock(j.GetModuleClock()); //JP: Check if memory leak ... should not due to shared_ptr usage, but confirm ...

  deltaT = j.deltaT;
  maxT = j.maxT;
  startT = j.startT;

  inP = nullptr;
  pShower = nullptr;

  SetUseIntialPartonShower(j.GetUseIntialPartonShower());
  AddPartonShowerGenerator(j.GetPartonShowerGenerator());

  VERBOSE(8) << "To be copied : # Subtasks = " << j.GetTaskList().size();
  for (auto it : j.GetTaskList()) {
    // Working via CRTP JetEnergyLossModule Clone function !
    auto st = dynamic_pointer_cast<JetEnergyLoss>(it)->Clone(); //shared ptr with clone !!????
    Add(st);
  }
}

void JetEnergyLoss::Clear() {
  VERBOSESHOWER(8);
  if (pShower)
    pShower->clear();

  this->final_Partons.clear();

  inP = nullptr;

  pIn.clear();
  vStartVec.clear();
}

void JetEnergyLoss::Init() {
  JetScapeModuleBase::Init();

  JSINFO << "Intialize JetEnergyLoss ...";

  deltaT = GetXMLElementDouble({"Eloss", "deltaT"});

  maxT = GetXMLElementDouble({"Eloss", "maxT"});

  if (GetActive())
    JSINFO << "Eloss shower with deltaT = " << deltaT << " and maxT = " << maxT;
  else
    JSINFO << "Eloss shower via Main Clock ...";

  std::string mutexOnString = GetXMLElementText({"Eloss", "mutex"}, false);
  if (!mutexOnString.compare("ON"))
    //Check mutual exclusion of Eloss Modules
  {
    if (GetNumberOfTasks() > 1) {
      for (auto elossModule : GetTaskList()) {
        shared_ptr<JetScapeModuleMutex> mutex_ptr = elossModule->GetMutex();
        if (mutex_ptr) {
          if (!(mutex_ptr->CheckMutex(GetTaskList()))) {
            JSWARN << "Mutual exclusive Energy-Loss modules attached together!";
            throw std::runtime_error("Fix it by attaching one of them.");
          }
        }
      }
    }
  }

  if (GetNumberOfTasks() < 1) {
    JSWARN << " : No valid Energy Loss modules found ...";
    exit(-1);
  }

  inP = nullptr;
  pShower = nullptr;

  JSINFO << "Found " << GetNumberOfTasks()
         << " Eloss Tasks/Modules Initialize them ... ";

  JetScapeTask::InitTasks();
}

// JP: with all changes check for memory leaks ...

void JetEnergyLoss::DoShower() {
  double tStart = startT;
  double currentTime = 0;

  VERBOSESHOWER(8) << "Hard Parton from Initial Hard Process ...";
  VERBOSEPARTON(6, *GetShowerInitiatingParton());

  // consider pointers for speed up ...
  //vector<Parton> pIn;
  // DEBUG this guy isn't linked to anything - put in test particle for now
  pIn.push_back(*GetShowerInitiatingParton());

  //vector<node> vStartVec;
  // Add here the Hard Shower emitting parton ...
  vStart = pShower->new_vertex(make_shared<Vertex>());
  vEnd = pShower->new_vertex(make_shared<Vertex>());
  // Add original parton later, after it had a chance to acquire virtuality
  // pShower->new_parton(vStart,vEnd,make_shared<Parton>(*GetShowerInitiatingParton()));

  // start then the recursive shower ...
  vStartVec.push_back(vEnd);

  // cerr << " ---------------------------------------------- " << endl;
  // cerr << "Start with " << *GetShowerInitiatingParton()
  //      << "  -> " << GetShowerInitiatingParton()->t() << endl;
  /*
  bool foundchangedorig = false;
  int droplet_stat = -11;
  int miss_stat = -13;
  int neg_stat = -17;
  */
  if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
    droplet_stat = liquefier_ptr.lock()->get_drop_stat();
    miss_stat = liquefier_ptr.lock()->get_miss_stat();
    neg_stat = liquefier_ptr.lock()->get_neg_stat();
  }
  do {
    vector<Parton> pOut;
    vector<Parton> pInTemp;

    vector<node> vStartVecOut;
    vector<node> vStartVecTemp;

    VERBOSESHOWER(7) << "Current time = " << currentTime << " with #Input "
                     << pIn.size();
    currentTime += deltaT;

    for (int i = 0; i < pIn.size(); i++) {
      vector<Parton> pInTempModule;
      vector<Parton> pOutTemp;
      // JSINFO << pIn.at(i).edgeid();
      pInTempModule.push_back(pIn[i]);
      SentInPartons(deltaT, currentTime, pIn[i].pt(), pInTempModule, pOutTemp);

      // apply liquefier
      if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
        liquefier_ptr.lock()->add_hydro_sources(pInTempModule, pOutTemp);
      }

      // stuffs related to vertex
      if (!foundchangedorig) {
        // cerr << " End with "<< pInTempModule.at(0) << "  -> "
        //      << pInTempModule.at(0).t() << endl;
        // cerr << " ---------------------------------------------- "
        //      << endl;
        pShower->new_parton(vStart, vEnd,
                            make_shared<Parton>(pInTempModule.at(0)));
        foundchangedorig = true;
      }

      vStart = vStartVec[i];
      if (pOutTemp.size() == 0) {
        // no need to generate a vStart for photons and liquefied
        // partons
        if (pInTempModule[0].pstat() != droplet_stat &&
            pInTempModule[0].pstat() != miss_stat &&
            pInTempModule[0].pstat() != neg_stat &&
            !pInTempModule[0].isPhoton(pInTempModule[0].pid())) {
          vStartVecTemp.push_back(vStart);
        }
      } else if (pOutTemp.size() == 1) {
        // no need to generate a vStart for photons and liquefied
        // partons
        if (pOutTemp[0].pstat() != droplet_stat &&
            pOutTemp[0].pstat() != miss_stat &&
            pOutTemp[0].pstat() != neg_stat &&
            !pOutTemp[0].isPhoton(pOutTemp[0].pid())) {
          vStartVecTemp.push_back(vStart);
        }
      } else {
        for (int k = 0; k < pOutTemp.size(); k++) {
          int edgeid = 0;
          if (pOutTemp[k].pstat() == neg_stat) {
            node vNewRootNode = pShower->new_vertex(
                                  make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
            edgeid = pShower->new_parton(vNewRootNode, vStart,
                                         make_shared<Parton>(pOutTemp[k]));
          } else {
            vEnd =
              pShower->new_vertex(make_shared<Vertex>(0, 0, 0, currentTime));
            edgeid = pShower->new_parton(vStart, vEnd,
                                         make_shared<Parton>(pOutTemp[k]));
          }
          pOutTemp[k].set_shower(pShower);
          pOutTemp[k].set_edgeid(edgeid);

          // no need to generate a vStart for photons and liquefied
          // partons
          if (pOutTemp[k].pstat() != droplet_stat &&
              pOutTemp[k].pstat() != miss_stat &&
              pOutTemp[k].pstat() != neg_stat &&
              !pOutTemp[k].isPhoton(pOutTemp[k].pid())) {
            vStartVecOut.push_back(vEnd);
          }

          // --------------------------------------------
          // Add new roots from ElossModules ...
          // (maybe add for clarity a new vector in the signal!???)
          // Otherwise keep track of input size (so far always 1
          // and check if size > 1 and create additional root nodes to that vertex ...
          // Simple Test here below:
          // DEBUG:
          //cout<<"In JetEnergyloss : "<<pInTempModule.size()<<end;
          if (pInTempModule.size() > 1) {
            VERBOSE(7) << pInTempModule.size() - 1
                       << " new root node(s) to be added ...";
            //cout << pInTempModule.size()-1
            //     << " new root node(s) to be added ..." << endl;

            for (int l = 1; l < pInTempModule.size(); l++) {
              node vNewRootNode = pShower->new_vertex(
                                    make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
              pShower->new_parton(vNewRootNode, vEnd,
                                  make_shared<Parton>(pInTempModule[l]));
            }
          }
        }
      }

      // update parton shower
      if (pOutTemp.size() == 0) {
        // this is the free-streaming case for MATTER
        // do not push back droplets
        if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
          if (pInTempModule[0].pstat() == droplet_stat)
            continue;
          if (pInTempModule[0].pstat() == miss_stat)
            continue;
          if (pInTempModule[0].pstat() == neg_stat)
            continue;
        }
        // do not push back photons
        if (pInTempModule[0].isPhoton(pInTempModule[0].pid()))
          continue;
        pInTemp.push_back(pInTempModule[0]);
      } else if (pOutTemp.size() == 1) {
        // this is the free-streaming case for MARTINI or LBT
        // do not push back droplets
        if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
          if (pOutTemp[0].pstat() == droplet_stat)
            continue;
          if (pOutTemp[0].pstat() == miss_stat)
            continue;
          if (pOutTemp[0].pstat() == neg_stat)
            continue;
        }
        // do not push back photons
        if (pOutTemp[0].isPhoton(pOutTemp[0].pid()))
          continue;
        pInTemp.push_back(pOutTemp[0]);
      } else {
        for (int k = 0; k < pOutTemp.size(); k++) {
          // do not push back droplets
          if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
            if (pOutTemp[k].pstat() == droplet_stat)
              continue;
            if (pOutTemp[k].pstat() == miss_stat)
              continue;
            if (pOutTemp[k].pstat() == neg_stat)
              continue;
          }
          // do not push back missing (from AdSCFT)
          if (pOutTemp[k].pstat() == miss_stat)
            continue;
          // do not push back photons
          if (pOutTemp[k].isPhoton(pOutTemp[k].pid()))
            continue;

          pOut.push_back(pOutTemp[k]);
        }
      }
    }

    // one time step is finished, now update parton shower to pIn
    pIn.clear();
    pIn.insert(pIn.end(), pInTemp.begin(), pInTemp.end());
    pIn.insert(pIn.end(), pOut.begin(), pOut.end());

    // update vertex vector
    vStartVec.clear();
    vStartVec.insert(vStartVec.end(), vStartVecTemp.begin(),
                     vStartVecTemp.end());
    vStartVec.insert(vStartVec.end(), vStartVecOut.begin(), vStartVecOut.end());
  } while (currentTime < maxT); // other criteria (how to include; TBD)
}

void JetEnergyLoss::Exec() {
  VERBOSE(1) << "Run JetEnergyLoss ...";
  VERBOSE(1) << "Found " << GetNumberOfTasks()
             << " Eloss Tasks/Modules Execute them ... ";
  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" | Run JetEnergyLoss ...";
  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" | Found "<<GetNumberOfTasks()<<" Eloss Tasks/Modules Execute them ... ";

  //if (GetShowerInitiatingParton()) {
  if (GetShowerInitiatingParton() && !useShower) {
    pShower = make_shared<PartonShower>();

    if (!GetPartonShowerGenerator()) {
      //VERBOSE(8)<<"Use Default DoShower() to do Parton shower stored in PartonShower Graph class";
      JSDEBUG<<"--> Use Default DoShower() to do Parton shower stored in PartonShower Graph class";

      // Shower handled in this class ...
      DoShower();
    }
    else {
      JSDEBUG<<"--> Use PartonShowerGenerator to do Parton shower stored in PartonShower Graph class";

      GetPartonShowerGenerator()->DoShower(*dynamic_pointer_cast<JetEnergyLoss>(shared_from_this()));
    }
  }
  else if (GetInitialPartonShower() && useShower)
    {
      if (GetPartonShowerGenerator())
        {
          pShower=make_shared<PartonShower>();

          JSDEBUG<<"--> Use DoShower() provided from PSG reading in full intial shower to do Parton shower stored in PartonShower Graph class";

	        GetPartonShowerGenerator()->DoShower(*dynamic_pointer_cast<JetEnergyLoss>(shared_from_this()));
        }
      else
        {JSWARN<<"No proper external Parton Shower Generator attached ..."; exit(-1);}
    }
  else {
    JSWARN << "NO Initial Hard Parton/or (ISR) shower for Parton shower received ...";exit(-1);
  }

  pShower->PrintNodes();
  pShower->PrintEdges();

  //REMARK JP: No idea what and why code below is needed !!!!
  //           Discuss and clean up in the future !!!!

  weak_ptr<HardProcess> hproc =
    JetScapeSignalManager::Instance()->GetHardProcessPointer();

  for (unsigned int ipart = 0; ipart < pShower->GetNumberOfPartons();
       ipart++) {
    //   Uncomment to dump the whole parton shower into the parton container
    // auto hp = hproc.lock();
    // if ( hp ) hp->AddParton(pShower->GetPartonAt(ipart));
  }

  shared_ptr<PartonPrinter> pPrinter =
    JetScapeSignalManager::Instance()->GetPartonPrinterPointer().lock();
  if (pPrinter) {
    pPrinter->GetFinalPartons(pShower);
  }

  shared_ptr<JetEnergyLoss> pEloss =
    JetScapeSignalManager::Instance()->GetEnergyLossPointer().lock();
  if (pEloss) {
    pEloss->GetFinalPartonsForEachShower(pShower);
  }

  //DEBUGTHREAD<<"Task Id = "<<this_thread::get_id()<<" Finished!";
  //JetScapeTask::ExecuteTasks(); // prevent Further modules to be execute, everything done by JetEnergyLoss ... (also set the no active flag ...!?)
}

void JetEnergyLoss::InitPerEvent()
{
  VERBOSE(3) << "InitPerEvent() for usage per time step ...";

  if (GetShowerInitiatingParton())
  {
    pShower = make_shared<PartonShower>();
    pIn.push_back(*GetShowerInitiatingParton());

    vStart = pShower->new_vertex(make_shared<Vertex>());
    vEnd = pShower->new_vertex(make_shared<Vertex>());

    // start then the recursive shower ...
    vStartVec.push_back(vEnd);

    if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
    droplet_stat = liquefier_ptr.lock()->get_drop_stat();
    miss_stat = liquefier_ptr.lock()->get_miss_stat();
    neg_stat = liquefier_ptr.lock()->get_neg_stat();
    }

  } else {
    JSWARN << "NO Initial Hard Parton for Parton shower received ...";exit(-1);
  }

}

void JetEnergyLoss::FinishPerEvent()
{
  VERBOSE(3) << "FinishPerEvent() for usage per time step ...";

  pShower->PrintNodes();
  pShower->PrintEdges();

  weak_ptr<HardProcess> hproc =
    JetScapeSignalManager::Instance()->GetHardProcessPointer();


  //REMARK JP: No idea what and why code below is needed !!!!
  //           Discuss and clean up in the future !!!!

  for (unsigned int ipart = 0; ipart < pShower->GetNumberOfPartons();
       ipart++) {
    //   Uncomment to dump the whole parton shower into the parton container
    // auto hp = hproc.lock();
    // if ( hp ) hp->AddParton(pShower->GetPartonAt(ipart));
  }

  shared_ptr<PartonPrinter> pPrinter =
    JetScapeSignalManager::Instance()->GetPartonPrinterPointer().lock();
  if (pPrinter) {
    pPrinter->GetFinalPartons(pShower);
  }

  shared_ptr<JetEnergyLoss> pEloss =
    JetScapeSignalManager::Instance()->GetEnergyLossPointer().lock();
  if (pEloss) {
    pEloss->GetFinalPartonsForEachShower(pShower);
  }

  //JP: Quick fix, to be discussed, similar to writer, clear is only called for active tasks, so call here directly ...

  Clear();
}

void JetEnergyLoss::CalculateTime()
{
  //VERBOSE(3) << "Calculate JLoss per time step ... Current (Module) Time = " << GetModuleCurrentTime() << " dT = " << GetModuleDeltaT();
  //ClockInfo();
}

//Put shower generation in exec for now ... to be discussed if correct and/or
//if in general Calc and Exec is needed. Adavantage, if multi-threaded and exec after cals, no issues with
//data not availabe at exec time ...

//JP: Exact copy of DoShower() for now, for test purpose only, should be
//seperated out in an extra function so that no duplciate code needed!!!

// Think about how to do this modular with PSG !????

void JetEnergyLoss::ExecTime()
{
  VERBOSE(3) << "Execute JLoss per time step ... Current (Module) Time = " << GetModuleCurrentTime() << " dT = " << GetModuleDeltaT();

  double currentTime = GetModuleCurrentTime();
  double moduleDeltaT = GetModuleDeltaT();

  if (GetShowerInitiatingParton())
  {
    vector<Parton> pOut;
    vector<Parton> pInTemp;

    vector<node> vStartVecOut;
    vector<node> vStartVecTemp;

    VERBOSESHOWER(7) << "Current time = " << currentTime << " with #Input "
    << pIn.size();

    //JP: Check if this usage is consistent with master clock etc start time ... !!!!

    currentTime += moduleDeltaT;

    for (int i = 0; i < pIn.size(); i++) {
      vector<Parton> pInTempModule;
      vector<Parton> pOutTemp;
      // JSINFO << pIn.at(i).edgeid();
      pInTempModule.push_back(pIn[i]);
      SentInPartons(moduleDeltaT, currentTime, pIn[i].pt(), pInTempModule, pOutTemp);

      // apply liquefier
      if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
        liquefier_ptr.lock()->add_hydro_sources(pInTempModule, pOutTemp);
      }

      // stuffs related to vertex
      if (!foundchangedorig) {
        // cerr << " End with "<< pInTempModule.at(0) << "  -> "
        //      << pInTempModule.at(0).t() << endl;
        // cerr << " ---------------------------------------------- "
        //      << endl;
        pShower->new_parton(vStart, vEnd,
                            make_shared<Parton>(pInTempModule.at(0)));
        foundchangedorig = true;
      }

      vStart = vStartVec[i];
      if (pOutTemp.size() == 0) {
        // no need to generate a vStart for photons and liquefied
        // partons
        if (pInTempModule[0].pstat() != droplet_stat &&
            pInTempModule[0].pstat() != miss_stat &&
            pInTempModule[0].pstat() != neg_stat &&
            !pInTempModule[0].isPhoton(pInTempModule[0].pid())) {
          vStartVecTemp.push_back(vStart);
        }
      } else if (pOutTemp.size() == 1) {
        // no need to generate a vStart for photons and liquefied
        // partons
        if (pOutTemp[0].pstat() != droplet_stat &&
            pOutTemp[0].pstat() != miss_stat &&
            pOutTemp[0].pstat() != neg_stat &&
            !pOutTemp[0].isPhoton(pOutTemp[0].pid())) {
          vStartVecTemp.push_back(vStart);
        }
      } else {
        for (int k = 0; k < pOutTemp.size(); k++) {
          int edgeid = 0;
          if (pOutTemp[k].pstat() == neg_stat) {
            node vNewRootNode = pShower->new_vertex(
                                  make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
            edgeid = pShower->new_parton(vNewRootNode, vStart,
                                         make_shared<Parton>(pOutTemp[k]));
          } else {
            vEnd =
              pShower->new_vertex(make_shared<Vertex>(0, 0, 0, currentTime));
            edgeid = pShower->new_parton(vStart, vEnd,
                                         make_shared<Parton>(pOutTemp[k]));
          }
          pOutTemp[k].set_shower(pShower);
          pOutTemp[k].set_edgeid(edgeid);

          // no need to generate a vStart for photons and liquefied
          // partons
          if (pOutTemp[k].pstat() != droplet_stat &&
              pOutTemp[k].pstat() != miss_stat &&
              pOutTemp[k].pstat() != neg_stat &&
              !pOutTemp[k].isPhoton(pOutTemp[k].pid())) {
            vStartVecOut.push_back(vEnd);
          }

          // --------------------------------------------
          // Add new roots from ElossModules ...
          // (maybe add for clarity a new vector in the signal!???)
          // Otherwise keep track of input size (so far always 1
          // and check if size > 1 and create additional root nodes to that vertex ...
          // Simple Test here below:
          // DEBUG:
          //cout<<"In JetEnergyloss : "<<pInTempModule.size()<<end;
          if (pInTempModule.size() > 1) {
            VERBOSE(7) << pInTempModule.size() - 1
                       << " new root node(s) to be added ...";
            //cout << pInTempModule.size()-1
            //     << " new root node(s) to be added ..." << endl;

            for (int l = 1; l < pInTempModule.size(); l++) {
              node vNewRootNode = pShower->new_vertex(
                                    make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
              pShower->new_parton(vNewRootNode, vEnd,
                                  make_shared<Parton>(pInTempModule[l]));
            }
          }
        }
      }

      // update parton shower
      if (pOutTemp.size() == 0) {
        // this is the free-streaming case for MATTER
        // do not push back droplets
        if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
          if (pInTempModule[0].pstat() == droplet_stat)
            continue;
          if (pInTempModule[0].pstat() == miss_stat)
            continue;
          if (pInTempModule[0].pstat() == neg_stat)
            continue;
        }
        // do not push back photons
        if (pInTempModule[0].isPhoton(pInTempModule[0].pid()))
          continue;
        pInTemp.push_back(pInTempModule[0]);
      } else if (pOutTemp.size() == 1) {
        // this is the free-streaming case for MARTINI or LBT
        // do not push back droplets
        if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
          if (pOutTemp[0].pstat() == droplet_stat)
            continue;
          if (pOutTemp[0].pstat() == miss_stat)
            continue;
          if (pOutTemp[0].pstat() == neg_stat)
            continue;
        }
        // do not push back photons
        if (pOutTemp[0].isPhoton(pOutTemp[0].pid()))
          continue;
        pInTemp.push_back(pOutTemp[0]);
      } else {
        for (int k = 0; k < pOutTemp.size(); k++) {
          // do not push back droplets
          if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
            if (pOutTemp[k].pstat() == droplet_stat)
              continue;
            if (pOutTemp[k].pstat() == miss_stat)
              continue;
            if (pOutTemp[k].pstat() == neg_stat)
              continue;
          }
          // do not push back missing (from AdSCFT)
          if (pOutTemp[k].pstat() == miss_stat)
            continue;
          // do not push back photons
          if (pOutTemp[k].isPhoton(pOutTemp[k].pid()))
            continue;

          pOut.push_back(pOutTemp[k]);
        }
      }
    }

    // one time step is finished, now update parton shower to pIn
    pIn.clear();
    pIn.insert(pIn.end(), pInTemp.begin(), pInTemp.end());
    pIn.insert(pIn.end(), pOut.begin(), pOut.end());

    // update vertex vector
    vStartVec.clear();
    vStartVec.insert(vStartVec.end(), vStartVecTemp.begin(),
                     vStartVecTemp.end());
    vStartVec.insert(vStartVec.end(), vStartVecOut.begin(), vStartVecOut.end());
  } else {
    JSWARN << "NO Initial Hard Parton for Parton shower received ...";
  }
}

void JetEnergyLoss::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  VERBOSE(4) << "In JetEnergyLoss::WriteTask";
  auto f = w.lock();
  if (!f)
    return;

  f->WriteComment("Energy loss Shower Initating Parton: " + GetId());
  f->Write(inP);

  VERBOSE(4) << " writing partons... found " << pShower->GetNumberOfPartons();
  f->Write(pShower);
}

void JetEnergyLoss::PrintShowerInitiatingParton() {
  //JSDEBUG<<inP->pid();
}

void JetEnergyLoss::GetFinalPartonsForEachShower(
  shared_ptr<PartonShower> shower) {

  this->final_Partons.push_back(shower.get()->GetFinalPartons());
}

} // end namespace Jetscape
