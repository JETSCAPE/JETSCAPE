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

// Create a pythia collision at a specified point and return the two inital hard partons

#include "InitialStateRadiationTest.h"
#include <sstream>
#include "QueryHistory.h"

#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<InitialStateRadiationTest> InitialStateRadiationTest::reg("InitialStateRadiationTest");

InitialStateRadiationTest::~InitialStateRadiationTest() { VERBOSE(8); }

void InitialStateRadiationTest::InitTask() {

  //VERBOSE(8);

  std::string s = GetXMLElementText({"Hard", "InitialStateRadiationTest", "name"});
  SetId(s);
  JSDEBUG << "Initializing " << s;
  cout << "Initializing InitialStateRadiationTest" << endl;;
}

void InitialStateRadiationTest::Exec() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();

  // dummy incoming initial hard-scattering partons to be evolve backward in time
  vector<FourVector> p_init;
  p_init.push_back(FourVector(1., 0., 19., 20.));
  p_init.push_back(FourVector(-1., 0., -19., 20.));
  FourVector x_init(0., 0., 0., 0.);
  FourVector x_next(0., 0., 0., -deltaT);

  // create two dummy initial hard-scattering partons
  // and push back to master graph structure
  for (unsigned int i=0; i<2; i++){
    // create an initial hard parton
    Parton p_hard_scat(0, 21, 0, p_init[i], x_init);

    // create a graph structure for each initial parton
    shared_ptr<PartonShower> initShower = make_shared<PartonShower>();

    // create starting and ending node (vertex) for the parton
    // vEnd is one time-step earlier than vStart
    vEnd = initShower->new_vertex(make_shared<Vertex>(x_init));
    vStart = initShower->new_vertex(make_shared<Vertex>(x_next));

    // create a edge that connects the two nodes using the hard parton
    // in the graph structure
    initShower->new_parton(vStart, vEnd, make_shared<Parton>(p_hard_scat));

    // push back the edge to the vector of the master graph structure
    pShowerMaster.push_back(initShower);
  }

  timeVec.resize(pShowerMaster.size());

  // Backward Shower
  BackwardISR();
  // Forward Shower
  ForwardISR();

  int ab=0;
  JSDEBUG << "Add fianl partons";
  // Send final state partons to the framework
  for (auto pShower : pShowerMaster) {
    JSDEBUG << "Shower:" << ab;
    ab++;
    unsigned int NumberOfPartons = pShower->GetNumberOfPartons();
    for (unsigned int ipart=0; ipart < NumberOfPartons; ipart++) {
      if (pShower->GetNumberOfChilds(ipart) == 0 &&
          pShower->GetPartonAt(ipart)->edgeid() > 0) {
        AddParton(pShower->GetPartonAt(ipart));
      }
    }
  }
  JSDEBUG << "End of InitialStateRadiationTest";

  VERBOSE(8) << "GetNHardPartons():" << GetNHardPartons();
}

void InitialStateRadiationTest::BackwardISR() {
  JSDEBUG << "Beginning of Backward Shower";

  vector<Parton> pIn, pOut;
  vector<node> vEndVec, vEndVecTemp;
  int edgeid;

  map<node, Parton> nodePartonPair;

  // iterate over master graph structure
  for (unsigned int ip=0; ip<pShowerMaster.size(); ip++) {

    shared_ptr<PartonShower> pShower = pShowerMaster[ip];
    unsigned int n_parton = pShower->GetNumberOfPartons();

    // time when hard scattering occurs
    currentTime = 0;

    pIn.clear();
    vEndVec.clear();

    // put the initial parton into pIn
    pIn.push_back(*pShower->GetPartonAt(0));
    vEndVec.push_back(pShower->GetNodeAt(1));
    JSDEBUG << "pShowerMaster[" << ip << "]";

    nodePartonPair.clear();

    // dummy iteration over time-steps in backward direction
    for (unsigned int i_timestep=0; i_timestep<n_timeStep; i_timestep++) {

      currentTime -= deltaT;
      JSDEBUG << "currentTime:" << currentTime;


      pOut.clear();
      vEndVecTemp.clear();

      // iterate over pIn -- the size of pIn is always 1?
      for (unsigned int i=0; i<pIn.size(); i++) {

        // create dummy space-like and time-like partons
        FourVector p_val = pIn[i].p_in();
        FourVector p_new(0.1*p_val.x(), 0.1*p_val.y(), 0.1*p_val.z(), 0.1*p_val.t());
        p_val += p_new;
        FourVector x_new(0., 0., 0., currentTime);

        // time-like parton
        Parton p_tlike = Parton(0, 21, timeLike_stat, p_new, x_new);
        vEnd = vEndVec[i];
        node vNewChildNode = pShower->new_vertex(
                              make_shared<Vertex>(0, 0, 0, currentTime + deltaT));
        edgeid = pShower->new_parton(vEnd, vNewChildNode,
                                     make_shared<Parton>(p_tlike));
        JSDEBUG << "time-like vEnd->vNewChildNode:" << vEnd << " " << vNewChildNode
             << " edgeid:" << edgeid;
        n_parton++;
        pShower->GetPartonAt(n_parton-1)->set_edgeid(edgeid);
        pShower->GetPartonAt(n_parton-1)->set_shower(pShower);
        nodePartonPair.insert({vNewChildNode, p_tlike});

        // space-like parton
        Parton p_slike = Parton(0, 21, spaceLike_stat, p_val, x_new);
        vStart = pShower->new_vertex(
                              make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
        edgeid = pShower->new_parton(vStart, vEnd,
                                     make_shared<Parton>(p_slike));
        JSDEBUG << "space-like vStart->vEnd:" << vStart << " " << vEnd
             << " edgeid:" << edgeid;
        pOut.push_back(p_slike);
        vEndVecTemp.push_back(vStart);
        n_parton++;
        pShower->GetPartonAt(n_parton-1)->set_edgeid(edgeid);
        pShower->GetPartonAt(n_parton-1)->set_shower(pShower);
      }

      // update pIn and vEndVec for next time step -- deep copy
      pIn.clear();
      pIn = pOut;
      vEndVec.clear();
      vEndVec = vEndVecTemp;

      // update the earliest time for a given parton shower
      timeVec[ip] = currentTime - deltaT;
    }

    // push back node-parton (time-like) pairs for forward shower
    nodePartonPairVec.push_back(nodePartonPair);
  }
}

void InitialStateRadiationTest::ForwardISR() {

  JSDEBUG << "Beginning of Forward Shower";
  vector<Parton> pIn, pOut;
  vector<node> vStartVec, vStartVecTemp;

  int edgeid;

  // iterate over master graph structure
  for (unsigned int ip=0; ip<pShowerMaster.size(); ip++) {

    shared_ptr<PartonShower> pShower = pShowerMaster[ip];
    unsigned int n_parton = pShower->GetNumberOfPartons();

    pIn.clear();
    vStartVec.clear();

    // put the time-like partons into pIn
    for (const auto& p : nodePartonPairVec[ip]) {
      vStartVec.push_back(p.first);
      pIn.push_back(p.second);
    }

    // current time
    currentTime = timeVec[ip];
    JSDEBUG << "pShowerMaster[" << ip << "]";

    while (currentTime < -deltaT-eps) {

      currentTime += deltaT;
      JSDEBUG << "currentTime:" << currentTime;

      pOut.clear();
      vStartVecTemp.clear();

      for (unsigned int i=0; i<pIn.size(); i++) {

        double parton_time = pIn[i].time();
        JSDEBUG << "parton_time:" << parton_time << " "
             << vStartVec[i];

        // skip if this parton was created later than current time
        if (parton_time > currentTime + eps) {
          vStartVecTemp.push_back(vStartVec[i]);
          pOut.push_back(pIn[i]);
        } else {

          vector<Parton> pInModule, pOutModule;
          pInModule.push_back(pIn[i]);

          // call MATTER -- add later
          //SentInPartons(moduleDeltaT, currentTime, parton.pt(), 
          //              pInModule, pOutModule);

          // instead do dummy forward shower
          FourVector p_val = pIn[i].p_in();
          FourVector p_new(0.5*p_val.x(), 0.5*p_val.y(), 0.5*p_val.z(), 0.5*p_val.t());
          FourVector x_new(0., 0., 0., currentTime);
          pOutModule.push_back(Parton(0, 21, 0, p_new, x_new));
          pOutModule.push_back(Parton(0, 21, 0, p_new, x_new));
          if (i == 0) {
          pOutModule.push_back(Parton(0, 21, -17, p_new, x_new));
          }

          // apply liquefier -- add later
          //if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
          //  liquefier_ptr.lock()->add_hydro_sources(pInModule, pOutModule);
          //}

          vStart = vStartVec[i];

          if (pOutModule.size() == 0) {
            // no need to generate a vStart for photons and liquefied
            // partons
            //if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
              int pstat = pInModule[0].pstat();
              if (pstat == droplet_stat || pstat == miss_stat || pstat == neg_stat)
                continue;
            //}
            if (pInModule[0].isPhoton(pInModule[0].pid()))
              continue;

            vStartVecTemp.push_back(vStart);
            pOut.push_back(pInModule[0]);

          } else if (pOutModule.size() == 1) {
            // no need to generate a vStart for photons and liquefied
            // partons
            //if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
              int pstat = pOutModule[0].pstat();
              if (pstat == droplet_stat || pstat == miss_stat || pstat == neg_stat)
                continue;
            //}
            if (pOutModule[0].isPhoton(pOutModule[0].pid()))
              continue;

            vStartVecTemp.push_back(vStart);
            pOut.push_back(pOutModule[0]);

          } else {
            for (int k = 0; k < pOutModule.size(); k++) {
              int edgeid = 0;
              if (pOutModule[k].pstat() == neg_stat) {
                node vNewRootNode = pShower->new_vertex(
                                  make_shared<Vertex>(0, 0, 0, currentTime - deltaT));
                edgeid = pShower->new_parton(vNewRootNode, vStart,
                                             make_shared<Parton>(pOutModule[k]));
                JSDEBUG << "negative vNewRootNode->vStart:"
                     << vNewRootNode << " " << vStart << " edgeid:" << edgeid;
              } else {
                vEnd = pShower->new_vertex(
                                  make_shared<Vertex>(0, 0, 0, currentTime + deltaT));
                edgeid = pShower->new_parton(vStart, vEnd,
                                             make_shared<Parton>(pOutModule[k]));
                JSDEBUG << "positive vStart->vEnd:"
                     << vStart << " " << vEnd << " edgeid:" << edgeid;
              }
              n_parton++;
              pShower->GetPartonAt(n_parton-1)->set_edgeid(edgeid);
              pShower->GetPartonAt(n_parton-1)->set_shower(pShower);

              // no need to generate a vStart for photons and liquefied
              // partons
              //if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
                int pstat = pOutModule[k].pstat();
                if (pstat == droplet_stat || pstat == miss_stat || pstat == neg_stat ||
                    pstat == -1)
                  continue;
              //}
              if (pOutModule[0].isPhoton(pOutModule[0].pid()))
              continue;

              vStartVecTemp.push_back(vEnd);
              pOut.push_back(pOutModule[k]);
            }
          }
        }
      }

      // update pIn and vStartVec for next time step -- deep copy
      pIn.clear();
      pIn = pOut;
      vStartVec.clear();
      vStartVec = vStartVecTemp;
    }
  }
}
