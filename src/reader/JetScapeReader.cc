/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "JetScapeReader.h"

#include <sstream>

namespace Jetscape {

template <class T>
JetScapeReader<T>::JetScapeReader()
    : currentEvent{-1},
      sigmaGen{-1},
      sigmaErr{-1},
      eventWeight{-1},
      EventPlaneAngle{0.0} {
  VERBOSE(8);
}

template <class T>
JetScapeReader<T>::~JetScapeReader() {
  VERBOSE(8);
}

template <class T>
void JetScapeReader<T>::Clear() {
  nodeVec.clear();
  edgeVec.clear();
  // pShower->clear();//pShower=nullptr; //check ...
  pShowers.clear();
  hadrons.clear();

  sigmaGen = -1;
  sigmaErr = -1;
  eventWeight = -1;
  EventPlaneAngle = 0.0;
}

template <class T>
void JetScapeReader<T>::AddNode(string s) {
  string token;
  // int counter=0;
  strT.set(s);

  vector<string> vS;

  while (!strT.done()) {
    token = strT.next();
    if (token.compare("V") != 0)
      vS.push_back(token);
  }

  nodeVec.push_back(pShower->new_vertex(
      make_shared<Vertex>(stod(vS[1]), stod(vS[2]), stod(vS[3]), stod(vS[4]))));
}

template <class T>
void JetScapeReader<T>::AddEdge(string s) {
  if (nodeVec.size() > 1) {
    string token;
    // int counter=0;
    strT.set(s);

    vector<string> vS;

    while (!strT.done()) {
      token = strT.next();
      if (token.compare("P") != 0)
        vS.push_back(token);
    }

    pShower->new_parton(
        nodeVec[stoi(vS[0])], nodeVec[stoi(vS[1])],
        make_shared<Parton>(stoi(vS[2]), stoi(vS[3]), stoi(vS[4]), stod(vS[5]),
                            stod(vS[6]), stod(vS[7]),
                            stod(vS[8])));  // use different constructor wit
                                            // true spatial posiiton ...
  } else
    JSWARN << "Node vector not filled, can not add edges/partons!";
}

template <class T>
void JetScapeReader<T>::AddHadron(string s) {
  string token;
  strT.set(s);

  vector<string> vS;
  double x[4];
  x[0] = x[1] = x[2] = x[3] = 0.0;
  while (!strT.done()) {
    token = strT.next();
    if (token.compare("H") != 0)
      vS.push_back(token);
  }
  hadrons.push_back(make_shared<Hadron>(stoi(vS[1]), stoi(vS[2]), stoi(vS[3]),
                                        stod(vS[4]), stod(vS[5]), stod(vS[6]),
                                        stod(vS[7]), x));
}

template <class T>
void JetScapeReader<T>::Next() {
  if (currentEvent > 0)
    Clear();

  // ReadEvent(currentPos);
  string line;
  string token;

  JSINFO << "Current Event = " << currentEvent;

  pShowers.push_back(make_shared<PartonShower>());
  pShower = pShowers[0];
  currentShower = 1;

  int nodeZeroCounter = 0;
  std::string EPAngleStr = "EventPlaneAngle";
  while (getline(inFile, line)) {
    strT.set(line);

    if (strT.isCommentEntry()) {
      // Cross section
      if (line.find("sigmaGen") != std::string::npos) {
        std::stringstream data(line);
        std::string dummy;
        data >> dummy >> dummy >> dummy >> sigmaGen;
        JSDEBUG << " sigma gen=" << sigmaGen;
      }
      // Cross section error
      if (line.find("sigmaErr") != std::string::npos) {
        std::stringstream data(line);
        std::string dummy;
        data >> dummy >> dummy >> dummy >> sigmaErr;
        JSDEBUG << " sigma err=" << sigmaErr;
      }
      // Event weight
      if (line.find("weight") != std::string::npos) {
        std::stringstream data(line);
        std::string dummy;
        data >> dummy >> dummy >> dummy >> eventWeight;
        JSDEBUG << " Event weight=" << eventWeight;
      }
      // EP angle
      if (line.find(EPAngleStr) != std::string::npos) {
        std::stringstream data(line);
        std::string dummy;
        data >> dummy >> dummy >> dummy >> EventPlaneAngle;
        JSDEBUG << " EventPlaneAngle=" << EventPlaneAngle;
      }
      continue;
    }

    if (strT.isEventEntry()) {
      int newEvent = stoi(strT.next());
      if (currentEvent != newEvent && currentEvent > -1) {
        currentEvent++;
        break;
      }
      currentEvent = newEvent;
      continue;
    }

    // not an event header -- done?
    if (!strT.isGraphEntry())
      continue;

    // node?
    if (strT.isNodeEntry()) {
      // catch starting node
      if (strT.isNodeZero()) {
        nodeZeroCounter++;
        if (nodeZeroCounter > currentShower) {
          nodeVec.clear();
          edgeVec.clear();
          pShowers.push_back(make_shared<PartonShower>());
          pShower = pShowers.back();
          currentShower++;
        }
      }
      AddNode(line);
      continue;
    }

    // edge?
    if (strT.isEdgeEntry()) {
      AddEdge(line);
      continue;
    }

    // rest is list entry == hadron entry
    // Some questionable nomenclature here - identifying all that begins with
    // "[" as "GraphEntry" oh well
    AddHadron(line);
  }

  if (Finished())
    currentEvent++;
}

template <class T>
vector<fjcore::PseudoJet> JetScapeReader<T>::GetHadronsForFastJet() {
  vector<fjcore::PseudoJet> forFJ;

  for (auto &h : hadrons) {
    forFJ.push_back(h->GetPseudoJet());
  }

  return forFJ;
}

template <class T>
void JetScapeReader<T>::Init() {
  VERBOSE(8) << "Open Input File = " << file_name_in;
  JSINFO << "Open Input File = " << file_name_in;

  inFile.open(file_name_in.c_str());

  if (!inFile.good()) {
    JSWARN << "Corrupt input file!";
    exit(-1);
  } else
    JSINFO << "File opened";

  currentEvent = 0;
}

template class JetScapeReader<ifstream>;

#ifdef USE_GZIP
template class JetScapeReader<igzstream>;
#endif

}  // end namespace Jetscape
