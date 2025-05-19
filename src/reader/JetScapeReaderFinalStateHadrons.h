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

#ifndef JETSCAPEREADERFINALSTATEHADRONS_H
#define JETSCAPEREADERFINALSTATEHADRONS_H

#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
#include "JetScapeParticles.h"
#include "JetScapeLogger.h"
#include "StringTokenizer.h"
#include "PartonShower.h"
#include <fstream>

using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;

namespace Jetscape {

class JetScapeReaderFinalStateHadrons {

public:
  JetScapeReaderFinalStateHadrons();
  JetScapeReaderFinalStateHadrons(string m_file_name_in) {
    file_name_in = m_file_name_in;
    InitTask();
  }
  virtual ~JetScapeReaderFinalStateHadrons();

  void Close() { inFile.close(); }
  void ClearTask();

  void Next();
  bool Finished() { return inFile.eof(); }

  int GetCurrentEvent() { return currentEvent - 1; }

  vector<shared_ptr<Hadron>> GetHadrons() { return hadrons; }
  vector<fjcore::PseudoJet> GetHadronsForFastJet();
  double GetSigmaGen() const { return sigmaGen; }
  double GetSigmaErr() const { return sigmaErr; }
  double GetEventWeight() const { return eventWeight; }
  double GetEventPlaneAngle() const { return EventPlaneAngle; }
  int TotalEventCount();

private:
  StringTokenizer strT;

  void InitTask();
  //void MakeGraph();
  void AddHadron(string s);
  string file_name_in;
  ifstream inFile;

  int currentEvent;
  vector<shared_ptr<Hadron>> hadrons;
  double sigmaGen;
  double sigmaErr;
  double eventWeight;
  double EventPlaneAngle;
};

} // end namespace Jetscape

// ---------------------

#endif

// ---------------------
