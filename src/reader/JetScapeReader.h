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

#ifndef JETSCAPEREADER_H
#define JETSCAPEREADER_H

#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include <fstream>

#include "GTL/graph.h"
#include "JetClass.h"
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "PartonShower.h"
#include "StringTokenizer.h"
#ifdef USE_GZIP
#include "gzstream.h"
#endif

using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;

namespace Jetscape {

template <class T>
class JetScapeReader {
 public:
  JetScapeReader();
  JetScapeReader(string m_file_name_in) {
    file_name_in = m_file_name_in;
    Init();
  }
  virtual ~JetScapeReader();

  void Close() { inFile.close(); }
  void Clear();

  void Next();
  bool Finished() { return inFile.eof(); }

  int GetCurrentEvent() { return currentEvent - 1; }
  int GetCurrentNumberOfPartonShowers() { return pShowers.size(); }

  // shared_ptr<PartonShower> GetPartonShower() {return pShower;}
  vector<shared_ptr<PartonShower>> GetPartonShowers() { return pShowers; }

  vector<shared_ptr<Hadron>> GetHadrons() { return hadrons; }
  vector<fjcore::PseudoJet> GetHadronsForFastJet();
  double GetSigmaGen() const { return sigmaGen; }
  double GetSigmaErr() const { return sigmaErr; }
  double GetEventWeight() const { return eventWeight; }
  double GetEventPlaneAngle() const { return EventPlaneAngle; }

 private:
  StringTokenizer strT;

  void Init();
  void AddNode(string s);
  void AddEdge(string s);
  // void MakeGraph();
  void AddHadron(string s);
  string file_name_in;
  T inFile;

  int currentEvent;
  int currentShower;

  shared_ptr<PartonShower> pShower;
  vector<shared_ptr<PartonShower>> pShowers;

  vector<node> nodeVec;
  vector<edge> edgeVec;
  vector<shared_ptr<Hadron>> hadrons;
  double sigmaGen;
  double sigmaErr;
  double eventWeight;
  double EventPlaneAngle;
};

typedef JetScapeReader<ifstream> JetScapeReaderAscii;
#ifdef USE_GZIP
typedef JetScapeReader<igzstream> JetScapeReaderAsciiGZ;
#endif

}  // end namespace Jetscape

// ---------------------

#endif

// ---------------------
