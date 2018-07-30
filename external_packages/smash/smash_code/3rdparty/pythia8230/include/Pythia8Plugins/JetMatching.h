// JetMatching.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Richard Corke (implementation of MLM matching as
// in Alpgen for Alpgen input)
// and Stephen Mrenna (implementation of MLM-style matching as
// in Madgraph for Alpgen or Madgraph 5 input.)
// and Simon de Visscher, Stefan Prestel (implementation of shower-kT
// MLM-style matching and flavour treatment for Madgraph input)
// and Stefan Prestel (FxFx NLO jet matching with aMC@NLO.)
// This file provides the classes to perform MLM matching of
// Alpgen or MadGraph 5 input.
// Example usage is shown in main32.cc, and further details
// can be found in the 'Jet Matching Style' manual page.

#ifndef Pythia8_JetMatching_H
#define Pythia8_JetMatching_H

// Includes
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/GeneratorInput.h"

namespace Pythia8 {

//==========================================================================

class HJSlowJet: public SlowJet {

 public:
  HJSlowJet(int powerIn, double Rin, double pTjetMinIn = 0.,
            double etaMaxIn = 25., int selectIn = 1, int massSetIn = 2,
            SlowJetHook* sjHookPtrIn = 0, bool useFJcoreIn = false,
            bool useStandardRin = true) :
  SlowJet(powerIn, Rin, pTjetMinIn, etaMaxIn, selectIn, massSetIn,
          sjHookPtrIn, useFJcoreIn, useStandardRin) {}

 protected:

  void findNext();

};

//--------------------------------------------------------------------------

// Find next cluster pair to join.

void HJSlowJet::findNext() {

  // Find smallest of diB, dij.
  if (clSize > 0) {
    iMin =  0;
    jMin = -1;
    dMin = 1.0/TINY;
    // Remove the possibility of choosing a beam clustering
    for (int i = 1; i < clSize; ++i) {
      for (int j = 0; j < i; ++j) {
        if (dij[i*(i-1)/2 + j] < dMin) {
          iMin = i;
          jMin = j;
          dMin = dij[i*(i-1)/2 + j];
        }
      }
    }

  // If no clusters left then instead default values.
  } else {
    iMin = -1;
    jMin = -1;
    dMin = 0.;
  }

}

//==========================================================================

// Declaration of main JetMatching class to perform MLM matching.
// Note that it is defined with virtual inheritance, so that it can
// be combined with other UserHooks classes, see e.g. main33.cc.

class JetMatching : virtual public UserHooks {

public:

  // Constructor and destructor
 JetMatching() : cellJet(NULL), slowJet(NULL), slowJetHard(NULL),
    hjSlowJet(NULL) {}
 ~JetMatching() {
    if (cellJet) delete cellJet;
    if (slowJet) delete slowJet;
    if (slowJetHard) delete slowJetHard;
    if (hjSlowJet) delete hjSlowJet;
  }

  // Initialisation
  virtual bool initAfterBeams() = 0;

  // Process level vetos
  bool canVetoProcessLevel() { return doMerge; }
  bool doVetoProcessLevel(Event& process) {
    eventProcessOrig = process;
    return false;
  }

  // Parton level vetos (before beam remnants and resonance decays)
  bool canVetoPartonLevelEarly() { return doMerge; }
  bool doVetoPartonLevelEarly(const Event& event);

  // Shower step vetoes (after the first emission, for Shower-kT scheme)
  int  numberVetoStep() {return 1;}
  bool canVetoStep() { return false; }
  bool doVetoStep(int,  int, int, const Event& ) { return false; }

protected:

  // Constants to be changed for debug printout or extra checks.
  static const bool MATCHINGDEBUG, MATCHINGCHECK;

  // Different steps of the matching algorithm.
  virtual void sortIncomingProcess(const Event &)=0;
  virtual void jetAlgorithmInput(const Event &, int)=0;
  virtual void runJetAlgorithm()=0;
  virtual bool matchPartonsToJets(int)=0;
  virtual int  matchPartonsToJetsLight()=0;
  virtual int  matchPartonsToJetsHeavy()=0;

  enum vetoStatus { NONE, LESS_JETS, MORE_JETS, HARD_JET,
                    UNMATCHED_PARTON, INCLUSIVE_VETO };
  enum partonTypes { ID_CHARM=4, ID_BOT=5, ID_TOP=6, ID_LEPMIN=11,
    ID_LEPMAX=16, ID_GLUON=21, ID_PHOTON=22 };

  // Master switch for merging
  bool   doMerge;
  // Switch for merging in the shower-kT scheme. Needed here because
  // the scheme uses different UserHooks functionality.
  bool   doShowerKt;

  // Maximum and current number of jets
  int    nJetMax, nJet;

  // Jet algorithm parameters
  int    jetAlgorithm;
  double eTjetMin, coneRadius, etaJetMax, etaJetMaxAlgo;

  // Internal jet algorithms
  CellJet* cellJet;
  SlowJet* slowJet;
  SlowJet* slowJetHard;
  HJSlowJet* hjSlowJet;

  // SlowJet specific
  int    slowJetPower;

  // Event records to store original incoming process, final-state of the
  // incoming process and what will be passed to the jet algorithm.
  // Not completely necessary to store all steps, but makes tracking the
  // steps of the algorithm a lot easier.
  Event eventProcessOrig, eventProcess, workEventJet;

  // Sort final-state of incoming process into light/heavy jets and 'other'
  vector<int> typeIdx[3];
  set<int>    typeSet[3];

  // Momenta output of jet algorithm (to provide same output regardless of
  // the selected jet algorithm)
  vector<Vec4> jetMomenta;

  // CellJet specific
  int    nEta, nPhi;
  double eTseed, eTthreshold;

  // Merging procedure parameters
  int    jetAllow, jetMatch, exclusiveMode;
  double coneMatchLight, coneRadiusHeavy, coneMatchHeavy;
  bool   exclusive;

  // Store the minimum eT/pT of matched light jets
  double eTpTlightMin;

};

//==========================================================================

// Declaration of main UserHooks class to perform Alpgen matching.

class JetMatchingAlpgen : virtual public JetMatching {

public:

  // Constructor and destructor
  JetMatchingAlpgen() { }
  ~JetMatchingAlpgen() { }

  // Initialisation
  bool initAfterBeams();

private:

  // Different steps of the matching algorithm.
  void sortIncomingProcess(const Event &);
  void jetAlgorithmInput(const Event &, int);
  void runJetAlgorithm();
  bool matchPartonsToJets(int);
  int  matchPartonsToJetsLight();
  int  matchPartonsToJetsHeavy();

  // Sorting utility
  void sortTypeIdx(vector < int > &vecIn);

  // Constants
  static const double GHOSTENERGY, ZEROTHRESHOLD;

};

//==========================================================================

// Declaration of main UserHooks class to perform Madgraph matching.

class JetMatchingMadgraph : virtual public JetMatching {

public:

  // Constructor and destructor
  JetMatchingMadgraph() : slowJetDJR(NULL) { }
  ~JetMatchingMadgraph() { if (slowJetDJR) delete slowJetDJR; }

  // Initialisation
  bool initAfterBeams();

  // Process level vetos
  bool canVetoProcessLevel() { return doMerge; }
  bool doVetoProcessLevel(Event& process);

  // Shower step vetoes (after the first emission, for Shower-kT scheme)
  int  numberVetoStep() {return 1;}
  bool canVetoStep() { return doShowerKt; }
  bool doVetoStep(int,  int, int, const Event& );

  // Jet algorithm to access the jet separations in the cleaned event
  // after showering.
  SlowJet* slowJetDJR;
  // Functions to return the jet clustering scales and number of ME partons.
  // These are useful to investigate the matching systematics.
  vector<double> getDJR() { return DJR;}
  pair<int,int> nMEpartons() { return nMEpartonsSave;}

  // For systematic variations of the jet matching parameters, it is helpful
  // to decouple the jet matching veto from the internal book-keeping. The
  // veto can then be applied in hindsight by an expert user. The functions
  // below return all the necessary information to do this.
  Event getWorkEventJet() { return workEventJetSave; }
  Event getProcessSubset() { return processSubsetSave; }
  bool  getExclusive() { return exclusive; }
  double getPTfirst() { return pTfirstSave; }

protected:

  // Stored values of all inputs necessary to perform the jet matching, as
  // needed when the veto is applied externally.
  Event processSubsetSave;
  Event workEventJetSave;
  double pTfirstSave;

  // Save if code should apply the veto, or simply store the things necessary
  // to perform the veto externally.
  bool performVeto;

  // Different steps of the matching algorithm.
  void sortIncomingProcess(const Event &);
  void jetAlgorithmInput(const Event &, int);
  void runJetAlgorithm();
  bool matchPartonsToJets(int);
  int  matchPartonsToJetsLight();
  int  matchPartonsToJetsHeavy();
  int  matchPartonsToJetsOther();
  bool doShowerKtVeto(double pTfirst);

  // Functions to clear and set the jet clustering scales.
  void clearDJR() { DJR.resize(0);}
  void setDJR( const Event& event);
  // Functions to clear and set the jet clustering scales.
  void clear_nMEpartons() { nMEpartonsSave.first = nMEpartonsSave.second =-1;}
  void set_nMEpartons( const int nOrig, const int nMatch) {
    clear_nMEpartons();
    nMEpartonsSave.first  = nOrig;
    nMEpartonsSave.second = nMatch;
  };

  // Variables.
  vector<int> origTypeIdx[3];
  int    nQmatch;
  double qCut, qCutSq, clFact;
  bool   doFxFx;
  int    nPartonsNow;
  double qCutME, qCutMESq;

  // Vector to store the jet clustering scales.
  vector<double> DJR;
  // Pair of integers giving the number of ME partons read from LHEF and used
  // in the matching (can be different if some partons should not be matched)
  pair<int,int> nMEpartonsSave;

  // Function to get the current number of partons in the Born state, as
  // read from LHE.
  int npNLO();

};

//==========================================================================

// Main implementation of JetMatching class.
// This may be split out to a separate C++ file if desired,
// but currently included here for ease of use.

//--------------------------------------------------------------------------

// Constants to be changed for debug printout or extra checks.
const bool JetMatching::MATCHINGDEBUG = false;
const bool JetMatching::MATCHINGCHECK = false;

//--------------------------------------------------------------------------

// Early parton level veto (before beam remnants and resonance showers)

inline bool JetMatching::doVetoPartonLevelEarly(const Event& event) {

  // 1) Sort the original incoming process. After this step is performed,
  //    the following assignments have been made:
  //    eventProcessOrig - the original incoming process
  //    eventProcess     - the final-state of the incoming process with
  //                       resonance decays removed (and resonances
  //                       themselves now with positive status code)
  //    typeIdx[0/1/2]   - Indices into 'eventProcess' of
  //                       light jets/heavy jets/other
  //    typeSet[0/1/2]   - Indices into 'event' of light jets/heavy jets/other
  //    workEvent        - partons from the hardest subsystem + ISR + FSR only
  sortIncomingProcess(event);

  // For the shower-kT scheme, do not perform any veto here, as any vetoing
  // will already have taken place in doVetoStep.
  if ( doShowerKt ) return false;

  // Debug printout.
  if (MATCHINGDEBUG) {
    // Begin
    cout << endl << "-------- Begin Madgraph Debug --------" << endl;
    // Original incoming process
    cout << endl << "Original incoming process:";
    eventProcessOrig.list();
    // Final-state of original incoming process
    cout << endl << "Final-state incoming process:";
    eventProcess.list();
    // List categories of sorted particles
    for (size_t i = 0; i < typeIdx[0].size(); i++)
      cout << ((i == 0) ? "Light jets: " : ", ")   << setw(3) << typeIdx[0][i];
    if( typeIdx[0].size()== 0 )
      cout << "Light jets: None";

    for (size_t i = 0; i < typeIdx[1].size(); i++)
      cout << ((i == 0) ? "\nHeavy jets: " : ", ") << setw(3) << typeIdx[1][i];
    for (size_t i = 0; i < typeIdx[2].size(); i++)
      cout << ((i == 0) ? "\nOther:      " : ", ") << setw(3) << typeIdx[2][i];
    // Full event at this stage
    cout << endl << endl << "Event:";
    event.list();
    // Work event (partons from hardest subsystem + ISR + FSR)
    cout << endl << "Work event:";
    workEvent.list();
  }

  // 2) Light/heavy jets: iType = 0 (light jets), 1 (heavy jets)
  int iTypeEnd = (typeIdx[2].empty()) ? 2 : 3;
  for (int iType = 0; iType < iTypeEnd; iType++) {

    // 2a) Find particles which will be passed from the jet algorithm.
    //     Input from 'workEvent' and output in 'workEventJet'.
    jetAlgorithmInput(event, iType);

    // Debug printout.
    if (MATCHINGDEBUG) {
      // Jet algorithm event
      cout << endl << "Jet algorithm event (iType = " << iType << "):";
      workEventJet.list();
    }

    // 2b) Run jet algorithm on 'workEventJet'.
    //     Output is stored in jetMomenta.
    runJetAlgorithm();

    // 2c) Match partons to jets and decide if veto is necessary
    if (matchPartonsToJets(iType) == true) {
      // Debug printout.
      if (MATCHINGDEBUG) {
        cout << endl << "Event vetoed" << endl
             << "----------  End MLM Debug  ----------" << endl;
      }
      return true;
    }
  }

  // Debug printout.
  if (MATCHINGDEBUG) {
    cout << endl << "Event accepted" << endl
         << "----------  End MLM Debug  ----------" << endl;
  }

  // If we reached here, then no veto
  return false;

}

//==========================================================================

// Main implementation of Alpgen UserHooks class.
// This may be split out to a separate C++ file if desired,
// but currently included here for ease of use.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// The energy of ghost particles. For technical reasons, this cannot be
// set arbitrarily low, see 'Particle::TINY' in 'Event.cc' for details.
const double JetMatchingAlpgen::GHOSTENERGY   = 1e-15;

// A zero threshold value for double comparisons.
const double JetMatchingAlpgen::ZEROTHRESHOLD = 1e-10;

//--------------------------------------------------------------------------

// Function to sort typeIdx vectors into descending eT/pT order.
// Uses a selection sort, as number of partons generally small
// and so efficiency not a worry.

inline void JetMatchingAlpgen::sortTypeIdx(vector < int > &vecIn) {
  for (size_t i = 0; i < vecIn.size(); i++) {
    size_t jMax = i;
    double vMax = (jetAlgorithm == 1) ?
      eventProcess[vecIn[i]].eT() :
      eventProcess[vecIn[i]].pT();
    for (size_t j = i + 1; j < vecIn.size(); j++) {
      double vNow = (jetAlgorithm == 1)
        ? eventProcess[vecIn[j]].eT() : eventProcess[vecIn[j]].pT();
      if (vNow > vMax) {
        vMax = vNow;
        jMax = j;
      }
    }
    if (jMax != i) swap(vecIn[i], vecIn[jMax]);
  }
}

//--------------------------------------------------------------------------

// Initialisation routine automatically called from Pythia::init().
// Setup all parts needed for the merging.

inline bool JetMatchingAlpgen::initAfterBeams() {

  // Read in parameters
  doMerge         = settingsPtr->flag("JetMatching:merge");
  jetAlgorithm    = settingsPtr->mode("JetMatching:jetAlgorithm");
  nJet            = settingsPtr->mode("JetMatching:nJet");
  nJetMax         = settingsPtr->mode("JetMatching:nJetMax");
  eTjetMin        = settingsPtr->parm("JetMatching:eTjetMin");
  coneRadius      = settingsPtr->parm("JetMatching:coneRadius");
  etaJetMax       = settingsPtr->parm("JetMatching:etaJetMax");
  doShowerKt      = settingsPtr->flag("JetMatching:doShowerKt");

  // Use etaJetMax + coneRadius in input to jet algorithms
  etaJetMaxAlgo   = etaJetMax + coneRadius;

  // CellJet specific
  nEta            = settingsPtr->mode("JetMatching:nEta");
  nPhi            = settingsPtr->mode("JetMatching:nPhi");
  eTseed          = settingsPtr->parm("JetMatching:eTseed");
  eTthreshold     = settingsPtr->parm("JetMatching:eTthreshold");

  // SlowJet specific
  slowJetPower    = settingsPtr->mode("JetMatching:slowJetPower");
  coneMatchLight  = settingsPtr->parm("JetMatching:coneMatchLight");
  coneRadiusHeavy = settingsPtr->parm("JetMatching:coneRadiusHeavy");
  if (coneRadiusHeavy < 0.) coneRadiusHeavy = coneRadius;
  coneMatchHeavy  = settingsPtr->parm("JetMatching:coneMatchHeavy");

  // Matching procedure
  jetAllow        = settingsPtr->mode("JetMatching:jetAllow");
  jetMatch        = settingsPtr->mode("JetMatching:jetMatch");
  exclusiveMode   = settingsPtr->mode("JetMatching:exclusive");

  // If not merging, then done
  if (!doMerge) return true;

  // Exclusive mode; if set to 2, then set based on nJet/nJetMax
  if (exclusiveMode == 2) {

    // No nJet or nJetMax, so default to exclusive mode
    if (nJet < 0 || nJetMax < 0) {
      infoPtr->errorMsg("Warning in JetMatchingAlpgen:init: "
          "missing jet multiplicity information; running in exclusive mode");
      exclusive = true;

    // Inclusive if nJet == nJetMax, exclusive otherwise
    } else {
      exclusive = (nJet == nJetMax) ? false : true;
    }

  // Otherwise, just set as given
  } else {
    exclusive = (exclusiveMode == 0) ? false : true;
  }

  // Initialise chosen jet algorithm. CellJet.
  if (jetAlgorithm == 1) {

    // Extra options for CellJet. nSel = 1 means that all final-state
    // particles are taken and we retain control of what to select.
    // smear/resolution/upperCut are not used and are set to default values.
    int    nSel = 2, smear = 0;
    double resolution = 0.5, upperCut = 2.;
    cellJet = new CellJet(etaJetMaxAlgo, nEta, nPhi, nSel,
                          smear, resolution, upperCut, eTthreshold);

  // SlowJet
  } else if (jetAlgorithm == 2) {
    slowJet = new SlowJet(slowJetPower, coneRadius, eTjetMin, etaJetMaxAlgo);
  }

  // Check the jetMatch parameter; option 2 only works with SlowJet
  if (jetAlgorithm == 1 && jetMatch == 2) {
    infoPtr->errorMsg("Warning in JetMatchingAlpgen:init: "
        "jetMatch = 2 only valid with SlowJet algorithm. "
        "Reverting to jetMatch = 1");
    jetMatch = 1;
  }

  // Setup local event records
  eventProcessOrig.init("(eventProcessOrig)", particleDataPtr);
  eventProcess.init("(eventProcess)", particleDataPtr);
  workEventJet.init("(workEventJet)", particleDataPtr);

  // Print information
  string jetStr  = (jetAlgorithm ==  1) ? "CellJet" :
                   (slowJetPower == -1) ? "anti-kT" :
                   (slowJetPower ==  0) ? "C/A"     :
                   (slowJetPower ==  1) ? "kT"      : "unknown";
  string modeStr = (exclusive)         ? "exclusive" : "inclusive";
  stringstream nJetStr, nJetMaxStr;
  if (nJet >= 0)    nJetStr    << nJet;    else nJetStr    << "unknown";
  if (nJetMax >= 0) nJetMaxStr << nJetMax; else nJetMaxStr << "unknown";
  cout << endl
       << " *-------  MLM matching parameters  -------*" << endl
       << " |  nJet                |  " << setw(14)
       << nJetStr.str() << "  |" << endl
       << " |  nJetMax             |  " << setw(14)
       << nJetMaxStr.str() << "  |" << endl
       << " |  Jet algorithm       |  " << setw(14)
       << jetStr << "  |" << endl
       << " |  eTjetMin            |  " << setw(14)
       << eTjetMin << "  |" << endl
       << " |  coneRadius          |  " << setw(14)
       << coneRadius << "  |" << endl
       << " |  etaJetMax           |  " << setw(14)
       << etaJetMax << "  |" << endl
       << " |  jetAllow            |  " << setw(14)
       << jetAllow << "  |" << endl
       << " |  jetMatch            |  " << setw(14)
       << jetMatch << "  |" << endl
       << " |  coneMatchLight      |  " << setw(14)
       << coneMatchLight << "  |" << endl
       << " |  coneRadiusHeavy     |  " << setw(14)
       << coneRadiusHeavy << "  |" << endl
       << " |  coneMatchHeavy      |  " << setw(14)
       << coneMatchHeavy << "  |" << endl
       << " |  Mode                |  " << setw(14)
       << modeStr << "  |" << endl
       << " *-----------------------------------------*" << endl;

  return true;
}

//--------------------------------------------------------------------------

// Step (1): sort the incoming particles

inline void JetMatchingAlpgen::sortIncomingProcess(const Event &event) {

  // Remove resonance decays from original process and keep only final
  // state. Resonances will have positive status code after this step.
  omitResonanceDecays(eventProcessOrig, true);
  eventProcess = workEvent;

  // Sort original process final state into light/heavy jets and 'other'.
  // Criteria:
  //   1 <= ID <= 5 and massless, or ID == 21 --> light jet (typeIdx[0])
  //   4 <= ID <= 6 and massive               --> heavy jet (typeIdx[1])
  //   All else                               --> other     (typeIdx[2])
  // Note that 'typeIdx' stores indices into 'eventProcess' (after resonance
  // decays are omitted), while 'typeSet' stores indices into the original
  // process record, 'eventProcessOrig', but these indices are also valid
  // in 'event'.
  for (int i = 0; i < 3; i++) {
    typeIdx[i].clear();
    typeSet[i].clear();
  }
  for (int i = 0; i < eventProcess.size(); i++) {
    // Ignore nonfinal and default to 'other'
    if (!eventProcess[i].isFinal()) continue;
    int idx = 2;

    // Light jets
    if (eventProcess[i].id() == ID_GLUON
      || (eventProcess[i].idAbs() <= ID_BOT
      && abs(eventProcess[i].m()) < ZEROTHRESHOLD)) idx = 0;

    // Heavy jets
    else if (eventProcess[i].idAbs() >= ID_CHARM
      && eventProcess[i].idAbs() <= ID_TOP) idx = 1;

    // Store
    typeIdx[idx].push_back(i);
    typeSet[idx].insert(eventProcess[i].daughter1());
  }

  // Extract partons from hardest subsystem + ISR + FSR only into
  // workEvent. Note no resonance showers or MPIs.
  subEvent(event);
}

//--------------------------------------------------------------------------

// Step (2a): pick which particles to pass to the jet algorithm

inline void JetMatchingAlpgen::jetAlgorithmInput(const Event &event,
  int iType) {

  // Take input from 'workEvent' and put output in 'workEventJet'
  workEventJet = workEvent;

  // Loop over particles and decide what to pass to the jet algorithm
  for (int i = 0; i < workEventJet.size(); ++i) {
    if (!workEventJet[i].isFinal()) continue;

    // jetAllow option to disallow certain particle types
    if (jetAllow == 1) {

      // Original AG+Py6 algorithm explicitly excludes tops,
      // leptons and photons.
      int id = workEventJet[i].idAbs();
      if ( (id >= ID_LEPMIN && id <= ID_LEPMAX) || id == ID_TOP
        || id == ID_PHOTON) {
        workEventJet[i].statusNeg();
        continue;
      }
    }

    // Get the index of this particle in original event
    int idx = workEventJet[i].daughter1();

    // Start with particle idx, and afterwards track mothers
    while (true) {

      // Light jets
      if (iType == 0) {

        // Do not include if originates from heavy jet or 'other'
        if (typeSet[1].find(idx) != typeSet[1].end() ||
            typeSet[2].find(idx) != typeSet[2].end()) {
          workEventJet[i].statusNeg();
          break;
        }

        // Made it to start of event record so done
        if (idx == 0) break;
        // Otherwise next mother and continue
        idx = event[idx].mother1();

      // Heavy jets
      } else if (iType == 1) {

        // Only include if originates from heavy jet
        if (typeSet[1].find(idx) != typeSet[1].end()) break;

        // Made it to start of event record with no heavy jet mother,
        // so DO NOT include particle
        if (idx == 0) {
          workEventJet[i].statusNeg();
          break;
        }

        // Otherwise next mother and continue
        idx = event[idx].mother1();

      } // if (iType)
    } // while (true)
  } // for (i)

  // For jetMatch = 2, insert ghost particles corresponding to
  // each hard parton in the original process
  if (jetMatch == 2) {
    for (int i = 0; i < int(typeIdx[iType].size()); i++) {
      // Get y/phi of the parton
      Vec4   pIn = eventProcess[typeIdx[iType][i]].p();
      double y   = pIn.rap();
      double phi = pIn.phi();

      // Create a ghost particle and add to the workEventJet
      double e   = GHOSTENERGY;
      double e2y = exp(2. * y);
      double pz  = e * (e2y - 1.) / (e2y + 1.);
      double pt  = sqrt(e*e - pz*pz);
      double px  = pt * cos(phi);
      double py  = pt * sin(phi);
      workEventJet.append( ID_GLUON, 99, 0, 0, 0, 0, 0, 0, px, py, pz, e);

      // Extra check on reconstructed y/phi values. If many warnings
      // of this type, GHOSTENERGY may be set too low.
      if (MATCHINGCHECK) {
      int lastIdx = workEventJet.size() - 1;
      if (abs(y   - workEventJet[lastIdx].y())   > ZEROTHRESHOLD ||
          abs(phi - workEventJet[lastIdx].phi()) > ZEROTHRESHOLD)
        infoPtr->errorMsg("Warning in JetMatchingAlpgen:jetAlgorithmInput: "
            "ghost particle y/phi mismatch");
      }

    } // for (i)
  } // if (jetMatch == 2)
}

//--------------------------------------------------------------------------

// Step (2b): run jet algorithm and provide common output

inline void JetMatchingAlpgen::runJetAlgorithm() {

  // Run the jet clustering algorithm
  if (jetAlgorithm == 1)
    cellJet->analyze(workEventJet, eTjetMin, coneRadius, eTseed);
  else
    slowJet->analyze(workEventJet);

  // Extract four-momenta of jets with |eta| < etaJetMax and
  // put into jetMomenta. Note that this is done backwards as
  // jets are removed with SlowJet.
  jetMomenta.clear();
  int iJet = (jetAlgorithm == 1) ? cellJet->size() - 1:
                                   slowJet->sizeJet() - 1;
  for (int i = iJet; i > -1; i--) {
    Vec4 jetMom = (jetAlgorithm == 1) ? cellJet->pMassive(i) :
                                        slowJet->p(i);
    double eta = jetMom.eta();

    if (abs(eta) > etaJetMax) {
      if (jetAlgorithm == 2) slowJet->removeJet(i);
      continue;
    }
    jetMomenta.push_back(jetMom);
  }

  // Reverse jetMomenta to restore eT/pT ordering
  reverse(jetMomenta.begin(), jetMomenta.end());
}

//--------------------------------------------------------------------------

// Step (2c): veto decision (returning true vetoes the event)

inline bool JetMatchingAlpgen::matchPartonsToJets(int iType) {

  // Use two different routines for light/heavy jets as
  // different veto conditions and for clarity
  if (iType == 0) return (matchPartonsToJetsLight() > 0);
  else if (iType == 1) return (matchPartonsToJetsHeavy() > 0);
  return true;
}

//--------------------------------------------------------------------------

// Step(2c): light jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto
//   1 = veto as number of jets less than number of partons
//   2 = veto as exclusive mode and number of jets greater than
//       number of partons
//   3 = veto as inclusive mode and there would be an extra jet
//       that is harder than any matched soft jet
//   4 = veto as there is a parton which does not match a jet

inline int JetMatchingAlpgen::matchPartonsToJetsLight() {

  // Always veto if number of jets is less than original number of jets
  if (jetMomenta.size() < typeIdx[0].size()) return LESS_JETS;
  // Veto if in exclusive mode and number of jets bigger than original
  if (exclusive && jetMomenta.size() > typeIdx[0].size()) return MORE_JETS;

  // Sort partons by eT/pT
  sortTypeIdx(typeIdx[0]);

  // Number of hard partons
  int nParton = typeIdx[0].size();

  // Keep track of which jets have been assigned a hard parton
  vector < bool > jetAssigned;
  jetAssigned.assign(jetMomenta.size(), false);

  // Jet matching procedure: (1) deltaR between partons and jets
  if (jetMatch == 1) {

    // Loop over light hard partons and get 4-momentum
    for (int i = 0; i < nParton; i++) {
      Vec4 p1 = eventProcess[typeIdx[0][i]].p();

      // Track which jet has the minimal dR measure with this parton
      int    jMin  = -1;
      double dRmin = 0.;

      // Loop over all jets (skipping those already assigned).
      for (int j = 0; j < int(jetMomenta.size()); j++) {
        if (jetAssigned[j]) continue;

        // DeltaR between parton/jet and store if minimum
        double dR = (jetAlgorithm == 1)
          ? REtaPhi(p1, jetMomenta[j]) : RRapPhi(p1, jetMomenta[j]);
        if (jMin < 0 || dR < dRmin) {
          dRmin = dR;
          jMin  = j;
        }
      } // for (j)

      // Check for jet-parton match
      if (jMin >= 0 && dRmin < coneRadius * coneMatchLight) {

        // If the matched jet is not one of the nParton hardest jets,
        // the extra left over jet would be harder than some of the
        // matched jets. This is disallowed, so veto.
        if (jMin >= nParton) return HARD_JET;

        // Mark jet as assigned.
        jetAssigned[jMin] = true;

      // If no match, then event will be vetoed in all cases
      } else return UNMATCHED_PARTON;

    } // for (i)

  // Jet matching procedure: (2) ghost particles in SlowJet
  } else {

    // Loop over added 'ghost' particles and find if assigned to a jet
    for (int i = workEventJet.size() - nParton;
        i < workEventJet.size(); i++) {
      int jMin = slowJet->jetAssignment(i);

      // Veto if:
      //  1) not one of nParton hardest jets
      //  2) not assigned to a jet
      //  3) jet has already been assigned
      if (jMin >= nParton)               return HARD_JET;
      if (jMin < 0 || jetAssigned[jMin]) return UNMATCHED_PARTON;

      // Mark jet as assigned
      jetAssigned[jMin] = true;

    } // for (i)
  } // if (jetMatch)

  // Minimal eT/pT (CellJet/SlowJet) of matched light jets. Needed
  // later for heavy jet vetos in inclusive mode.
  if (nParton > 0)
    eTpTlightMin = (jetAlgorithm == 1) ? jetMomenta[nParton - 1].eT()
                                       : jetMomenta[nParton - 1].pT();
  else
    eTpTlightMin = -1.;

  // No veto
  return NONE;
}

//--------------------------------------------------------------------------

// Step(2c): heavy jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto as there are no extra jets present
//   1 = veto as in exclusive mode and extra jets present
//   2 = veto as in inclusive mode and extra jets were harder
//       than any matched light jet

inline int JetMatchingAlpgen::matchPartonsToJetsHeavy() {

  // If there are no extra jets, then accept
  if (jetMomenta.empty()) return NONE;

  // Number of hard partons
  int nParton = typeIdx[1].size();

  // Remove jets that are close to heavy quarks
  set < int > removeJets;

  // Jet matching procedure: (1) deltaR between partons and jets
  if (jetMatch == 1) {

    // Loop over heavy hard partons and get 4-momentum
    for (int i = 0; i < nParton; i++) {
      Vec4 p1 = eventProcess[typeIdx[1][i]].p();

      // Loop over all jets, find dR and mark for removal if match
      for (int j = 0; j < int(jetMomenta.size()); j++) {
        double dR = (jetAlgorithm == 1) ?
            REtaPhi(p1, jetMomenta[j]) : RRapPhi(p1, jetMomenta[j]);
        if (dR < coneRadiusHeavy * coneMatchHeavy)
          removeJets.insert(j);

      } // for (j)
    } // for (i)

  // Jet matching procedure: (2) ghost particles in SlowJet
  } else {

    // Loop over added 'ghost' particles and if assigned to a jet
    // then mark this jet for removal
    for (int i = workEventJet.size() - nParton;
        i < workEventJet.size(); i++) {
      int jMin = slowJet->jetAssignment(i);
      if (jMin >= 0) removeJets.insert(jMin);
    }

  }

  // Remove jets (backwards order to not disturb indices)
  for (set < int >::reverse_iterator it  = removeJets.rbegin();
                                     it != removeJets.rend(); it++)
    jetMomenta.erase(jetMomenta.begin() + *it);

  // Handle case if there are still extra jets
  if (!jetMomenta.empty()) {

    // Exclusive mode, so immediate veto
    if (exclusive) return MORE_JETS;

    // Inclusive mode; extra jets must be softer than any matched light jet
    else if (eTpTlightMin >= 0.)
      for (size_t j = 0; j < jetMomenta.size(); j++) {
        // CellJet uses eT, SlowJet uses pT
        if ( (jetAlgorithm == 1 && jetMomenta[j].eT() > eTpTlightMin) ||
             (jetAlgorithm == 2 && jetMomenta[j].pT() > eTpTlightMin) )
          return HARD_JET;
      }

  } // if (!jetMomenta.empty())

  // No extra jets were present so no veto
  return NONE;
}

//==========================================================================

// Main implementation of Madgraph UserHooks class.
// This may be split out to a separate C++ file if desired,
// but currently included here for ease of use.

//--------------------------------------------------------------------------

// Initialisation routine automatically called from Pythia::init().
// Setup all parts needed for the merging.

inline bool JetMatchingMadgraph::initAfterBeams() {

  // Initialise values for stored jet matching veto inputs.
  pTfirstSave = -1.;
  processSubsetSave.init("(eventProcess)", particleDataPtr);
  workEventJetSave.init("(workEventJet)", particleDataPtr);

  // Read in Madgraph specific configuration variables
  bool setMad    = settingsPtr->flag("JetMatching:setMad");

  // If Madgraph parameters are present, then parse in MadgraphPar object
  MadgraphPar par(infoPtr);
  string parStr = infoPtr->header("MGRunCard");
  if (!parStr.empty()) {
    par.parse(parStr);
    par.printParams();
  }

  // Set Madgraph merging parameters from the file if requested
  if (setMad) {
    if ( par.haveParam("xqcut")    && par.haveParam("maxjetflavor")
      && par.haveParam("alpsfact") && par.haveParam("ickkw") ) {
      settingsPtr->flag("JetMatching:merge", par.getParam("ickkw"));
      settingsPtr->parm("JetMatching:qCut", par.getParam("xqcut"));
      settingsPtr->mode("JetMatching:nQmatch",
        par.getParamAsInt("maxjetflavor"));
      settingsPtr->parm("JetMatching:clFact",
        clFact = par.getParam("alpsfact"));
      if (par.getParamAsInt("ickkw") == 0)
        infoPtr->errorMsg("Error in JetMatchingMadgraph:init: "
          "Madgraph file parameters are not set for merging");

    // Warn if setMad requested, but one or more parameters not present
    } else {
       infoPtr->errorMsg("Warning in JetMatchingMadgraph:init: "
          "Madgraph merging parameters not found");
       if (!par.haveParam("xqcut")) infoPtr->errorMsg("Warning in "
          "JetMatchingMadgraph:init: No xqcut");
       if (!par.haveParam("ickkw")) infoPtr->errorMsg("Warning in "
          "JetMatchingMadgraph:init: No ickkw");
       if (!par.haveParam("maxjetflavor")) infoPtr->errorMsg("Warning in "
          "JetMatchingMadgraph:init: No maxjetflavor");
       if (!par.haveParam("alpsfact")) infoPtr->errorMsg("Warning in "
          "JetMatchingMadgraph:init: No alpsfact");
    }
  }

  // Read in FxFx matching parameters
  doFxFx       = settingsPtr->flag("JetMatching:doFxFx");
  nPartonsNow  = settingsPtr->mode("JetMatching:nPartonsNow");
  qCutME       = settingsPtr->parm("JetMatching:qCutME");
  qCutMESq     = pow(qCutME,2);

  // Read in Madgraph merging parameters
  doMerge      = settingsPtr->flag("JetMatching:merge");
  doShowerKt   = settingsPtr->flag("JetMatching:doShowerKt");
  qCut         = settingsPtr->parm("JetMatching:qCut");
  nQmatch      = settingsPtr->mode("JetMatching:nQmatch");
  clFact       = settingsPtr->parm("JetMatching:clFact");

  // Read in jet algorithm parameters
  jetAlgorithm   = settingsPtr->mode("JetMatching:jetAlgorithm");
  nJetMax        = settingsPtr->mode("JetMatching:nJetMax");
  eTjetMin       = settingsPtr->parm("JetMatching:eTjetMin");
  coneRadius     = settingsPtr->parm("JetMatching:coneRadius");
  etaJetMax      = settingsPtr->parm("JetMatching:etaJetMax");
  slowJetPower   = settingsPtr->mode("JetMatching:slowJetPower");

  // Matching procedure
  jetAllow       = settingsPtr->mode("JetMatching:jetAllow");
  exclusiveMode  = settingsPtr->mode("JetMatching:exclusive");
  qCutSq         = pow(qCut,2);
  etaJetMaxAlgo  = etaJetMax;

  // Read if veto should be performed internally.
  performVeto    = settingsPtr->flag("JetMatching:doVeto");

  // If not merging, then done
  if (!doMerge) return true;

  // Exclusive mode; if set to 2, then set based on nJet/nJetMax
  if (exclusiveMode == 2) {

    // No nJet or nJetMax, so default to exclusive mode
    if (nJetMax < 0) {
      infoPtr->errorMsg("Warning in JetMatchingMadgraph:init: "
        "missing jet multiplicity information; running in exclusive mode");
      exclusiveMode = 1;
    }
  }

  // Initialise chosen jet algorithm.
  // Currently, this only supports the kT-algorithm in SlowJet.
  // Use the QCD distance measure by default.
  jetAlgorithm = 2;
  slowJetPower = 1;
  slowJet = new SlowJet(slowJetPower, coneRadius, eTjetMin,
    etaJetMaxAlgo, 2, 2, NULL, false);

  // For FxFx, also initialise jet algorithm to define matrix element jets.
  // Currently, this only supports the kT-algorithm in SlowJet.
  // Use the QCD distance measure by default.
  slowJetHard = new SlowJet(slowJetPower, coneRadius, qCutME,
    etaJetMaxAlgo, 2, 2, NULL, false);

  // To access the DJR's
  slowJetDJR = new SlowJet(slowJetPower, coneRadius, qCutME,
    etaJetMaxAlgo, 2, 2, NULL, false);

  // A special version of SlowJet to handle heavy and other partons
  hjSlowJet = new HJSlowJet(slowJetPower, coneRadius, 0.0,
    100.0, 1, 2, NULL, false, true);

  // Setup local event records
  eventProcessOrig.init("(eventProcessOrig)", particleDataPtr);
  eventProcess.init("(eventProcess)", particleDataPtr);
  workEventJet.init("(workEventJet)", particleDataPtr);

  // Print information
  string jetStr  = (jetAlgorithm ==  1) ? "CellJet" :
                   (slowJetPower == -1) ? "anti-kT" :
                   (slowJetPower ==  0) ? "C/A"     :
                   (slowJetPower ==  1) ? "kT"      : "unknown";
  string modeStr = (exclusiveMode)         ? "exclusive" : "inclusive";
  cout << endl
       << " *-----  Madgraph matching parameters  -----*" << endl
       << " |  qCut                |  " << setw(14)
       << qCut << "  |" << endl
       << " |  nQmatch             |  " << setw(14)
       << nQmatch << "  |" << endl
       << " |  clFact              |  " << setw(14)
       << clFact << "  |" << endl
       << " |  Jet algorithm       |  " << setw(14)
       << jetStr << "  |" << endl
       << " |  eTjetMin            |  " << setw(14)
       << eTjetMin << "  |" << endl
       << " |  etaJetMax           |  " << setw(14)
       << etaJetMax << "  |" << endl
       << " |  jetAllow            |  " << setw(14)
       << jetAllow << "  |" << endl
       << " |  Mode                |  " << setw(14)
       << modeStr << "  |" << endl
       << " *-----------------------------------------*" << endl;

  return true;
}

//--------------------------------------------------------------------------

// Process level vetos

inline bool JetMatchingMadgraph::doVetoProcessLevel(Event& process) {

  eventProcessOrig = process;

  // Setup for veto if hard ME has too many partons.
  // This is done to achieve consistency with the Pythia6 implementation.

  // Clear the event of MPI systems and resonace decay products. Store trimmed
  // event in workEvent.
  sortIncomingProcess(process);

  // Veto in case the hard input matrix element already has too many partons.
  if ( !doFxFx && int(typeIdx[0].size()) > nJetMax )
    return true;
  if ( doFxFx && npNLO() < nJetMax && int(typeIdx[0].size()) > nJetMax )
    return true;

  // Done
  return false;

}

//--------------------------------------------------------------------------

inline bool JetMatchingMadgraph::doVetoStep(int iPos, int nISR, int nFSR,
  const Event& event)  {

  // Do not perform any veto if not in the Shower-kT scheme.
  if ( !doShowerKt ) return false;

  // Do nothing for emissions after the first one.
  if ( nISR + nFSR > 1 ) return false;

  // Do nothing in resonance decay showers.
  if (iPos == 5) return false;

  // Clear the event of MPI systems and resonace decay products. Store trimmed
  // event in workEvent.
  sortIncomingProcess(event);

  // Get (kinematical) pT of first emission
  double pTfirst = 0.;

  // Get weak bosons, for later checks if the emission is a "QCD emission".
  vector<int> weakBosons;
  for (int i = 0; i < event.size(); i++) {
    if ( event[i].id() == 22
      && event[i].id() == 23
      && event[i].idAbs() == 24)
      weakBosons.push_back(i);
  }

  for (int i =  workEvent.size()-1; i > 0; --i) {
    if ( workEvent[i].isFinal() && workEvent[i].colType() != 0
      && (workEvent[i].statusAbs() == 43 || workEvent[i].statusAbs() == 51)) {
      // Check if any of the EW bosons are ancestors of this parton. This
      // should never happen for the first non-resonance shower emission.
      // Check just to be sure.
      bool QCDemission = true;
      // Get position of this parton in the actual event (workEvent does
      // not contain right mother-daughter relations). Stored in daughters.
      int iPosOld = workEvent[i].daughter1();
      for (int j = 0; i < int(weakBosons.size()); ++i)
        if ( event[iPosOld].isAncestor(j)) {
          QCDemission = false;
          break;
        }
      // Done for a QCD emission.
      if (QCDemission){
        pTfirst = workEvent[i].pT();
        break;
      }
    }
  }

  // Store things that are necessary to perform the shower-kT veto externally.
  pTfirstSave   = pTfirst;
  // Done if only inputs for an external vetoing procedure should be stored.
  if (!performVeto) return false;

  // Check veto.
  if ( doShowerKtVeto(pTfirst) ) return true;

  // No veto if come this far.
  return false;

}

//--------------------------------------------------------------------------

inline bool JetMatchingMadgraph::doShowerKtVeto(double pTfirst) {

  // Only check veto in the shower-kT scheme.
  if ( !doShowerKt ) return false;

  // Reset veto code
  bool doVeto = false;

  // Find the (kinematical) pT of the softest (light) parton in the hard
  // process.
  int nParton = typeIdx[0].size();
  double pTminME=1e10;
  for ( int i = 0; i < nParton; ++i)
    pTminME = min(pTminME,eventProcess[typeIdx[0][i]].pT());

  // Veto if the softest hard process parton is below Qcut.
  if ( nParton > 0 && pow(pTminME,2) < qCutSq ) doVeto = true;

  // For non-highest multiplicity, veto if the hardest emission is harder
  // than Qcut.
  if ( exclusive && pow(pTfirst,2) > qCutSq ) {
    doVeto = true;
  // For highest multiplicity sample, veto if the hardest emission is harder
  // than the hard process parton.
  } else if ( !exclusive && nParton > 0 && pTfirst > pTminME ) {
    doVeto = true;
  }

  // Return veto
  return doVeto;

}

//--------------------------------------------------------------------------

// Function to set the jet clustering scales (to be used as output)

inline void JetMatchingMadgraph::setDJR( const Event& event) {

 // Clear members.
 clearDJR();
 vector<double> result;

  // Initialize SlowJetDJR jet algorithm with event
  if (!slowJetDJR->setup(event) ) {
    infoPtr->errorMsg("Warning in JetMatchingMadgraph:setDJR"
      ": the SlowJet algorithm failed on setup");
    return;
  }

  // Cluster in steps to find all hadronic jets
  while ( slowJetDJR->sizeAll() - slowJetDJR->sizeJet() > 0 ) {
    // Save the next clustering scale.
    result.push_back(sqrt(slowJetDJR->dNext()));
    // Perform step.
    slowJetDJR->doStep();
  }

  // Save clustering scales in reserve order.
  for (int i=int(result.size())-1; i >= 0; --i)
    DJR.push_back(result[i]);

}

//--------------------------------------------------------------------------

// Function to get the current number of partons in the Born state, as
// read from LHE.

inline int JetMatchingMadgraph::npNLO(){
  string npIn = infoPtr->getEventAttribute("npNLO",true);
  int np = (npIn != "") ? atoi((char*)npIn.c_str()) : -1;
  if ( np < 0 ) { ; }
  else return np;
  return nPartonsNow;
}

//--------------------------------------------------------------------------

// Step (1): sort the incoming particles

inline void JetMatchingMadgraph::sortIncomingProcess(const Event &event) {

  // Remove resonance decays from original process and keep only final
  // state. Resonances will have positive status code after this step.
  omitResonanceDecays(eventProcessOrig, true);
  clearDJR();
  clear_nMEpartons();

  // For FxFx, pre-cluster partons in the event into jets.
  if (doFxFx) {

    // Get final state partons
    eventProcess.clear();
    workEventJet.clear();
    for( int i=0; i < workEvent.size(); ++i) {
      // Original AG+Py6 algorithm explicitly excludes tops,
      // leptons and photons.
      int id = workEvent[i].idAbs();
      if ((id >= ID_LEPMIN && id <= ID_LEPMAX) || id == ID_TOP
        || id == ID_PHOTON || id == 23 || id == 24 || id == 25) {
        eventProcess.append(workEvent[i]);
      } else {
        workEventJet.append(workEvent[i]);
      }
    }

    // Initialize SlowJetHard jet algorithm with current working event
    if (!slowJetHard->setup(workEventJet) ) {
      infoPtr->errorMsg("Warning in JetMatchingMadgraph:sortIncomingProcess"
        ": the SlowJet algorithm failed on setup");
      return;
    }

    // Get matrix element cut scale.
    double localQcutSq = qCutMESq;
    // Cluster in steps to find all hadronic jets at the scale qCutME
    while ( slowJetHard->sizeAll() - slowJetHard->sizeJet() > 0 ) {
      // Done if next step is above qCut
      if( slowJetHard->dNext() > localQcutSq ) break;
      // Done if we're at or below the number of partons in the Born state.
      if( slowJetHard->sizeAll()-slowJetHard->sizeJet() <= npNLO()) break;
      slowJetHard->doStep();
    }

    // Construct a master copy of the event containing only the
    // hardest nPartonsNow hadronic clusters. While constructing the event,
    // the parton type (ID_GLUON) and status (98,99) are arbitrary.
    int nJets = slowJetHard->sizeJet();
    int nClus = slowJetHard->sizeAll();
    int nNow = 0;
    for (int i = nJets; i < nClus; ++i) {
      vector<int> parts;
      if (i < nClus-nJets) parts = slowJetHard->clusConstituents(i);
      else parts = slowJetHard->constituents(nClus-nJets-i);
      int flavour = ID_GLUON;
      for(int j=0; j < int(parts.size()); ++j)
        if (workEventJet[parts[j]].id() == ID_BOT)
          flavour = ID_BOT;
      eventProcess.append( flavour, 98,
        workEventJet[parts.back()].mother1(),
        workEventJet[parts.back()].mother2(),
        workEventJet[parts.back()].daughter1(),
        workEventJet[parts.back()].daughter2(),
        0, 0, slowJetHard->p(i).px(), slowJetHard->p(i).py(),
        slowJetHard->p(i).pz(), slowJetHard->p(i).e() );
      nNow++;
    }

    // Done. Clean-up
    workEventJet.clear();

  // For MLM matching, simply take hard process state from workEvent,
  // without any preclustering.
  } else {
    eventProcess = workEvent;
  }

  // Sort original process final state into light/heavy jets and 'other'.
  // Criteria:
  //   1 <= ID <= nQmatch, or ID == 21         --> light jet (typeIdx[0])
  //   nQMatch < ID                            --> heavy jet (typeIdx[1])
  //   All else that is colored                --> other     (typeIdx[2])
  // Note that 'typeIdx' stores indices into 'eventProcess' (after resonance
  // decays are omitted), while 'typeSet' stores indices into the original
  // process record, 'eventProcessOrig', but these indices are also valid
  // in 'event'.
  for (int i = 0; i < 3; i++) {
    typeIdx[i].clear();
    typeSet[i].clear();
    origTypeIdx[i].clear();
  }
  for (int i = 0; i < eventProcess.size(); i++) {
    // Ignore non-final state and default to 'other'
    if (!eventProcess[i].isFinal()) continue;
    int idx = -1;
    int orig_idx = -1;

    // Light jets: all gluons and quarks with id less than or equal to nQmatch
    if (eventProcess[i].isGluon()
      || (eventProcess[i].idAbs() <= nQmatch) ) {
      orig_idx = 0;
      // Crucial point: MG puts the scale of a non-QCD particle to eCM. For
      // such particles, we should keep the default "2"
      idx = ( eventProcess[i].scale() < 1.999 * sqrt(infoPtr->eA()
        * infoPtr->eB()) ) ? 0 : 2;
    }

    // Heavy jets:  all quarks with id greater than nQmatch
    else if (eventProcess[i].idAbs() > nQmatch
      && eventProcess[i].idAbs() <= ID_TOP) {
      idx = 1;
      orig_idx = 1;
    // Update to include non-SM colored particles
    } else if (eventProcess[i].colType() != 0
      && eventProcess[i].idAbs() > ID_TOP) {
      idx = 1;
      orig_idx = 1;
    }
    if( idx < 0 ) continue;
    // Store
    typeIdx[idx].push_back(i);
    typeSet[idx].insert(eventProcess[i].daughter1());
    origTypeIdx[orig_idx].push_back(i);
  }

  // Exclusive mode; if set to 2, then set based on nJet/nJetMax
  if (exclusiveMode == 2) {

    // Inclusive if nJet == nJetMax, exclusive otherwise
    int nParton = origTypeIdx[0].size();
    exclusive = (nParton == nJetMax) ? false : true;

  // Otherwise, just set as given
  } else {
    exclusive = (exclusiveMode == 0) ? false : true;
  }

  // Extract partons from hardest subsystem + ISR + FSR only into
  // workEvent. Note no resonance showers or MPIs.
  subEvent(event);

  // Store things that are necessary to perform the kT-MLM veto externally.
  int nParton = typeIdx[0].size();
  processSubsetSave.clear();
  for ( int i = 0; i < nParton; ++i)
    processSubsetSave.append( eventProcess[typeIdx[0][i]] );

}

//--------------------------------------------------------------------------

// Step (2a): pick which particles to pass to the jet algorithm

inline void JetMatchingMadgraph::jetAlgorithmInput(const Event &event,
  int iType) {

  // Take input from 'workEvent' and put output in 'workEventJet'
  workEventJet = workEvent;

  // Loop over particles and decide what to pass to the jet algorithm
  for (int i = 0; i < workEventJet.size(); ++i) {
    if (!workEventJet[i].isFinal()) continue;

    // jetAllow option to disallow certain particle types
    if (jetAllow == 1) {
      // Remove all non-QCD partons from veto list
      if( workEventJet[i].colType() == 0 ) {
        workEventJet[i].statusNeg();
        continue;
      }
    }

    // Get the index of this particle in original event
    int idx = workEventJet[i].daughter1();

    // Start with particle idx, and afterwards track mothers
    while (true) {

      // Light jets
      if (iType == 0) {

        // Do not include if originates from heavy jet or 'other'
        if (typeSet[1].find(idx) != typeSet[1].end() ||
            typeSet[2].find(idx) != typeSet[2].end()) {
          workEventJet[i].statusNeg();
          break;
        }

        // Made it to start of event record so done
        if (idx == 0) break;
        // Otherwise next mother and continue
        idx = event[idx].mother1();

      // Heavy jets
      } else if (iType == 1) {

        // Only include if originates from heavy jet
        if (typeSet[1].find(idx) != typeSet[1].end()) break;

        // Made it to start of event record with no heavy jet mother,
        // so DO NOT include particle
        if (idx == 0) {
          workEventJet[i].statusNeg();
          break;
        }

        // Otherwise next mother and continue
        idx = event[idx].mother1();

      // Other jets
      } else if (iType == 2) {

        // Only include if originates from other jet
        if (typeSet[2].find(idx) != typeSet[2].end()) break;

        // Made it to start of event record with no heavy jet mother,
        // so DO NOT include particle
        if (idx == 0) {
          workEventJet[i].statusNeg();
          break;
        }

        // Otherwise next mother and continue
        idx = event[idx].mother1();

      } // if (iType)
    } // while (true)
  } // for (i)
}

//--------------------------------------------------------------------------

// Step (2b): run jet algorithm and provide common output
// This does nothing, because the jet algorithm is run several times
//  in the matching algorithm.

inline void JetMatchingMadgraph::runJetAlgorithm() {; }

//--------------------------------------------------------------------------

// Step (2c): veto decision (returning true vetoes the event)

inline bool JetMatchingMadgraph::matchPartonsToJets(int iType) {

  // Use different routines for light/heavy/other jets as
  // different veto conditions and for clarity
  if (iType == 0) {
    // Record the jet separations here, also if matchPartonsToJetsLight
    // returns preemptively.
    setDJR(workEventJet);
    set_nMEpartons(origTypeIdx[0].size(), typeIdx[0].size());
    // Perform jet matching.
    return (matchPartonsToJetsLight() > 0);
  } else if (iType == 1) {
     return (matchPartonsToJetsHeavy() > 0);
  } else {
     return (matchPartonsToJetsOther() > 0);
  }

}

//--------------------------------------------------------------------------

// Step(2c): light jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto
//   1 = veto as number of jets less than number of partons
//   2 = veto as exclusive mode and number of jets greater than
//       number of partons
//   3 = veto as inclusive mode and there would be an extra jet
//       that is harder than any matched soft jet
//   4 = veto as there is a parton which does not match a jet

inline int JetMatchingMadgraph::matchPartonsToJetsLight() {

  // Store things that are necessary to perform the kT-MLM veto externally.
  workEventJetSave  = workEventJet;
  // Done if only inputs for an external vetoing procedure should be stored.
  if (!performVeto) return false;

  // Count the number of hard partons
  int nParton = typeIdx[0].size();

  // Initialize SlowJet with current working event
  if (!slowJet->setup(workEventJet) ) {
    infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
      "Light: the SlowJet algorithm failed on setup");
    return NONE;
  }
  double localQcutSq = qCutSq;
  double dOld = 0.0;
  // Cluster in steps to find all hadronic jets at the scale qCut
  while ( slowJet->sizeAll() - slowJet->sizeJet() > 0 ) {
    if( slowJet->dNext() > localQcutSq ) break;
    dOld = slowJet->dNext();
    slowJet->doStep();
  }
  int nJets = slowJet->sizeJet();
  int nClus = slowJet->sizeAll();

  // Debug printout.
  if (MATCHINGDEBUG) slowJet->list(true);

  // Count of the number of hadronic jets in SlowJet accounting
  int nCLjets = nClus - nJets;
  // Get number of partons. Different for MLM and FxFx schemes.
  int nRequested = (doFxFx) ? npNLO() : nParton;

  // Veto event if too few hadronic jets
  if ( nCLjets < nRequested ) return LESS_JETS;

  // In exclusive mode, do not allow more hadronic jets than partons
  if ( exclusive && !doFxFx ) {
    if ( nCLjets > nRequested ) return MORE_JETS;
  } else {

    // For FxFx, in the non-highest multipicity, all jets need to matched to
    // partons. For nCLjets > nRequested, this is not possible. Hence, we can
    // veto here already.
    if ( doFxFx && nRequested < nJetMax && nCLjets > nRequested )
      return MORE_JETS;

    // Now continue in inclusive mode.
    // In inclusive mode, there can be more hadronic jets than partons,
    // provided that all partons are properly matched to hadronic jets.
    // Start by setting up the jet algorithm.
    if (!slowJet->setup(workEventJet) ) {
      infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
        "Light: the SlowJet algorithm failed on setup");
      return NONE;
    }

    // For FxFx, continue clustering as long as the jet separation is above
    // qCut.
    if (doFxFx) {
      while ( slowJet->sizeAll() - slowJet->sizeJet() > 0 ) {
        if( slowJet->dNext() > localQcutSq ) break;
        slowJet->doStep();
      }
    // For MLM, cluster into hadronic jets until there are the same number as
    // partons.
    } else {
      while ( slowJet->sizeAll() - slowJet->sizeJet() > nParton )
        slowJet->doStep();
    }

    // Sort partons in pT.  Update local qCut value.
    //  Hadronic jets are already sorted in pT.
    localQcutSq = dOld;
    if ( clFact >= 0. && nParton > 0 ) {
       vector<double> partonPt;
       for (int i = 0; i < nParton; ++i)
         partonPt.push_back( eventProcess[typeIdx[0][i]].pT2() );
       sort( partonPt.begin(), partonPt.end());
       localQcutSq = max( qCutSq, partonPt[0]);
    }
    nJets = slowJet->sizeJet();
    nClus = slowJet->sizeAll();
  }
  // Update scale if clustering factor is non-zero
  if ( clFact != 0. ) localQcutSq *= pow2(clFact);

  Event tempEvent;
  tempEvent.init( "(tempEvent)", particleDataPtr);
  int nPass = 0;
  double pTminEstimate = -1.;
  // Construct a master copy of the event containing only the
  // hardest nParton hadronic clusters. While constructing the event,
  // the parton type (ID_GLUON) and status (98,99) are arbitrary.
  for (int i = nJets; i < nClus; ++i) {
    tempEvent.append( ID_GLUON, 98, 0, 0, 0, 0, 0, 0, slowJet->p(i).px(),
      slowJet->p(i).py(), slowJet->p(i).pz(), slowJet->p(i).e() );
    ++nPass;
    pTminEstimate = max( pTminEstimate, slowJet->pT(i));
    if(nPass == nRequested) break;
  }

  int tempSize = tempEvent.size();
  // This keeps track of which hadronic jets are matched to parton
  vector<bool> jetAssigned;
  jetAssigned.assign( tempSize, false);

  // This keeps track of which partons are matched to which hadronic
  // jets.
  vector< vector<bool> > partonMatchesJet;
  for (int i=0; i < nParton; ++i )
    partonMatchesJet.push_back( vector<bool>(tempEvent.size(),false) );

  // Begin matching.
  // Do jet matching for FxFx.
  // Make sure that the nPartonsNow hardest hadronic jets are matched to any
  // of the nPartonsNow (+1) partons. This matching is done by attaching a jet
  // from the list of unmatched hadronic jets, and appending a jet from the
  // list of partonic jets, one at a time. The partonic jet will be clustered
  // with the hadronic jet or the beam if the distance measure is below the
  // cut. The hadronic jet is matched once this happens. Otherwise, another
  // partonic jet is tried. When a hadronic jet is matched to a partonic jet,
  // it is removed from the list of unmatched hadronic jets. This process
  // continues until the nPartonsNow hardest hadronic jets are matched to
  // partonic jets, or it is not possible to make a match for a hadronic jet.
  int iNow = 0;
  int nMatched = 0;
  while ( doFxFx && iNow < tempSize ) {

    // Check if this shower jet matches any partonic jet.
    Event tempEventJet;
    tempEventJet.init("(tempEventJet)", particleDataPtr);
    for (int i=0; i < nParton; ++i ) {

      //// Only assign a parton once.
      //for (int j=0; j < tempSize; ++j )
      //  if ( partonMatchesJet[i][j]) continue;

      // Attach a single hadronic jet.
      tempEventJet.clear();
      tempEventJet.append( ID_GLUON, 98, 0, 0, 0, 0, 0, 0,
        tempEvent[iNow].px(), tempEvent[iNow].py(),
        tempEvent[iNow].pz(), tempEvent[iNow].e() );
      // Attach the current parton.
      Vec4 pIn = eventProcess[typeIdx[0][i]].p();
      tempEventJet.append( ID_GLUON, 99, 0, 0, 0, 0, 0, 0,
        pIn.px(), pIn.py(), pIn.pz(), pIn.e() );

      // Setup jet algorithm.
      if ( !slowJet->setup(tempEventJet) ) {
        infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
          "Light: the SlowJet algorithm failed on setup");
        return NONE;
      }

      // These are the conditions for the hadronic jet to match the parton
      //  at the local qCut scale
      if ( slowJet->iNext() == tempEventJet.size() - 1
        && slowJet->jNext() > -1 && slowJet->dNext() < localQcutSq ) {
        jetAssigned[iNow] = true;
        partonMatchesJet[i][iNow] = true;
      }

    } // End loop over hard partons.

    // Veto if the jet could not be assigned to any parton.
    if ( jetAssigned[iNow] ) nMatched++;

    // Continue;
    ++iNow;
  }

  // Jet matching veto for FxFx
  if (doFxFx) {
    if ( nRequested <  nJetMax && nMatched != nRequested )
      return UNMATCHED_PARTON;
    if ( nRequested == nJetMax && nMatched <  nRequested )
      return UNMATCHED_PARTON;
  }

  // Do jet matching for MLM.
  // Take the list of unmatched hadronic jets and append a parton, one at
  // a time. The parton will be clustered with the "closest" hadronic jet
  // or the beam if the distance measure is below the cut. When a hadronic
  // jet is matched to a parton, it is removed from the list of unmatched
  // hadronic jets. This process continues until all hadronic jets are
  // matched to partons or it is not possible to make a match.
  iNow = 0;
  while (!doFxFx && iNow < nParton ) {
    Event tempEventJet;
    tempEventJet.init("(tempEventJet)", particleDataPtr);
    for (int i = 0; i < tempSize; ++i) {
      if (jetAssigned[i]) continue;
      Vec4 pIn = tempEvent[i].p();
      // Append unmatched hadronic jets
      tempEventJet.append( ID_GLUON, 98, 0, 0, 0, 0, 0, 0,
        pIn.px(), pIn.py(), pIn.pz(), pIn.e() );
    }

    Vec4 pIn = eventProcess[typeIdx[0][iNow]].p();
    // Append the current parton
    tempEventJet.append( ID_GLUON, 99, 0, 0, 0, 0, 0, 0,
      pIn.px(), pIn.py(), pIn.pz(), pIn.e() );
    if ( !slowJet->setup(tempEventJet) ) {
      infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
        "Light: the SlowJet algorithm failed on setup");
      return NONE;
    }
    // These are the conditions for the hadronic jet to match the parton
    //  at the local qCut scale
    if ( slowJet->iNext() == tempEventJet.size() - 1
      && slowJet->jNext() > -1 && slowJet->dNext() < localQcutSq ) {
      int iKnt = -1;
      for (int i = 0; i != tempSize; ++i) {
        if (jetAssigned[i]) continue;
        ++iKnt;
        // Identify the hadronic jet that matches the parton
        if (iKnt == slowJet->jNext() ) jetAssigned[i] = true;
      }
    } else {
      return UNMATCHED_PARTON;
    }
    ++iNow;
  }

  // Minimal eT/pT (CellJet/SlowJet) of matched light jets.
  // Needed later for heavy jet vetos in inclusive mode.
  // This information is not used currently.
  if (nParton > 0 && pTminEstimate > 0) eTpTlightMin = pTminEstimate;
  else eTpTlightMin = -1.;

  // Record the jet separations.
  setDJR(workEventJet);

  // No veto
  return NONE;
}

//--------------------------------------------------------------------------

// Step(2c): heavy jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto as there are no extra jets present
//   1 = veto as in exclusive mode and extra jets present
//   2 = veto as in inclusive mode and extra jets were harder
//       than any matched light jet

inline int JetMatchingMadgraph::matchPartonsToJetsHeavy() {

  // Currently, heavy jets are unmatched
  // If there are no extra jets, then accept
  // jetMomenta is NEVER used by MadGraph and is always empty.
  //  This check does nothing.
  //  Rather, if there is any heavy flavor that is harder than
  //  what is present at the LHE level, then the event should
  //  be vetoed.

  // if (jetMomenta.empty()) return NONE;
  // Count the number of hard partons
  int nParton = typeIdx[1].size();

  Event tempEventJet(workEventJet);

  double scaleF(1.0);
  // Rescale the heavy partons that are from the hard process to
  //  have pT=collider energy.   Soft/collinear gluons will cluster
  //  onto them, leaving a remnant of hard emissions.
  for( int i=0; i<nParton; ++i) {
    scaleF = eventProcessOrig[0].e()/workEventJet[typeIdx[1][i]].pT();
    tempEventJet[typeIdx[1][i]].rescale5(scaleF);
  }

  if (!hjSlowJet->setup(tempEventJet) ) {
    infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
                        "Heavy: the SlowJet algorithm failed on setup");
    return NONE;
  }


  while ( hjSlowJet->sizeAll() - hjSlowJet->sizeJet() > 0 ) {
    if( hjSlowJet->dNext() > qCutSq ) break;
    hjSlowJet->doStep();
  }

  int nCLjets(0);
  // Count the number of clusters with pT>qCut.  This includes the
  //  original hard partons plus any hard emissions.
  for(int idx=0 ; idx< hjSlowJet->sizeAll(); ++idx) {
    if( hjSlowJet->pT(idx) > sqrt(qCutSq) ) nCLjets++;
  }

  // Debug printout.
  if (MATCHINGDEBUG) hjSlowJet->list(true);

  // Count of the number of hadronic jets in SlowJet accounting
  //  int nCLjets = nClus - nJets;
  // Get number of partons. Different for MLM and FxFx schemes.
  int nRequested = nParton;

  // Veto event if too few hadronic jets
  if ( nCLjets < nRequested ) {
    if (MATCHINGDEBUG) cout << "veto : hvy  LESS_JETS " << endl;
    if (MATCHINGDEBUG) cout << "nCLjets = " << nCLjets << "; nRequest = "
      << nRequested << endl;
    return LESS_JETS;
  }

  // In exclusive mode, do not allow more hadronic jets than partons
  if ( exclusive ) {
    if ( nCLjets > nRequested ) {
      if (MATCHINGDEBUG) cout << "veto : excl hvy  MORE_JETS " << endl;
      return MORE_JETS;
    }
  }

  // No extra jets were present so no veto
  return NONE;
}

//--------------------------------------------------------------------------

// Step(2c): other jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto as there are no extra jets present
//   1 = veto as in exclusive mode and extra jets present
//   2 = veto as in inclusive mode and extra jets were harder
//       than any matched light jet

inline int JetMatchingMadgraph::matchPartonsToJetsOther() {

  // Currently, heavy jets are unmatched
  // If there are no extra jets, then accept
  // jetMomenta is NEVER used by MadGraph and is always empty.
  //  This check does nothing.
  //  Rather, if there is any heavy flavor that is harder than
  //  what is present at the LHE level, then the event should
  //  be vetoed.

  // if (jetMomenta.empty()) return NONE;
  // Count the number of hard partons
  int nParton = typeIdx[2].size();

  Event tempEventJet(workEventJet);

  double scaleF(1.0);
  // Rescale the heavy partons that are from the hard process to
  //  have pT=collider energy.   Soft/collinear gluons will cluster
  //  onto them, leaving a remnant of hard emissions.
  for( int i=0; i<nParton; ++i) {
    scaleF = eventProcessOrig[0].e()/workEventJet[typeIdx[2][i]].pT();
    tempEventJet[typeIdx[2][i]].rescale5(scaleF);
  }

  if (!hjSlowJet->setup(tempEventJet) ) {
    infoPtr->errorMsg("Warning in JetMatchingMadgraph:matchPartonsToJets"
                        "Heavy: the SlowJet algorithm failed on setup");
    return NONE;
  }


  while ( hjSlowJet->sizeAll() - hjSlowJet->sizeJet() > 0 ) {
    if( hjSlowJet->dNext() > qCutSq ) break;
    hjSlowJet->doStep();
  }

  int nCLjets(0);
  // Count the number of clusters with pT>qCut.  This includes the
  //  original hard partons plus any hard emissions.
  for(int idx=0 ; idx< hjSlowJet->sizeAll(); ++idx) {
    if( hjSlowJet->pT(idx) > sqrt(qCutSq) ) nCLjets++;
  }

  // Debug printout.
  if (MATCHINGDEBUG) hjSlowJet->list(true);

  // Count of the number of hadronic jets in SlowJet accounting
  //  int nCLjets = nClus - nJets;
  // Get number of partons. Different for MLM and FxFx schemes.
  int nRequested = nParton;

  // Veto event if too few hadronic jets
  if ( nCLjets < nRequested ) {
    if (MATCHINGDEBUG) cout << "veto : other LESS_JETS " << endl;
    if (MATCHINGDEBUG) cout << "nCLjets = " << nCLjets << "; nRequest = "
      << nRequested << endl;
    return LESS_JETS;
  }

  // In exclusive mode, do not allow more hadronic jets than partons
  if ( exclusive ) {
    if ( nCLjets > nRequested ) {
      if (MATCHINGDEBUG) cout << "veto : excl other MORE_JETS" << endl;
      return MORE_JETS;
    }
  }

  // No extra jets were present so no veto
  return NONE;
}

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_JetMatching_H
