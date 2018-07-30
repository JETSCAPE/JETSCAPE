// TimeShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the timelike final-state showers.
// TimeDipoleEnd: data on a radiating dipole end.
// TimeShower: handles the showering description.

#ifndef Pythia8_TimeShower_H
#define Pythia8_TimeShower_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/WeakShowerMEs.h"

namespace Pythia8 {

//==========================================================================

// Data on radiating dipole ends; only used inside TimeShower class.

class TimeDipoleEnd {

public:

  // Constructors.
  TimeDipoleEnd() : iRadiator(-1), iRecoiler(-1), pTmax(0.), colType(0),
    chgType(0), gamType(0), weakType(0), isrType(0), system(0), systemRec(0),
    MEtype(0), iMEpartner(-1), weakPol(0), isOctetOnium(false),
    isHiddenValley(false), colvType(0), MEmix(0.), MEorder(true),
    MEsplit(true), MEgluinoRec(false), isFlexible(false) { }
 TimeDipoleEnd(int iRadiatorIn, int iRecoilerIn, double pTmaxIn = 0.,
    int colIn = 0, int chgIn = 0, int gamIn = 0, int weakTypeIn = 0,
    int isrIn = 0, int systemIn = 0, int MEtypeIn = 0, int iMEpartnerIn = -1,
    int weakPolIn = 0, bool isOctetOniumIn = false,
    bool isHiddenValleyIn = false, int colvTypeIn = 0, double MEmixIn = 0.,
    bool MEorderIn = true, bool MEsplitIn = true, bool MEgluinoRecIn = false,
    bool isFlexibleIn = false) :
    iRadiator(iRadiatorIn), iRecoiler(iRecoilerIn), pTmax(pTmaxIn),
    colType(colIn), chgType(chgIn), gamType(gamIn), weakType(weakTypeIn),
    isrType(isrIn), system(systemIn), systemRec(systemIn), MEtype(MEtypeIn),
    iMEpartner(iMEpartnerIn), weakPol(weakPolIn), isOctetOnium(isOctetOniumIn),
    isHiddenValley(isHiddenValleyIn), colvType(colvTypeIn), MEmix(MEmixIn),
    MEorder (MEorderIn), MEsplit(MEsplitIn), MEgluinoRec(MEgluinoRecIn),
    isFlexible(isFlexibleIn)  { }

  // Basic properties related to dipole and matrix element corrections.
  int    iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, gamType, weakType, isrType, system, systemRec,
         MEtype, iMEpartner, weakPol;
  bool   isOctetOnium, isHiddenValley;
  int    colvType;
  double MEmix;
  bool   MEorder, MEsplit, MEgluinoRec, isFlexible;

  // Properties specific to current trial emission.
  int    flavour, iAunt;
  double mRad, m2Rad, mRec, m2Rec, mDip, m2Dip, m2DipCorr,
         pT2, m2, z, mFlavour, asymPol, flexFactor, pAccept;

};

//==========================================================================

// The TimeShower class does timelike showers.

class TimeShower {

public:

  // Constructor.
  TimeShower() {beamOffset = 0;}

  // Destructor.
  virtual ~TimeShower() {}

  // Initialize various pointers.
  // (Separated from rest of init since not virtual.)
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
    CoupSM* coupSMPtrIn, PartonSystems* partonSystemsPtrIn,
    UserHooks* userHooksPtrIn, MergingHooks* mergingHooksPtrIn,
    PartonVertex* partonVertexPtrIn) { infoPtr = infoPtrIn;
    settingsPtr = settingsPtrIn; particleDataPtr = particleDataPtrIn;
    rndmPtr = rndmPtrIn; coupSMPtr = coupSMPtrIn;
    partonSystemsPtr = partonSystemsPtrIn; userHooksPtr = userHooksPtrIn;
    mergingHooksPtr = mergingHooksPtrIn; partonVertexPtr = partonVertexPtrIn;
  }

  // Initialize alphaStrong and related pTmin parameters.
  virtual void init( BeamParticle* beamAPtrIn = 0,
    BeamParticle* beamBPtrIn = 0);

  // New beams possible for handling of hard diffraction. (Not virtual.)
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int beamOffsetIn = 0) {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamOffset = beamOffsetIn;}

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0.,
    double Q2Ren = 0.);

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() {return pTmaxFudge;}

  // Top-level routine to do a full time-like shower in resonance decay.
  virtual int shower( int iBeg, int iEnd, Event& event, double pTmax,
    int nBranchMax = 0);

  // Top-level routine for QED radiation in hadronic decay to two leptons.
  virtual int showerQED( int i1, int i2, Event& event, double pTmax);

  // Provide the pT scale of the last branching in the above shower.
  double pTLastInShower() {return pTLastBranch;}

  // Global recoil: reset counters and store locations of outgoing partons.
  virtual void prepareGlobal( Event& event);

  // Prepare system for evolution after each new interaction; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after a multiparton interactions rescattering.
  virtual void rescatterUpdate( int iSys, Event& event);

  // Update dipole list after each ISR emission.
  virtual void update( int iSys, Event& event, bool hasWeakRad = false);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll,
    bool isFirstTrial = false, bool doTrialIn = false);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event, bool isInterleaved = false);

  // Initialize data members for calculation of uncertainty bands.
  bool initUncertainties();

  // Calculate uncertainty-band weights for accepted/rejected trial branching.
  void calcUncertainties(bool accept, double pAccept, double enhance,
    double vp, TimeDipoleEnd* dip, Particle* radPtr, Particle* emtPtr,
    Particle* recPtr);

  // Tell which system was the last processed one.
  virtual int system() const {return iSysSel;};

  // Tell whether FSR has done a weak emission.
  bool getHasWeaklyRadiated() {return hasWeaklyRadiated;}

  // Print dipole list; for debug mainly.
  virtual void list() const;

  // Functions to allow usage of shower kinematics, evolution variables,
  // and splitting probabilities outside of shower.
  // Virtual so that shower plugins can overwrite these functions.
  // This makes it possible for another piece of the code to request
  // these - which is very convenient for merging.
  // Function variable names are not included to avoid compiler warnings.
  // Please see the documentation under "Implement New Showers" for details.

  // Return clustering kinematics - as needed for merging.
  virtual Event clustered( const Event& , int , int , int , string )
    { return Event();}

  // Return the evolution variable(s).
  // Important note: this map must contain the following entries
  // - a key "t" for the value of the shower evolution variable;
  // - a key "tRS" for the value of the shower evolution variable
  //   from which the shower would be restarted after a branching;
  // - a key "scaleAS" for the argument of alpha_s used for the branching;
  // - a key "scalePDF" for the argument of the PDFs used for the branching.
  // Usage: getStateVariables( event, iRad, iEmt, iRec,  name)
  virtual map<string, double> getStateVariables (const Event& , int , int ,
    int , string ) { return map<string,double>();}

  // Check if attempted clustering is handled by timelike shower
  // Usage: isTimelike( event, iRad, iEmt, iRec, name)
  virtual bool isTimelike(const Event& , int , int , int , string )
    { return false; }

  // Return a string identifier of a splitting.
  // Usage: getSplittingName( event, iRad, iEmt, iRec)
  virtual vector<string> getSplittingName( const Event& , int, int , int)
    { return vector<string>();}

  // Return the splitting probability.
  // Usage: getSplittingProb( event, iRad, iEmt, iRec)
  virtual double getSplittingProb( const Event& , int , int , int , string )
    { return 0.;}

  virtual bool allowedSplitting( const Event& , int , int)
    { return true; }
  virtual vector<int> getRecoilers( const Event&, int, int, string)
    { return vector<int>(); }

  // Pointer to MergingHooks object for NLO merging.
  MergingHooks* mergingHooksPtr;

protected:

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the settings database.
  Settings*      settingsPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // Pointer to the random number generator.
  Rndm*          rndmPtr;

  // Pointer to Standard Model couplings.
  CoupSM*        coupSMPtr;

  // Pointers to the two incoming beams. Offset their location in event.
  BeamParticle*  beamAPtr;
  BeamParticle*  beamBPtr;
  int            beamOffset;

  // Pointer to information on subcollision parton locations.
  PartonSystems* partonSystemsPtr;

  // Pointer to userHooks object for user interaction with program.
  UserHooks*     userHooksPtr;

  // Pointer to assign space-time vertices during parton evolution.
  PartonVertex*  partonVertexPtr;

  // Weak matrix elements used for corrections both of ISR and FSR.
  WeakShowerMEs  weakShowerMEs;

  // Store properties to be returned by methods.
  int    iSysSel;
  double pTmaxFudge, pTLastBranch;

private:

  // Constants: could only be changed in the code itself.
  static const double MCMIN, MBMIN, SIMPLIFYROOT, XMARGIN, XMARGINCOMB,
         TINYPDF, LARGEM2, THRESHM2, LAMBDA3MARGIN, WEAKPSWEIGHT, WG2QEXTRA,
         REJECTFACTOR, PROBLIMIT;
  // Rescatter: try to fix up recoil between systems
  static const bool   FIXRESCATTER, VETONEGENERGY;
  static const double MAXVIRTUALITYFRACTION, MAXNEGENERGYFRACTION;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, doQEDshowerByOther,
         doQEDshowerByGamma, doWeakShower, doMEcorrections, doMEextended,
         doMEafterFirst, doPhiPolAsym, doPhiPolAsymHard, doInterleave,
         allowBeamRecoil, dampenBeamRecoil, recoilToColoured, useFixedFacScale,
         allowRescatter, canVetoEmission, doHVshower, brokenHVsym,
         globalRecoil, useLocalRecoilNow, doSecondHard, hasUserHooks,
         singleWeakEmission, alphaSuseCMW, vetoWeakJets, allowMPIdipole,
         weakExternal, recoilDeadCone, doUncertainties, uVarMuSoftCorr,
         uVarMPIshowers, doDipoleRecoil, doPartonVertex, noResVariations,
         noProcVariations;
  int    pTmaxMatch, pTdampMatch, alphaSorder, alphaSnfmax, nGluonToQuark,
         weightGluonToQuark, alphaEMorder, nGammaToQuark, nGammaToLepton,
         nCHV, idHV, alphaHVorder, nMaxGlobalRecoil, weakMode;
  double pTdampFudge, mc, mb, m2c, m2b, renormMultFac, factorMultFac,
         fixedFacScale2, alphaSvalue, alphaS2pi, Lambda3flav, Lambda4flav,
         Lambda5flav, Lambda3flav2, Lambda4flav2, Lambda5flav2,
         scaleGluonToQuark, extraGluonToQuark, pTcolCutMin, pTcolCut,
         pT2colCut, pTchgQCut, pT2chgQCut, pTchgLCut, pT2chgLCut,
         pTweakCut, pT2weakCut, mMaxGamma, m2MaxGamma, octetOniumFraction,
         octetOniumColFac, mZ, gammaZ, thetaWRat, mW, gammaW, CFHV, nFlavHV,
         alphaHVfix, LambdaHV, pThvCut, pT2hvCut, mHV, pTmaxFudgeMPI,
         weakEnhancement, vetoWeakDeltaR2, dASmax, cNSpTmin, uVarpTmin2,
         overFactor;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM     alphaEM;

  // Some current values.
  bool   dopTlimit1, dopTlimit2, dopTdamp, hasWeaklyRadiated;
  double pT2damp, kRad, kEmt, pdfScale2;

  // Bookkeeping of enhanced  actual or trial emissions (see EPJC (2013) 73).
  bool doTrialNow, canEnhanceEmission, canEnhanceTrial, canEnhanceET,
       doUncertaintiesNow;
  string splittingNameNow, splittingNameSel;
  map< double, pair<string,double> > enhanceFactors;
  void storeEnhanceFactor(double pT2, string name, double enhanceFactorIn)
    { enhanceFactors.insert(make_pair(pT2,make_pair(name,enhanceFactorIn)));}

  // All dipole ends and a pointer to the selected hardest dipole end.
  vector<TimeDipoleEnd> dipEnd;
  TimeDipoleEnd* dipSel;
  int iDipSel;

  // Setup a dipole end, either QCD, QED/photon, weak or Hidden Valley one.
  void setupQCDdip( int iSys, int i, int colTag,  int colSign, Event& event,
    bool isOctetOnium = false, bool limitPTmaxIn = true);
  void setupQEDdip( int iSys, int i, int chgType, int gamType, Event& event,
    bool limitPTmaxIn = true);
  void setupWeakdip( int iSys, int i,int weakType, Event& event,
    bool limitPTmaxIn = true);
  // Special setup for weak dipoles if already specified in info ptr.
  void setupWeakdipExternal(Event& event, bool limitPTmaxIn = true);
  void setupHVdip( int iSys, int i, Event& event, bool limitPTmaxIn = true);

  // Evolve a QCD dipole end.
  void pT2nextQCD( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a QED dipole end, either charged or photon.
  void pT2nextQED( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a weak dipole end.
  void pT2nextWeak( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a Hidden Valley dipole end.
  void pT2nextHV( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& );

  // Find kind of QCD ME correction.
  void findMEtype( Event& event, TimeDipoleEnd& dip);

  // Find type of particle; used by findMEtype.
  int findMEparticle( int id, bool isHiddenColour = false);

  // Find mixture of V and A in gamma/Z: energy- and flavour-dependent.
  double gammaZmix( Event& event, int iRes, int iDau1, int iDau2);

  // Set up to calculate QCD ME correction with calcMEcorr.
  double findMEcorr(TimeDipoleEnd* dip, Particle& rad, Particle& partner,
   Particle& emt, bool cutEdge = true);

  // Set up to calculate weak ME corrections for t-channel processes.
  double findMEcorrWeak(TimeDipoleEnd* dip, Vec4 rad, Vec4 rec,
   Vec4 emt, Vec4 p3, Vec4 p4, Vec4 radBef, Vec4 recBef);

  // Calculate value of QCD ME correction.
  double calcMEcorr( int kind, int combiIn, double mixIn, double x1,
    double x2, double r1, double r2, double r3 = 0., bool cutEdge = true);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, TimeDipoleEnd* dip);

  // Rescatter: propagate dipole recoil to internal lines connecting systems.
  bool rescatterPropagateRecoil( Event& event, Vec4& pNew);

  // Properties stored for (some) global recoil schemes.
  // Vectors of event indices defining the hard process.
  vector<int> hardPartons;
  // Number of partons in current hard event, number of partons in Born-type
  // hard event (to distinguish between S and H), maximally allowed number of
  // global recoil branchings.
  int nHard, nFinalBorn, nMaxGlobalBranch;
  // Number of proposed splittings in hard scattering systems.
  map<int,int> nProposed;
  // Number of splittings with global recoil (currently only 1).
  int nGlobal, globalRecoilMode;
  // Switch to constrain recoiling system.
  bool limitMUQ;

  // 2 -> 2 information needed for the external weak setup.
  vector<Vec4> weakMomenta;
  vector<int> weak2to2lines;
  int weakHardSize;

  // Store uncertainty variations relevant to TimeShower.
  int nUncertaintyVariations, nVarQCD, uVarNflavQ;
  map<int,double> varG2GGmuRfac, varQ2QGmuRfac, varG2QQmuRfac, varX2XGmuRfac;
  map<int,double> varG2GGcNS, varQ2QGcNS, varG2QQcNS, varX2XGcNS;
  map<int,double> varPDFplus, varPDFminus;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_TimeShower_H
