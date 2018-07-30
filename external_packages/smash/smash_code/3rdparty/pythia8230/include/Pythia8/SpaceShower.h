// SpaceShower.h is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the spacelike initial-state showers.
// SpaceDipoleEnd: radiating dipole end in ISR.
// SpaceShower: handles the showering description.

#ifndef Pythia8_SpaceShower_H
#define Pythia8_SpaceShower_H

#include "Pythia8/Basics.h"
#include "Pythia8/BeamParticle.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/PartonVertex.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/UserHooks.h"
#include "Pythia8/MergingHooks.h"
#include "Pythia8/WeakShowerMEs.h"


namespace Pythia8 {

//==========================================================================

// Data on radiating dipole ends, only used inside SpaceShower.

class SpaceDipoleEnd {

public:

  // Constructor.
  SpaceDipoleEnd( int systemIn = 0, int sideIn = 0, int iRadiatorIn = 0,
    int iRecoilerIn = 0, double pTmaxIn = 0., int colTypeIn = 0,
    int chgTypeIn = 0, int weakTypeIn = 0,  int MEtypeIn = 0,
    bool normalRecoilIn = true, int weakPolIn = 0,
    int iColPartnerIn = 0, int idColPartnerIn = 0) :
    system(systemIn), side(sideIn), iRadiator(iRadiatorIn),
    iRecoiler(iRecoilerIn), pTmax(pTmaxIn), colType(colTypeIn),
    chgType(chgTypeIn), weakType(weakTypeIn), MEtype(MEtypeIn),
    normalRecoil(normalRecoilIn), weakPol(weakPolIn),
    iColPartner(iColPartnerIn), idColPartner(idColPartnerIn),
    nBranch(0), pT2Old(0.), zOld(0.5) { }

  // Store values for trial emission.
  void store( int idDaughterIn, int idMotherIn, int idSisterIn,
    double x1In, double x2In, double m2DipIn, double pT2In, double zIn,
    double xMoIn, double Q2In, double mSisterIn, double m2SisterIn,
    double pT2corrIn, int iColPartnerIn, double m2IFIn, double mColPartnerIn)
    {idDaughter = idDaughterIn; idMother = idMotherIn;
    idSister = idSisterIn; x1 = x1In; x2 = x2In; m2Dip = m2DipIn;
    pT2 = pT2In; z = zIn; xMo = xMoIn; Q2 = Q2In; mSister = mSisterIn;
    m2Sister = m2SisterIn; pT2corr = pT2corrIn; iColPartner = iColPartnerIn;
    m2IF = m2IFIn; mColPartner = mColPartnerIn;}

  // Basic properties related to evolution and matrix element corrections.
  int    system, side, iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, weakType, MEtype;
  bool   normalRecoil;
  int    weakPol, iColPartner, idColPartner;

  // Properties specific to current trial emission.
  int    nBranch, idDaughter, idMother, idSister, iFinPol;
  double x1, x2, m2Dip, pT2, z, xMo, Q2, mSister, m2Sister, pT2corr,
         pT2Old, zOld, asymPol, m2IF, mColPartner;

  // Properties needed for the evaluation of parameter variations
  double pAccept;

} ;

//==========================================================================

// The SpaceShower class does spacelike showers.

class SpaceShower {

public:

  // Constructor.
  SpaceShower() {beamOffset = 0;}

  // Destructor.
  virtual ~SpaceShower() {}

  // Initialize various pointers.
  // (Separated from rest of init since not virtual.)
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn, CoupSM* coupSMPtrIn,
    PartonSystems* partonSystemsPtrIn, UserHooks* userHooksPtrIn,
    MergingHooks* mergingHooksPtrIn, PartonVertex* partonVertexPtrIn) {
    infoPtr = infoPtrIn; settingsPtr = settingsPtrIn;
    particleDataPtr = particleDataPtrIn; rndmPtr = rndmPtrIn;
    coupSMPtr = coupSMPtrIn; partonSystemsPtr = partonSystemsPtrIn;
    userHooksPtr = userHooksPtrIn; mergingHooksPtr = mergingHooksPtrIn;
    partonVertexPtr = partonVertexPtrIn; }

  // Initialize generation. Possibility to force re-initialization by hand.
  virtual void init(BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);

  // New beams possible for handling of hard diffraction. (Not virtual.)
  void reassignBeamPtrs( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn,
    int beamOffsetIn = 0) {beamAPtr = beamAPtrIn; beamBPtr = beamBPtrIn;
    beamOffset = beamOffsetIn;}

  // Find whether to limit maximum scale of emissions, and whether to dampen.
  virtual bool limitPTmax( Event& event, double Q2Fac = 0.,
    double Q2Ren = 0.);

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() const {return pTmaxFudge;}

  // Prepare system for evolution; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmaxIn = true);

  // Update dipole list after each FSR emission.
  // Usage: update( iSys, event).
  virtual void update( int , Event&, bool hasWeakRad = false);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll,
    int nRadIn = -1, bool doTrialIn = false);

  // ME corrections and kinematics that may give failure.
  virtual bool branch( Event& event);

  // Initialize data members for calculation of uncertainty bands.
  bool initUncertainties();

  // Calculate uncertainty-band weights for accepted/rejected trial branching.
  void calcUncertainties(bool accept, double pAcceptIn, double pT20in,
    double enhance, double vp, SpaceDipoleEnd* dip, Particle* motherPtr,
    Particle* sisterPtr);

  // Tell if latest scattering was a gamma->qqbar.
  bool wasGamma2qqbar() { return gamma2qqbar; }

  // Tell which system was the last processed one.
  virtual int system() const {return iSysSel;}

  // Flag for failure in branch(...) that will force a retry of parton level.
  bool doRestart() const {return rescatterFail;}

  // Tell whether ISR has done a weak emission.
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

  // Return clustering kinematics - as needed form merging.
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

  // Check if attempted clustering is handled by spacelike shower.
  // Usage: isSpacelike( event, iRad, iEmt, iRec, name)
  virtual bool isSpacelike(const Event&, int, int, int, string)
    { return false; }

  // Return a string identifier of a splitting.
  // Usage: getSplittingName( event, iRad, iEmt, iRec)
  virtual vector<string> getSplittingName( const Event& , int , int , int )
    { return vector<string>();}

  // Return the splitting probability.
  // Usage: getSplittingProb( event, iRad, iEmt, iRec)
  virtual double getSplittingProb( const Event& , int , int , int , string )
    { return 0.;}

  virtual bool allowedSplitting( const Event& , int , int)
    { return true;}
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
  bool   rescatterFail;
  int    iSysSel;
  double pTmaxFudge;

private:

  // Constants: could only be changed in the code itself.
  static const int    MAXLOOPTINYPDF;
  static const double MCMIN, MBMIN, CTHRESHOLD, BTHRESHOLD, EVALPDFSTEP,
         TINYPDF, TINYKERNELPDF, TINYPT2, HEAVYPT2EVOL, HEAVYXEVOL,
         EXTRASPACEQ, LAMBDA3MARGIN, PT2MINWARN, LEPTONXMIN, LEPTONXMAX,
         LEPTONPT2MIN, LEPTONFUDGE, WEAKPSWEIGHT, HEADROOMQ2Q, HEADROOMQ2G,
         HEADROOMG2G, HEADROOMG2Q, HEADROOMHQG, REJECTFACTOR, PROBLIMIT;

  // Initialization data, normally only set once.
  bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, useSamePTasMPI,
         doWeakShower, doMEcorrections, doMEafterFirst, doPhiPolAsym,
         doPhiPolAsymHard, doPhiIntAsym, doRapidityOrder, useFixedFacScale,
         doSecondHard, canVetoEmission, hasUserHooks, alphaSuseCMW,
         singleWeakEmission, vetoWeakJets, weakExternal, doRapidityOrderMPI,
         doUncertainties, uVarMuSoftCorr, uVarMPIshowers, doMPI, gamma2qqbar,
         doDipoleRecoil, doPartonVertex;
  int    pTmaxMatch, pTdampMatch, alphaSorder, alphaSnfmax, alphaEMorder,
         nQuarkIn, enhanceScreening, weakMode;
  double pTdampFudge, mc, mb, m2c, m2b, renormMultFac, factorMultFac,
         fixedFacScale2, alphaSvalue, alphaS2pi, Lambda3flav, Lambda4flav,
         Lambda5flav, Lambda3flav2, Lambda4flav2, Lambda5flav2, pT0Ref,
         ecmRef, ecmPow, pTmin, sCM, eCM, pT0, pTminChgQ, pTminChgL, pT20,
         pT2min, pT2minChgQ, pT2minChgL, pTweakCut, pT2weakCut, pTmaxFudgeMPI,
         strengthIntAsym, weakEnhancement, mZ, gammaZ, thetaWRat, mW, gammaW,
         weakMaxWt, vetoWeakDeltaR2, dASmax, cNSpTmin, uVarpTmin2, overFactor;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM alphaEM;

  // Some current values.
  bool   sideA, dopTlimit1, dopTlimit2, dopTdamp, hasWeaklyRadiated, tChannel,
         doUncertaintiesNow;
  int    iNow, iRec, idDaughter, nRad, idResFirst, idResSecond;
  double xDaughter, x1Now, x2Now, m2ColPair, mColPartner, m2ColPartner,
         m2Dip, m2Rec, pT2damp, pTbegRef, pdfScale2;

  // Bookkeeping of enhanced  actual or trial emissions (see EPJC (2013) 73).
  bool doTrialNow, canEnhanceEmission, canEnhanceTrial, canEnhanceET;
  string splittingNameNow, splittingNameSel;
  map< double, pair<string,double> > enhanceFactors;
  void storeEnhanceFactor(double pT2, string name, double enhanceFactorIn)
    { enhanceFactors.insert(make_pair(pT2,make_pair(name,enhanceFactorIn)));}

  // List of emissions in different sides in different systems:
  vector<int> nRadA,nRadB;

  // All dipole ends
  vector<SpaceDipoleEnd> dipEnd;

  // List of 2 -> 2 momenta for external weak setup.
  vector<Vec4> weakMomenta;

  // Pointers to the current and hardest (so far) dipole ends.
  int iDipNow, iSysNow;
  SpaceDipoleEnd* dipEndNow;
  int iDipSel;
  SpaceDipoleEnd* dipEndSel;

  // Evolve a QCD dipole end.
  void pT2nextQCD( double pT2begDip, double pT2endDip);

  // Evolve a QCD and QED dipole end near heavy quark threshold region.
  void pT2nearThreshold( BeamParticle& beam, double m2Massive,
    double m2Threshold, double xMaxAbs, double zMinAbs,
    double zMaxMassive, int iColPartner);

  // Evolve a QED dipole end.
  void pT2nextQED( double pT2begDip, double pT2endDip);

  // Evolve a Weak dipole end.
  void pT2nextWeak( double pT2begDip, double pT2endDip);

  // Find class of ME correction.
  int findMEtype( int iSys, Event& event, bool weakRadiation = false);

  // Provide maximum of expected ME weight; for preweighting of evolution.
  double calcMEmax( int MEtype, int idMother, int idDaughterIn);

  // Provide actual ME weight for current branching.
  double calcMEcorr(int MEtype, int idMother, int idDaughterIn, double M2,
    double z, double Q2,double m2Sister);

  // Provide actual ME weight for t-channel weak emissions.
  double calcMEcorrWeak(int MEtype, double m2, double z,
    double pT2, Vec4 pMother, Vec4 pB, Vec4 pDaughter,
    Vec4 pB0, Vec4 p1, Vec4 p2, Vec4 pSister);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, SpaceDipoleEnd* dip);

  // Store uncertainty variations relevant to TimeShower.
  int nUncertaintyVariations, nVarQCD, uVarNflavQ;
  map<int,double> varG2GGmuRfac, varQ2QGmuRfac, varQ2GQmuRfac, varG2QQmuRfac,
    varX2XGmuRfac;
  map<int,double> varG2GGcNS, varQ2QGcNS, varQ2GQcNS, varG2QQcNS, varX2XGcNS;
  map<int,double> varPDFplus, varPDFminus;

  // Find a possible colour partner in the case of dipole recoil.
  int findColPartner(Event& event, int iSideA, int iSideB, int iSystem);

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_SpaceShower_H
