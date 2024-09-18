/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2019
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

#include "HybridHadronization.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

#include "FluidCellInfo.h"
#include "FluidDynamics.h"
#include "FluidEvolutionHistory.h"
#include "JetScapeConstants.h"
#include "JetScapeLogger.h"
#include "JetScapeSignalManager.h"
#include "JetScapeXML.h"
#include "SurfaceCellInfo.h"
#include "ThermPtnSampler.h"
#include "tinyxml2.h"
//#include <cmath>

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<HybridHadronization> HybridHadronization::reg(
    "HybridHadronization");

// Initialize static helper here
Pythia8::Pythia HybridHadronization::pythia("IntentionallyEmpty", false);

// RNG - Mersenne Twist - 64 bit
// std::mt19937_64 eng(std::random_device{}());
// std::mt19937_64 eng(1);
// returns a random number between 0 and 1, based on above engine
double HybridHadronization::ran() {
  std::uniform_real_distribution<double> uniran(0.0, 1.0);
  return uniran(eng);
}

HybridHadronization::HybridHadronization() {
  SetId("HybridHadronization");
  VERBOSE(8);
}

HybridHadronization::~HybridHadronization() { VERBOSE(8); }

// meson width function
double HybridHadronization::SigM2_calc(double R2chg, double qm1, double qm2,
                                       double qq1, double qq2) {
  return R2chg * (2. / 3.) * (qm1 + qm2) * (qm1 + qm2) /
         (std::abs(qq1) * qm2 * qm2 + std::abs(qq2) * qm1 * qm1) *
         (std::abs(qq1) + std::abs(qq2));
}

// baryon width function
double HybridHadronization::SigBR2_calc(double R2chg, double qm1, double qm2,
                                        double qm3, double qq1, double qq2,
                                        double qq3) {
  return R2chg * (2. / 3.) * (qm1 + qm2 + qm3) /
         (std::abs(qq1) * qm2 * (qm2 + qm3) / (qm1 + qm2) +
          std::abs(qq2) * qm1 * (qm1 + qm3) / (qm1 + qm2) +
          std::abs(qq3) * (qm1 * qm2) / qm3) *
         (std::abs(qq1) + std::abs(qq2) + std::abs(qq3));
}

double HybridHadronization::SigBL2_calc(double SigBR2, double qm1, double qm2,
                                        double qm3) {
  return SigBR2 * ((qm1 * qm2) / (qm1 + qm2)) /
         (qm3 * (qm1 + qm2) / (qm1 + qm2 + qm3));
}

void HybridHadronization::Init() {
  tinyxml2::XMLElement* hadronization = GetXMLElement({"JetHadronization"});

  if (!hadronization) {
    JSWARN << "Couldn't find tag Jet Hadronization";
    throw std::runtime_error("Couldn't find tag Jet Hadronization");
  }
  if (hadronization) {
    string s = hadronization->FirstChildElement("name")->GetText();
    // std::string s = GetXMLElementText({"JetHadronization", "name"});
    JSDEBUG << s << " to be initialized ...";
    JSINFO << "Initialize Hybrid Hadronization ...";

    JSDEBUG << "Initialize HybridHadronization";
    VERBOSE(8);

    maxM_level = 1;  // maximum energy level considered for meson recombination;
                     // maximum: 4
    maxB_level =
        1;  // maximum energy level considered for baryon recombination; no
            // maximum but no physical baryons beyond 0
    goldstonereco = false;  // don't allow recombination of Goldstone bosons
    gmax = 1.25;  // maximum allowed mass of the gluon (for q-qbar split), in
                  // GeV
    xmq = 0.33;   // light quark mass, in GeV
    xms = 0.5;    // strange quark mass, in GeV
    xmc = 1.5;    // charm quark mass in GeV
    xmb = 4.8;    // bottom quark mass in GeV
    hbarc = 0.197327;  // GeV*fm - maybe just set this as a constant in common?
    dist2cut = 25.;    // maximum distance [fm] squared for recombination
                       // (involving thermal partons) - in lab frame
    sh_recofactor =
        1.;  // suppression/enhancement factor for shower-shower recombination
    th_recofactor =
        1.;  // suppression/enhancement factor for shower-thermal recombination
    attempts_max = 15;  // maximum number of failed attempts to hadronize a
                        // single event before we give up.
    p_fake = 0.;        // momentum used for fake parton, if needed
    rand_seed = 0;  // seed for RNGs used - 0 means use a randomly determined
                    // random seed (from system time or std::random_device{}())
    had_prop =
        0.;  // propagation of hadrons after formation by this time in lab frame
    part_prop = 0.;  // minimum propagation time of partons after last split
    reco_hadrons_pythia = 0;  // flag to put recombination hadrons into pythia
                              // for decays (position would be lost)

    // xml read in to alter settings...
    double xml_doublein = -1.;
    int xml_intin = -1;
    unsigned int xml_uintin = std::numeric_limits<unsigned int>::max();

    xml_doublein =
        GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
    if (xml_doublein >= 0.) {
      p_fake = xml_doublein / 6.;
    }
    xml_doublein =
        -1.;  // for colliders set fake momentum to valence quark energy

    xml_doublein =
        GetXMLElementDouble({"JetHadronization", "thermreco_distmax"});
    if (xml_doublein >= 0.) {
      dist2cut = xml_doublein * xml_doublein;
    }
    xml_doublein = -1.;

    xml_doublein =
        GetXMLElementDouble({"JetHadronization", "shower_recofactor"});
    if (xml_doublein >= 0.) {
      sh_recofactor = xml_doublein;
    }
    xml_doublein = -1.;

    xml_doublein =
        GetXMLElementDouble({"JetHadronization", "thermal_recofactor"});
    if (xml_doublein >= 0.) {
      th_recofactor = xml_doublein;
    }
    xml_doublein = -1.;

    xml_intin = GetXMLElementInt({"JetHadronization", "reco_Mlevelmax"});
    if (xml_intin >= 0) {
      maxM_level = xml_intin;
    }
    xml_intin = -1;

    xml_intin = GetXMLElementInt({"JetHadronization", "reco_Blevelmax"});
    if (xml_intin >= -1) {
      maxB_level = xml_intin;
    }
    xml_intin = -1;

    xml_intin = GetXMLElementInt({"JetHadronization", "reco_goldstone"});
    if (xml_intin >= 0) {
      goldstonereco = xml_intin;
    }
    xml_intin = -1;

    torder_reco = false;
    xml_intin = GetXMLElementInt({"JetHadronization", "recobias_t"});
    if (xml_intin == 1) {
      torder_reco = true;
    }
    xml_intin = -1;

    xml_doublein = GetXMLElementDouble({"JetHadronization", "hydro_Tc"});
    if (xml_doublein >= 0) {
      hydro_Tc = xml_doublein;
    }
    xml_doublein = -1;

    xml_doublein = GetXMLElementDouble({"JetHadronization", "had_postprop"});
    if (xml_doublein >= 0) {
      had_prop = xml_doublein;
    }
    xml_doublein = -1;

    xml_doublein = GetXMLElementDouble({"JetHadronization", "part_prop"});
    if (xml_doublein >= 0) {
      part_prop = xml_doublein;
    }
    xml_doublein = -1;

    xml_intin =
        GetXMLElementInt({"JetHadronization", "reco_hadrons_in_pythia"});
    if (xml_intin == 0 || xml_intin == 1) {
      reco_hadrons_pythia = xml_intin;
    }
    xml_intin = -1;

    xml_intin =
        GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"});
    if (xml_intin == 1) {
      afterburner_frag_hadrons = true;
    }
    xml_intin = -1;

    if (maxM_level > 4) {
      maxM_level = 4;
      JSWARN << "Requested maximum energy level for mesons too large. Set it "
                "to 4.";
    }
    // No, really, the maximum is four

    // random seed
    // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
    tinyxml2::XMLElement* RandomXmlDescription = GetXMLElement({"Random"});
    if (RandomXmlDescription) {
      tinyxml2::XMLElement* xmle;
      xmle = RandomXmlDescription->FirstChildElement("seed");
      xmle->QueryUnsignedText(&xml_uintin);
      if (xml_uintin < std::numeric_limits<unsigned int>::max()) {
        rand_seed = xml_uintin;
      }
      xml_uintin = std::numeric_limits<unsigned int>::max();
    } else {
      JSWARN << "No <Random> element found in xml, seeding to 0";
    }
    VERBOSE(7) << "Seeding PYTHIA(hadronization) to " << rand_seed;

    // not sure if we can seed with a negative integer...
    if (rand_seed != 0) {
      eng.seed(rand_seed);
    }
    // else{eng.seed(std::random_device{}());}
    else {  // seeding the mt19937_64 object 'eng' PROPERLY!
      std::random_device rd;
      std::array<int, std::mt19937_64::state_size> seedarray;
      std::generate_n(seedarray.data(), seedarray.size(), std::ref(rd));
      std::seed_seq seeds(std::begin(seedarray), std::end(seedarray));
      eng.seed(seeds);
    }

    // Since PYTHIA has no spacetime information, we shouldn't use recombination
    // as it is necessary to calculate recombination probabilities later, this
    // will instead update to set a flag to attempt to artificially generate
    // this information for now, we just print a warning - but still try to run
    // recombination.  It will just be unphysically enhanced, esp. for certain
    // configurations tinyxml2::XMLElement
    // *XmlPythiaGun=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard"
    // )->FirstChildElement("PythiaGun");
    /*	tinyxml2::XMLElement *XmlPythiaGun = GetXMLElement({"Hard",
       "PythiaGun"}); bool PYgun_FSRon = false; if ( XmlPythiaGun ){
                  tinyxml2::XMLElement *xmle; xmle =
       XmlPythiaGun->FirstChildElement( "FSR_on" );
                  xmle->QueryIntText(&xml_intin);
                  if(xml_intin == 1){PYgun_FSRon = true;} xml_intin = -1;
          }
          if(PYgun_FSRon && (sh_recofactor > 0.0000000001)){JSWARN <<
       "Recombination with a PYTHIA FSR shower is not fully implemented.";}
    */
    xml_intin = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
    if (xml_intin && (sh_recofactor > 0.0000000001)) {
      JSWARN
          << "Recombination with a PYTHIA FSR shower is not fully implemented.";
    }
    xml_intin = -1;

    // quark masses/charges used to get widths from charged radii (use Pythia
    // masses)
    double Qm_ud = xmq;
    double Qm_s = xms;
    double Qm_c = xmc;
    double Qm_b = xmb;
    double chg_u = 2. / 3.;
    double chg_d = -1. / 3.;

    // rms charge radii
    // mesons
    double R2chg_Pi = 0.42;
    double R2chg_Phi = 0.21;
    double R2chg_K = 0.34;
    double R2chg_Jpi = 0.04;
    double R2chg_Ds = 0.09;
    double R2chg_D = 0.165;
    double R2chg_Ups =
        0.032;  //????? setting to a linear extrapolation from Jpsi -> Bc -> Ups
    double R2chg_Bc = 0.036;
    double R2chg_B = 0.273;
    // baryons
    double R2chg_Nuc = 0.69;
    double R2chg_Omg = 0.355;
    double R2chg_Xi = 0.52;
    double R2chg_Sig = 0.61;
    double R2chg_Occc = 0.179;  //?????
    double R2chg_Occ = 0.043;   //?
    double R2chg_Xicc = 0.049;  //?
    double R2chg_Oc = 0.1;      //?????
    double R2chg_Xic = 0.24;    //?
    double R2chg_Sigc = 0.27;   //?
    double R2chg_Obbb = 0.001;  //?????
    double R2chg_Obbc = 0.001;  //?????
    double R2chg_Obb = 0.001;   //?????
    double R2chg_Xibb = 0.001;  //?????
    double R2chg_Obcc = 0.001;  //?????
    double R2chg_Obc = 0.001;   //?????
    double R2chg_Xibc = 0.001;  //?????
    double R2chg_Ob = 0.6;      //?????
    double R2chg_Xib = 0.63;    //?????
    double R2chg_Sigb = 0.66;   //?????

    // meson width calculations (r2) - recalc if r2chg is changed on command
    // line...
    SigPi2 = SigM2_calc(R2chg_Pi, Qm_ud, Qm_ud, chg_d, chg_u);
    SigPhi2 = SigM2_calc(R2chg_Phi, Qm_s, Qm_s, chg_d, -chg_d);  // normalizing
    SigK2 = SigM2_calc(R2chg_K, Qm_s, Qm_ud, chg_d, chg_u);
    SigJpi2 = SigM2_calc(R2chg_Jpi, Qm_c, Qm_c, chg_u, -chg_u);  // normalizing
    SigDs2 = SigM2_calc(R2chg_Ds, Qm_c, Qm_s, chg_u, chg_d);
    SigD2 = SigM2_calc(R2chg_D, Qm_c, Qm_ud, chg_u, chg_d);
    SigUps2 = SigM2_calc(R2chg_Ups, Qm_b, Qm_b, chg_d, -chg_d);  // normalizing
    SigBc2 = SigM2_calc(R2chg_Bc, Qm_b, Qm_c, chg_d, chg_u);
    SigB2 =
        SigM2_calc(R2chg_B, Qm_b, Qm_ud, chg_d, chg_u);  // (treating B_s as B)

    // baryon width calculations (r2) - recalc if r2chg is changed on command
    // line... light/strange baryons
    SigNucR2 = SigBR2_calc(R2chg_Nuc, Qm_ud, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
    SigNucL2 = SigBL2_calc(SigNucR2, Qm_ud, Qm_ud, Qm_ud);
    SigOmgR2 = SigBR2_calc(R2chg_Omg, Qm_s, Qm_s, Qm_s, chg_d, chg_d, chg_d);
    SigOmgL2 = SigBL2_calc(SigOmgR2, Qm_s, Qm_s, Qm_s);
    SigXiR2 = SigBR2_calc(R2chg_Xi, Qm_s, Qm_s, Qm_ud, chg_d, chg_d, chg_d);
    SigXiL2 = SigBL2_calc(SigXiR2, Qm_s, Qm_s, Qm_ud);
    SigSigR2 = SigBR2_calc(R2chg_Sig, Qm_s, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
    SigSigL2 = SigBL2_calc(SigSigR2, Qm_s, Qm_ud, Qm_ud);

    // charm baryons
    SigOcccR2 = SigBR2_calc(R2chg_Occc, Qm_c, Qm_c, Qm_c, chg_d, chg_d,
                            chg_d);  // ! maybe need to normalize? (just setting
                                     // all to -1/3 for now)
    SigOcccL2 = SigBL2_calc(SigOcccR2, Qm_c, Qm_c, Qm_c);
    SigOccR2 = SigBR2_calc(R2chg_Occ, Qm_c, Qm_c, Qm_s, chg_u, chg_u, chg_d);
    SigOccL2 = SigBL2_calc(SigOccR2, Qm_c, Qm_c, Qm_s);
    SigXiccR2 = SigBR2_calc(R2chg_Xicc, Qm_c, Qm_c, Qm_ud, chg_u, chg_u, chg_d);
    SigXiccL2 = SigBL2_calc(SigXiccR2, Qm_c, Qm_c, Qm_ud);
    SigOcR2 = SigBR2_calc(R2chg_Oc, Qm_c, Qm_s, Qm_s, chg_d, chg_d,
                          chg_d);  // ! setting all quark charges to -1/3
    SigOcL2 = SigBL2_calc(SigOcR2, Qm_c, Qm_s, Qm_s);
    SigXicR2 = SigBR2_calc(R2chg_Xic, Qm_c, Qm_s, Qm_ud, chg_u, chg_d, chg_u);
    SigXicL2 = SigBL2_calc(SigXicR2, Qm_c, Qm_s, Qm_ud);
    SigSigcR2 =
        SigBR2_calc(R2chg_Sigc, Qm_c, Qm_ud, Qm_ud, chg_u, chg_d, chg_u);
    SigSigcL2 = SigBL2_calc(SigSigcR2, Qm_c, Qm_ud, Qm_ud);

    // bottom baryons
    SigObbbR2 = SigBR2_calc(R2chg_Obbb, Qm_b, Qm_b, Qm_b, chg_d, chg_d, chg_d);
    SigObbbL2 = SigBL2_calc(SigObbbR2, Qm_b, Qm_b, Qm_b);
    SigObbcR2 = SigBR2_calc(R2chg_Obbc, Qm_b, Qm_b, Qm_c, chg_d, chg_d,
                            chg_d);  // ! setting all quark charges to -1/3
    SigObbcL2 = SigBL2_calc(SigObbcR2, Qm_b, Qm_b, Qm_c);
    SigObbR2 = SigBR2_calc(R2chg_Obb, Qm_b, Qm_b, Qm_s, chg_d, chg_d, chg_d);
    SigObbL2 = SigBL2_calc(SigObbR2, Qm_b, Qm_b, Qm_s);
    SigXibbR2 = SigBR2_calc(R2chg_Xibb, Qm_b, Qm_b, Qm_ud, chg_d, chg_d, chg_d);
    SigXibbL2 = SigBL2_calc(SigXibbR2, Qm_b, Qm_b, Qm_ud);
    SigObccR2 = SigBR2_calc(R2chg_Obcc, Qm_b, Qm_c, Qm_c, chg_d, chg_u, chg_u);
    SigObccL2 = SigBL2_calc(SigObccR2, Qm_b, Qm_c, Qm_c);
    SigObcR2 = SigBR2_calc(R2chg_Obc, Qm_b, Qm_c, Qm_s, chg_d, chg_d,
                           chg_d);  // ! flipping c quark charge (all to -1/3)
    SigObcL2 = SigBL2_calc(SigObcR2, Qm_b, Qm_c, Qm_s);
    SigXibcR2 = SigBR2_calc(R2chg_Xibc, Qm_b, Qm_c, Qm_ud, chg_d, chg_u, chg_u);
    SigXibcL2 = SigBL2_calc(SigXibcR2, Qm_b, Qm_c, Qm_ud);
    SigObR2 = SigBR2_calc(R2chg_Ob, Qm_b, Qm_s, Qm_s, chg_d, chg_d, chg_d);
    SigObL2 = SigBL2_calc(SigObR2, Qm_b, Qm_s, Qm_s);
    SigXibR2 = SigBR2_calc(R2chg_Xib, Qm_b, Qm_s, Qm_ud, chg_d, chg_d, chg_d);
    SigXibL2 = SigBL2_calc(SigXibR2, Qm_b, Qm_s, Qm_ud);
    SigSigbR2 =
        SigBR2_calc(R2chg_Sigb, Qm_b, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
    SigSigbL2 = SigBL2_calc(SigSigbR2, Qm_b, Qm_ud, Qm_ud);

    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedParticleData = off");

    // Standard settings
    pythia.readString("ProcessLevel:all = off");
    // pythia.readString("PartonLevel:FSR=off"); //is this necessary?

    // General settings for hadron decays
    pythia_decays = GetXMLElementText({"JetHadronization", "pythia_decays"});
    double tau0Max = 10.0;
    double tau0Max_xml = GetXMLElementDouble({"JetHadronization", "tau0Max"});
    if (tau0Max_xml >= 0) {
      tau0Max = tau0Max_xml;
    } else {
      JSWARN << "tau0Max should be larger than 0. Set it to 10.";
    }
    if (pythia_decays == "on") {
      JSINFO << "Pythia decays are turned on for tau0Max < " << tau0Max;
      pythia.readString("HadronLevel:Decay = on");
      pythia.readString("ParticleDecays:limitTau0 = on");
      pythia.readString("ParticleDecays:tau0Max = " + std::to_string(tau0Max));
    } else {
      JSINFO << "Pythia decays are turned off";
      pythia.readString("HadronLevel:Decay = off");
    }

    // Settings for decays (old flag, will be depracted at some point)
    // This overwrites the previous settings if the user xml file contains the
    // flag
    std::string weak_decays =
        GetXMLElementText({"JetHadronization", "weak_decays"});
    if (weak_decays == "off") {
      JSINFO << "Hadron decays are turned off.";
      JSWARN << "This parameter will be depracted at some point. Use "
                "'pythia_decays' instead.\nOverwriting 'pythia_decays'.";
      pythia.readString("HadronLevel:Decay = off");
    } else if (weak_decays == "on") {
      JSINFO << "Hadron decays inside a range of 10 mm/c are turned on.";
      JSWARN << "This parameter will be depracted at some point. Use "
                "'pythia_decays' and 'tau0Max' for more control on "
                "decays.\nOverwriting 'pythia_decays' and fix 'tau0Max' to 10.";
      pythia.readString("HadronLevel:Decay = on");
      pythia.readString("ParticleDecays:limitTau0 = on");
      pythia.readString("ParticleDecays:tau0Max = 10.0");
    }

    // setting seed, or using random seed
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(rand_seed));

    // additional settings
    // turning off pythia checks for runtime decrease (can be turned back on if
    // necessary, but it shouldn't make much of a difference)
    pythia.readString(
        "Check:event = off");  // is probably a bad idea, but shouldn't really
                               // be necessary... will use a bit of runtime on
                               // event checks...
    pythia.readString(
        "Check:history = off");  // might be a good idea to set 'off' as it
                                 // saves runtime - provided we know that we've
                                 // set up mother/daughter relations
                                 // correctly...

    // making the pythia event checks a little less stringent (PYTHIA
    // documentation already states that LHC events will occasionally violate
    // default constraint, without concern) pythia.readString("Check:epTolWarn =
    // 1e-4");   // setting E/P conservation violation constraint somewhat
    // weaker, just for ease pythia.readString("Check:epTolErr  = 1e-2");   //
    // setting E/P conservation violation constraint somewhat weaker, just for
    // ease pythia.readString("Check:mTolWarn  = 1e-2");   // setting EP/M
    // conservation violation constraint somewhat weaker, just for ease
    // pythia.readString("Check:mTolErr   = 1e-1");   // setting EP/M
    // conservation violation constraint somewhat weaker, just for ease

    // allowing for partonic space-time information to be used by PYTHIA
    // pythia.readString("PartonVertex:setVertex = on");        //this might
    // allow PYTHIA to keep track of partonic space-time information (default
    // was for 'rope hadronization')

    // setting hadron color tags so that spacetime information can be
    // reconstructed
    pythia.readString("StringFragmentation:TraceColours = on");

    // using QCD based color reconnection (original PYTHIA MPI based CR can't be
    // used at hadron level)
    pythia.readString(
        "ColourReconnection:reconnect = off");  // allowing color reconnections
                                                // (should have been default on,
                                                // but doing it here for
                                                // clarity)
    /*pythia.readString("ColourReconnection:mode = 1");                //sets
       the color reconnection scheme to 'new' QCD based scheme (TODO: make sure
       this is better than (2)gluon move)
          pythia.readString("ColourReconnection:forceHadronLevelCR = on");
       //allowing color reconnections for these constructed strings!
                                                                            //a
       few params for the QCD based color reconnection scheme are set below.
          pythia.readString("MultipartonInteractions:pT0Ref = 2.15"); //not sure
       if this is needed for this setup, but is part of the recommended
       'default' pythia.readString("ColourReconnection:allowDoubleJunRem =
       off");  //default on - allows directly connected double junction systems
       to split into two strings
          pythia.readString("ColourReconnection:junctionCorrection = 1.15");
          pythia.readString("ColourReconnection:timeDilationMode = 3"); //allow
       reconnection if single pair of dipoles are in causal contact (maybe try 5
       as well?) pythia.readString("ColourReconnection:timeDilationPar = 0.18");
       //parameter used in causal interaction between strings (mode set
       above)(maybe try 0.073?)
    */

    std::stringstream lines;
    lines << GetXMLElementText({"JetHadronization", "LinesToRead"}, false);
    while (std::getline(lines, s, '\n')) {
      if (s.find_first_not_of(" \t\v\f\r") == s.npos)
        continue;  // skip empty lines
      JSINFO << "Also reading in: " << s;
      pythia.readString(s);
    }

    // optional input of another pythia particle data xml file (higher excited
    // states,...)
    xml_intin =
        GetXMLElementInt({"JetHadronization", "additional_pythia_particles"});
    if (xml_intin == 0 || xml_intin == 1) {
      additional_pythia_particles = xml_intin;
    }
    xml_intin = -1;

    if (additional_pythia_particles == 1) {
      std::string additional_pythia_particle_file = GetXMLElementText(
          {"JetHadronization", "additional_pythia_particles_path"});
      pythia.particleData.readXML(additional_pythia_particle_file, false);
    }

    // And initialize
    pythia.init();

    // reading in info for thermal partons
    inbrick = false;
    brickL = -1.;
    inhydro = false;
    nreusehydro = 1;

    xml_intin = GetXMLElementInt({"Eloss", "Matter", "brick_med"});
    if (xml_intin == 1) {
      inbrick = true;
    }
    xml_intin = -1;
    xml_doublein = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});
    if (inbrick && (xml_doublein >= 0.)) {
      brickL = xml_doublein;
    }
    xml_doublein = -1.;
    xml_intin = GetXMLElementInt({"nReuseHydro"});
    if (xml_intin > 0) {
      nreusehydro = xml_intin;
    }
    xml_intin = -1;
    xml_doublein = GetXMLElementDouble({"Eloss", "deltaT"});
    if (xml_doublein >= 0.) {
      delta_t = xml_doublein;
    }
    xml_doublein = -1.;

    // Check if hydro is used in user xml file
    tinyxml2::XMLElement* elementXML =
        (tinyxml2::XMLElement*)JetScapeXML::Instance()
            ->GetXMLRootUser()
            ->FirstChildElement();
    while (elementXML) {
      std::string elementName = elementXML->Name();
      if (elementName == "Hydro") {
        inhydro = true;
      }
      elementXML = elementXML->NextSiblingElement();
    }  // at this point the hydro could still be a brick -> check this in the
       // following part

    if (inhydro) {
      tinyxml2::XMLElement* element =
          (tinyxml2::XMLElement*)JetScapeXML::Instance()
              ->GetXMLRootUser()
              ->FirstChildElement();
      while (element) {
        std::string elementName = element->Name();
        if (elementName == "Hydro") {
          inhydro = true;
          tinyxml2::XMLElement* childElement =
              (tinyxml2::XMLElement*)element->FirstChildElement();
          while (childElement) {
            std::string childElementName = childElement->Name();
            if (childElementName == "Brick") {
              inbrick = true;
              inhydro = false;
            }
            childElement = childElement->NextSiblingElement();
          }
        }
        element = element->NextSiblingElement();
      }
    }

    // this is only important if a boost invariant 2+1d hydro is used
    xml_doublein =
        GetXMLElementDouble({"JetHadronization", "eta_max_boost_inv"});
    if (inhydro && (xml_doublein >= 0.)) {
      eta_max_boost_inv = xml_doublein;
    }
    xml_doublein = -1.;
  }
}

void HybridHadronization::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("Hadronization Module : " + GetId());
}

// TODO: Junction Strings, Thermal Partons
void HybridHadronization::DoHadronization(
    vector<vector<shared_ptr<Parton>>>& shower,
    vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut) {
  number_p_fake = 0;  // reset counter for the number of fake partons
  double energy_hadrons =
      0.;  // needed for the kinetic scaling of the negative hadrons

  // JSINFO<<"Start Hybrid Hadronization using both Recombination and PYTHIA
  // Lund string model.";
  pythia.event.reset();
  HH_shower.clear();
  parton_collection neg_ptns;

  // pointer to get hydro temperature info
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

  bool boost_invariant = true;
  bool Cartesian_hydro = false;  // not properly initialized by all hydro
                                 // modules at the moment, don't use it
  if (inhydro) {
    std::shared_ptr<FluidDynamics> hydro_ptr =
        JetScapeSignalManager::Instance()->GetHydroPointer().lock();
    const EvolutionHistory& bulk_info = hydro_ptr->get_bulk_info();

    boost_invariant = bulk_info.is_boost_invariant();
    Cartesian_hydro = bulk_info.is_Cartesian();
    Cartesian_hydro = false;
  }

  // negative brickL, then we will not run a brick sampler
  if (brickL < 0.) {
    brickL = 0.;
    inbrick = false;
  }

  // The framework uses unsigned int for partons, while HHpartons have int
  // anti-/color tags. Changing the HH partons to unsigned int is not that easy
  // as some workflow, e.g. in recomb(), needs them to be integer in some
  // structures. Let's have a workaround and set the tags of partons which are
  // too large for an int value to a new unused value.
  convert_color_tags_to_int_type(shower);

  double total_shower_energy = 0.;      // used to scale the negative hadrons
  double total_shower_energy_neg = 0.;  // used to scale the negative hadrons
  // sort positive partons into HH_shower and negative partons from LBT to
  // neg_ptns
  for (unsigned int ishower = 0; ishower < shower.size(); ++ishower) {
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ++ipart) {
      HHparton sh_parton;
      sh_parton.is_shower(true);
      sh_parton.id(shower.at(ishower).at(ipart)->pid());
      sh_parton.orig(0);  // sh_parton.string_id(str)
      sh_parton.px(shower.at(ishower).at(ipart)->px());
      sh_parton.py(shower.at(ishower).at(ipart)->py());
      sh_parton.pz(shower.at(ishower).at(ipart)->pz());
      sh_parton.e(shower.at(ishower).at(ipart)->e());
      sh_parton.x(shower.at(ishower).at(ipart)->x_in().x());
      sh_parton.y(shower.at(ishower).at(ipart)->x_in().y());
      sh_parton.z(shower.at(ishower).at(ipart)->x_in().z());
      sh_parton.x_t(shower.at(ishower).at(ipart)->x_in().t());
      sh_parton.mass(
          sh_parton.e() * sh_parton.e() - sh_parton.px() * sh_parton.px() -
          sh_parton.py() * sh_parton.py() - sh_parton.pz() * sh_parton.pz());
      sh_parton.mass((sh_parton.mass() >= 0.) ? sqrt(sh_parton.mass())
                                              : sqrt(-sh_parton.mass()));
      sh_parton.col(shower.at(ishower).at(ipart)->color());
      sh_parton.acol(shower.at(ishower).at(ipart)->anti_color());

      if (shower.at(ishower).at(ipart)->pstat() > -1) {
        HH_shower.add(sh_parton);
        total_shower_energy += HH_shower[HH_shower.num() - 1].e();
      } else {
        neg_ptns.add(sh_parton);
        total_shower_energy_neg += neg_ptns[neg_ptns.num() - 1].e();
      }
    }
    JSDEBUG << "Shower#" << ishower + 1
            << ". Number of partons to hadronize: " << HH_shower.num();
    JSDEBUG << "Shower#" << ishower + 1
            << ". Number of (neg) partons to hadronize: " << neg_ptns.num();
  }

  // sample thermal partons, brick or hydro hypersurface
  // checking to see if we need to sample thermal partons, if either in brick or
  // hydro
  bool runsampler = (((inhydro || inbrick) && HH_thermal.num() == 0) ||
                     (GetCurrentEvent() % nreusehydro == 0))
                        ? true
                        : false;

  if (runsampler && inbrick) {
    HH_thermal.clear();  // emptying the thermal partons if we're resampling
    ThermalPartonSampler brick(rand_seed, hydro_Tc);  // creating a thermal
                                                      // brick
    brick.brick_length_width(brickL, brickL);
    brick.brick_flow(0., 0., 0.);
    brick.sample_brick();

    JSINFO << "A " << brickL << " fm brick was sampled, generating "
           << brick.nTot() << " partons (" << brick.th_nL() << " light, "
           << brick.th_nS() << " strange).";

    for (int ith = 0; ith < brick.nTot(); ++ith) {
      HHparton thparton;
      thparton.is_thermal(true);
      thparton.id(brick.th_pid(ith));
      thparton.orig(1);  // sh_parton.string_id(str)
      thparton.col(0);
      thparton.acol(0);
      thparton.px(brick.th_px(ith));
      thparton.py(brick.th_py(ith));
      thparton.pz(brick.th_pz(ith));
      thparton.e(brick.th_e(ith));
      thparton.x(brick.th_x(ith));
      thparton.y(brick.th_y(ith));
      thparton.z(brick.th_z(ith));
      thparton.x_t(brick.th_t(ith));
      thparton.mass(
          thparton.e() * thparton.e() - thparton.px() * thparton.px() -
          thparton.py() * thparton.py() - thparton.pz() * thparton.pz());
      thparton.mass((thparton.mass() >= 0.) ? sqrt(thparton.mass())
                                            : sqrt(-thparton.mass()));
      thparton.pos_str(1);
      HH_thermal.add(thparton);  // adding this parton to thermal collection
    }
    // read in thermal partons, THEN do sibling setup...
    for (int i = 0; i < HH_thermal.num(); ++i) {
      HH_thermal[i].sibling(i);
      HH_thermal[i].string_id(-i);
      HH_thermal[i].is_used(true);
      HH_thermal[i].sibling(findthermalsibling(i, HH_thermal));
      HH_thermal[i].is_used(false);
    }
  } else if (runsampler && inhydro) {
    HH_thermal.clear();  // emptying the thermal partons if we're resampling
    std::vector<SurfaceCellInfo> surface_cells;
    GetHydroHyperSurface(hydro_Tc, surface_cells);
    /*std::cout << "\n\n     Surface Cells\n";
    for(int icel=0; icel<surface_cells.size(); ++icel){
        if(icel % int(surface_cells.size()/100) != 0){continue;}
        std::cout << surface_cells[icel].tau << ":" << surface_cells[icel].eta
    << ":" << surface_cells[icel].x << ":" << surface_cells[icel].y << ",   ";
        std::cout << surface_cells[icel].vx << ":" << surface_cells[icel].vy <<
    ":" << surface_cells[icel].vz << ",   "; std::cout <<
    surface_cells[icel].qgp_fraction << ":" << surface_cells[icel].temperature
    << ":" << surface_cells[icel].pressure << "\n";
    }
    std::cout << "\n\n";*/

    // Preallocate memory for surface
    std::vector<std::vector<double>> surface;
    surface.reserve(surface_cells.size());
    for (const auto& cell_info : surface_cells) {
      std::vector<double> cell = {
          cell_info.tau,
          cell_info.x,
          cell_info.y,
          cell_info.eta,
          cell_info.d3sigma_mu[0],
          cell_info.d3sigma_mu[1],
          cell_info.d3sigma_mu[2],
          cell_info.d3sigma_mu[3],
          cell_info.temperature,
          cell_info.umu[1] / cell_info.umu[0],  // vx
          cell_info.umu[2] / cell_info.umu[0],  // vy
          cell_info.umu[3] / cell_info.umu[0]   // vz
      };
      surface.emplace_back(std::move(cell));
    }

    ThermalPartonSampler part_samp(
        rand_seed, hydro_Tc);  // initializing sampler with random seed
    part_samp.set_hypersurface(surface);
    if (boost_invariant) {
      part_samp.sample_2p1d(eta_max_boost_inv);
    } else {
      part_samp.sample_3p1d(Cartesian_hydro);
    }
    JSINFO << "Hydro was sampled, generating " << part_samp.nTot()
           << " partons (" << part_samp.th_nL() << " light, "
           << part_samp.th_nS() << " strange).";

    for (int ith = 0; ith < part_samp.nTot(); ++ith) {
      HHparton thparton;
      thparton.is_thermal(true);
      thparton.id(part_samp.th_pid(ith));
      thparton.orig(1);  // sh_parton.string_id(str)
      thparton.px(part_samp.th_px(ith));
      thparton.py(part_samp.th_py(ith));
      thparton.pz(part_samp.th_pz(ith));
      thparton.e(part_samp.th_e(ith));
      thparton.x(part_samp.th_x(ith));
      thparton.y(part_samp.th_y(ith));
      thparton.z(part_samp.th_z(ith));
      thparton.x_t(part_samp.th_t(ith));
      thparton.mass(
          thparton.e() * thparton.e() - thparton.px() * thparton.px() -
          thparton.py() * thparton.py() - thparton.pz() * thparton.pz());
      thparton.mass((thparton.mass() >= 0.) ? sqrt(thparton.mass())
                                            : sqrt(-thparton.mass()));
      thparton.pos_str(1);
      HH_thermal.add(thparton);  // adding this parton to thermal collection
    }
    // read in thermal partons, THEN do sibling setup...
    for (int i = 0; i < HH_thermal.num(); ++i) {
      HH_thermal[i].sibling(i);
      HH_thermal[i].string_id(-i);
      HH_thermal[i].is_used(true);
      HH_thermal[i].sibling(findthermalsibling(i, HH_thermal));
      HH_thermal[i].is_used(false);
    }
  }

  if (!inbrick && !inhydro) {  // propagate positive and negative partons for
                               // time part_prop in vacuum
    for (int i_sh = 0; i_sh < HH_shower.num(); ++i_sh) {
      double max_t_dif = part_prop;
      double vel[3];
      vel[0] = HH_shower[i_sh].px() / HH_shower[i_sh].e();
      vel[1] = HH_shower[i_sh].py() / HH_shower[i_sh].e();
      vel[2] = HH_shower[i_sh].pz() / HH_shower[i_sh].e();
      HH_shower[i_sh].x(HH_shower[i_sh].x() + vel[0] * max_t_dif);
      HH_shower[i_sh].y(HH_shower[i_sh].y() + vel[1] * max_t_dif);
      HH_shower[i_sh].z(HH_shower[i_sh].z() + vel[2] * max_t_dif);
      HH_shower[i_sh].x_t(HH_shower[i_sh].x_t() + max_t_dif);
    }
    for (int i_sh = 0; i_sh < neg_ptns.num(); ++i_sh) {
      double max_t_dif = part_prop;
      double vel[3];
      vel[0] = neg_ptns[i_sh].px() / neg_ptns[i_sh].e();
      vel[1] = neg_ptns[i_sh].py() / neg_ptns[i_sh].e();
      vel[2] = neg_ptns[i_sh].pz() / neg_ptns[i_sh].e();
      neg_ptns[i_sh].x(neg_ptns[i_sh].x() + vel[0] * max_t_dif);
      neg_ptns[i_sh].y(neg_ptns[i_sh].y() + vel[1] * max_t_dif);
      neg_ptns[i_sh].z(neg_ptns[i_sh].z() + vel[2] * max_t_dif);
      neg_ptns[i_sh].x_t(neg_ptns[i_sh].x_t() + max_t_dif);
    }
  }

  // checking to see if partons are inside medium, and if so then propagate all
  // partons to hypersurface! also can propagate by part_prop time
  for (int i_sh = 0; i_sh < HH_shower.num(); ++i_sh) {
    double vel[3];
    vel[0] = HH_shower[i_sh].px() / HH_shower[i_sh].e();
    vel[1] = HH_shower[i_sh].py() / HH_shower[i_sh].e();
    vel[2] = HH_shower[i_sh].pz() / HH_shower[i_sh].e();
    double t_dif = 0.;
    if (inbrick && HH_shower[i_sh].x_t() < brickL) {
      t_dif = brickL - HH_shower[i_sh].x_t();
    } else if (inhydro) {
      int i = 0;
      double temp = 999999.;
      while (temp > hydro_Tc) {
        double tnow = HH_shower[i_sh].x_t() + delta_t * double(i);
        double xnow = HH_shower[i_sh].x() + vel[0] * delta_t * double(i);
        double ynow = HH_shower[i_sh].y() + vel[1] * delta_t * double(i);
        double znow = HH_shower[i_sh].z() + vel[2] * delta_t * double(i);
        GetHydroCellSignal(tnow, xnow, ynow, znow, check_fluid_info_ptr);
        temp = check_fluid_info_ptr->temperature;
        if (temp <= hydro_Tc) {
          t_dif = delta_t * double(i);
        }
        ++i;
      }
    }
    double max_t_dif = std::max(t_dif, part_prop);
    HH_shower[i_sh].x(HH_shower[i_sh].x() + vel[0] * max_t_dif);
    HH_shower[i_sh].y(HH_shower[i_sh].y() + vel[1] * max_t_dif);
    HH_shower[i_sh].z(HH_shower[i_sh].z() + vel[2] * max_t_dif);
    HH_shower[i_sh].x_t(HH_shower[i_sh].x_t() + max_t_dif);
    // JSINFO<<"Parton propagated to: " << HH_shower[i_sh].x() << ", " <<
    // HH_shower[i_sh].y() << ", " << HH_shower[i_sh].z() << ", " <<
    // HH_shower[i_sh].x_t();
  }
  // do the same for the negative partons
  for (int i_sh = 0; i_sh < neg_ptns.num(); ++i_sh) {
    double vel[3];
    vel[0] = neg_ptns[i_sh].px() / neg_ptns[i_sh].e();
    vel[1] = neg_ptns[i_sh].py() / neg_ptns[i_sh].e();
    vel[2] = neg_ptns[i_sh].pz() / neg_ptns[i_sh].e();
    double t_dif = 0.;
    if (inbrick && neg_ptns[i_sh].x_t() < brickL) {
      t_dif = brickL - neg_ptns[i_sh].x_t();
    } else if (inhydro) {
      int i = 0;
      double temp = 999999.;
      while (temp > hydro_Tc) {
        double tnow = neg_ptns[i_sh].x_t() + delta_t * double(i);
        double xnow = neg_ptns[i_sh].x() + vel[0] * delta_t * double(i);
        double ynow = neg_ptns[i_sh].y() + vel[1] * delta_t * double(i);
        double znow = neg_ptns[i_sh].z() + vel[2] * delta_t * double(i);
        GetHydroCellSignal(tnow, xnow, ynow, znow, check_fluid_info_ptr);
        temp = check_fluid_info_ptr->temperature;
        if (temp <= hydro_Tc) {
          t_dif = delta_t * double(i);
        }
        ++i;
      }
    }
    double max_t_dif = std::max(t_dif, part_prop);
    neg_ptns[i_sh].x(neg_ptns[i_sh].x() + vel[0] * max_t_dif);
    neg_ptns[i_sh].y(neg_ptns[i_sh].y() + vel[1] * max_t_dif);
    neg_ptns[i_sh].z(neg_ptns[i_sh].z() + vel[2] * max_t_dif);
    neg_ptns[i_sh].x_t(neg_ptns[i_sh].x_t() + max_t_dif);
    // JSINFO<<"Parton propagated to: " << neg_ptns[i_sh].x() << ", " <<
    // neg_ptns[i_sh].y() << ", " << neg_ptns[i_sh].z() << ", " <<
    // neg_ptns[i_sh].x_t();
  }

  /*for(int i=0; i<HH_thermal.num(); i++){
    std::cout << HH_thermal[i].id() << "," << HH_thermal[i].x_t() << "," <<
  HH_thermal[i].x()
    << "," << HH_thermal[i].y() << "," << HH_thermal[i].z() << "," <<
  HH_thermal[i].e()
    << "," << HH_thermal[i].px() << "," << HH_thermal[i].py() << "," <<
  HH_thermal[i].pz()
    << "," << HH_thermal[i].mass() << std::endl;
  }*/

  // separately running over positive and negative partons
  for (int pos_ptn = 1; pos_ptn >= 0; --pos_ptn) {
    double tmp_threco =
        th_recofactor;  // need to keep track of th_recofactor (strength) to
                        // disable it for negative partons
    double tmp_maxB_level =
        maxB_level;  // need to keep track of maxB_level to disable baryon
                     // recombination for negative partons
    int tmp_reco_hadrons_pythia = reco_hadrons_pythia;
    if (pos_ptn == 0) {  // handling negative partons by wiping shower partons &
                         // hadrons, then rerunning with negative partons
      pythia.event.reset();
      HH_shower.clear();
      HH_shower = neg_ptns;  // hadrons, remnants, pyremn cleared below

      // if the fragmentation hadrons are given to the afterburner module
      // let the negative partons decay here, as they are not propagated in the
      // afterburner
      if (afterburner_frag_hadrons) {
        reco_hadrons_pythia = 1;
        pythia.readString("HadronLevel:Decay = on");
        pythia.init();
      }

      // add holes left by used thermal partons to HH_shower
      for (int i = 0; i < HH_thermal.num(); ++i) {
        if (HH_thermal[i].is_used()) {
          HH_thermal[i].is_used(false);
          HH_shower.add(HH_thermal[i]);
        }
      }

      // add Extrapartons from recomb() of hadrons to HH_shower of negative
      // partons
      for (int i = 0; i < HH_recomb_extrapartons.num(); ++i) {
        HH_shower.add(HH_recomb_extrapartons[i]);
      }

      if (HH_shower.num() != 0) {
        th_recofactor = 0.;  // killing thermal + negative parton recombination.
        maxB_level = -1;  // disable baryon recombination for negative partons
      }  // do this only if there are negative partons, otherwise the baryon
         // recombination would be switched off after the first event without
         // negative partons
    }
    if (HH_shower.num() == 0) {
      th_recofactor = tmp_threco;
      maxB_level = tmp_maxB_level;
      reco_hadrons_pythia = tmp_reco_hadrons_pythia;
      if (pythia_decays == "off") {
        pythia.readString("HadronLevel:Decay = off");
        pythia.init();
      }
      continue;
    }  // attempting to handle events/configurations with 0 partons will result
       // in a crash

    int attempt_num = 0;
    bool run_successfully = false;
    while ((attempt_num < attempts_max) && (!run_successfully)) {
      HH_showerptns = HH_shower;
      // clearing hadrons and remnants collections
      HH_hadrons.clear();
      HH_remnants.clear();
      HH_pyremn.clear();
      HH_pythia_hadrons.clear();
      HH_recomb_extrapartons.clear();

      // since we 'might' need to reset the thermal partons (if present!)...
      // because we alter the thermal parton collection in the string repair
      // routine
      for (int i = 0; i < HH_thermal.num(); ++i) {
        HH_thermal[i].is_used(false);
        HH_thermal[i].is_decayedglu(false);
        HH_thermal[i].is_remnant(false);
        HH_thermal[i].used_reco(false);
        HH_thermal[i].used_str(false);
        HH_thermal[i].is_thermal(true);
        HH_thermal[i].orig(1);
        HH_thermal[i].par(-1);
        HH_thermal[i].status(0);
        HH_thermal[i].col(0);
        HH_thermal[i].acol(0);
        HH_thermal[i].is_fakeparton(false);
        HH_thermal[i].PY_par1(-1);
        HH_thermal[i].PY_par2(-1);
        HH_thermal[i].PY_dau1(-1);
        HH_thermal[i].PY_dau2(-1);
        HH_thermal[i].PY_stat(23);
        HH_thermal[i].string_id(-i);
        HH_thermal[i].is_strendpt(false);
        HH_thermal[i].pos_str(1);
        HH_thermal[i].endpt_id(0);

        if (pos_ptn == 0) {
          HH_thermal[i].is_used(
              true);  // set all the thermal partons to is used for negative
                      // partons, such that they are not used in the
                      // findthermalsibling routine
        }
      }

      // checking the shower for any color singlet particles to just dump into
      // PYTHIA also handling colored non-partonic particles
      int i_show = 0;
      while (i_show < HH_showerptns.num()) {
        if (HH_showerptns[i_show].id() == 21) {
          ++i_show;
          continue;
        }  // is gluon
        else if (std::abs(HH_showerptns[i_show].id()) <= 6) {
          ++i_show;
          continue;
        }  // is quark
        else if ((std::abs(HH_showerptns[i_show].id()) >= 1103) &&
                 (std::abs(HH_showerptns[i_show].id()) <= 5503) &&
                 ((HH_showerptns[i_show].id() / 10) % 10 == 0)) {
          ++i_show;
          continue;
        }  // is diquark
        else if (pythia.particleData.colType(HH_showerptns[i_show].id()) ==
                 2) {  // this is a non-gluon color octet...
          HH_showerptns[i_show].PY_origid(HH_showerptns[i_show].id());
          HH_showerptns[i_show].id(21);
          ++i_show;
          continue;
        }  // this is a non-(simple)quark color triplet... (pid 42(scalar
           // leptoquark), 4000001(1-6, excited quarks), and SUSY particles)
        else if (std::abs(pythia.particleData.colType(
                     HH_showerptns[i_show].id())) == 1) {
          HH_showerptns[i_show].PY_origid(HH_showerptns[i_show].id());
          HH_showerptns[i_show].id(
              1 * (2 * std::signbit(-pythia.particleData.colType(
                           HH_showerptns[i_show].id())) -
                   1));
          ++i_show;
          continue;
        }  // and the two above should catch 'colored technihadrons' too!
        else if (HH_showerptns[i_show].id() == 90) {
          HH_showerptns.remove(i_show);
          continue;
        }  // PYTHIA specific pid for the 'system' as a whole, shouldn't
           // actually be here.  Dumping it.
        // is none of the above (a colorless object of some sort (a hadron,
        // photon, lepton...))
        HHhadron had;
        had.id(HH_showerptns[i_show].id());
        had.orig(HH_showerptns[i_show].orig());
        had.mass(HH_showerptns[i_show].mass());
        had.pos(HH_showerptns[i_show].pos());
        had.P(HH_showerptns[i_show].P());
        HH_hadrons.add(had);
        HH_showerptns.remove(i_show);
      }

      // setting up the strings appropriately for showers - assumes that color
      // tags are set. if there are colored particles without set color tags, it
      // will dump those partons in particular into a single string
      int num_strings = 0;
      stringform();
      num_strings = HH_showerptns[HH_showerptns.num() - 1].string_id();

      // last attempt in hadronization is done without recombination
      double sh_recofactor_store = sh_recofactor;
      if (attempt_num == attempts_max - 1) {
        sh_recofactor = 0.;
        if (pos_ptn == 1) {
          JSWARN << "Hadronization failed, try without recombination one more "
                    "time";
        } else {
          JSWARN << "Hadronization of negative partons failed, try without "
                    "recombination one more time";
        }
      }

      // running recombination
      recomb();
      sh_recofactor = sh_recofactor_store;

      // function to prepare strings for input into PYTHIA8
      // will need a final reindexing for py_remn (all PY_par#, PY_dau# will
      // need to be += 1) AFTER this function (either here, or in invoke_py)
      // when recursive 'workaround' for large strings is properly
      // handled/removed, then can put this reindexing inside the function
      if (HH_remnants.num()) {
        bool cut = (attempt_num > -1) ? true : false;
        stringprep(HH_remnants, HH_pyremn, cut);
      }

      // temporary workaround to force all formed hadrons to be final state - so
      // pythia won't overwrite spacetime info
      for (int i = 0; i < HH_hadrons.num(); ++i) {
        HH_hadrons[i].is_final(true);
      }

      // if negative fakepartons, we set the momenta and energy to zero before
      // giving them to pythia
      if (pos_ptn == 0) {
        for (int ipart = 0; ipart < HH_pyremn.num(); ipart++) {
          if (HH_pyremn[ipart].is_fakeparton()) {
            if (HH_pyremn[ipart].id() == 21 &&
                HH_pyremn[ipart].mass() >
                    (2. * pythia.particleData.m0(211) + 0.01)) {
              HH_pyremn[ipart].mass(2. * pythia.particleData.m0(211) + 0.01);
            }
            HH_pyremn[ipart].e(HH_pyremn[ipart].mass());
            HH_pyremn[ipart].px(0.);
            HH_pyremn[ipart].py(0.);
            HH_pyremn[ipart].pz(0.);
          }
        }
      }

      // running remaining partons through PYTHIA8 string fragmentation
      run_successfully = invoke_py();

      // for a successful run, go though final hadrons here and set parton
      // parents up
      // bring the hadrons to their corresponding pythia mass shell
      if (run_successfully) {
        // remove all hadrons from HH_hadrons if reco hadrons are put into
        // pythia in this case they are already contained in HH_pythia_hadrons
        if (reco_hadrons_pythia) {
          HH_hadrons.clear();
        }

        // add all fragmentation hadrons from pythia to HH_hadrons, now that
        // pythia hadronization was successfull
        for (int iHad = 0; iHad < HH_pythia_hadrons.num(); iHad++) {
          HH_hadrons.add(HH_pythia_hadrons[iHad]);
        }

        bring_hadrons_to_mass_shell(HH_hadrons);

        // hadron propagation (neglect the photons)
        for (int iHad = 0; iHad < HH_hadrons.num(); ++iHad) {
          if (std::abs(had_prop) > 0.0001 &&
              HH_hadrons[iHad].id() !=
                  22) {  // assumes that hadron will be propagated by more than
                         // 0.0001 fm/c in lab frame - can even propagate
                         // backwards, if that's something wanted...
            double vel[3];
            vel[0] = HH_hadrons[iHad].px() / HH_hadrons[iHad].e();
            vel[1] = HH_hadrons[iHad].py() / HH_hadrons[iHad].e();
            vel[2] = HH_hadrons[iHad].pz() / HH_hadrons[iHad].e();
            double t_prop = had_prop / sqrt(1. - vel[0] * vel[0] -
                                            vel[1] * vel[1] - vel[2] * vel[2]);
            HH_hadrons[iHad].x(HH_hadrons[iHad].x() + vel[0] * t_prop);
            HH_hadrons[iHad].y(HH_hadrons[iHad].y() + vel[1] * t_prop);
            HH_hadrons[iHad].z(HH_hadrons[iHad].z() + vel[2] * t_prop);
            HH_hadrons[iHad].x_t(HH_hadrons[iHad].x_t() + had_prop);
            // JSINFO<<"Hadron propagated to: " << HH_hadrons[iHad].x() << ", "
            // << HH_hadrons[iHad].y() << ", " << HH_hadrons[iHad].z() << ", "
            // << HH_hadrons[iHad].x_t();
          }
        }

        double energy_check = 0.;
        for (int iHad = 0; iHad < HH_hadrons.num(); iHad++) {
          energy_check += HH_hadrons[iHad].e();
        }
        if (pos_ptn == 1 && (energy_check < 0.99 * total_shower_energy)) {
          run_successfully = false;
        }
      }
      ++attempt_num;
    }

    if (!run_successfully) {
      HH_hadrons.clear();
      if (pos_ptn == 1) {
        JSWARN << "This event could not be hadronized.";
      } else {
        JSWARN << "The negative partons of this event could not be hadronized.";
      }
    }

    // scale the negative hadron energies/momenta such that the energy is
    // conserved when subtracting negative hadrons from the positive ones
    if (pos_ptn == 0) {
      if (HH_hadrons.num() > 0) {
        scale_kinematics_negative_hadrons(
            HH_hadrons, total_shower_energy - total_shower_energy_neg,
            energy_hadrons);
      }
    }

    for (unsigned int iHad = 0; iHad < HH_hadrons.num(); ++iHad) {
      if (HH_hadrons[iHad].is_final()) {
        // setting status flag => xy -> x=1 reco, x=2 frag; y=1 sh-sh, y=2
        // sh-th; neg sign means negative hadron
        // all colorless particles get fragmentation status labels
        int stat = (HH_hadrons[iHad].is_recohad()) ? 810 : 820;
        stat += (HH_hadrons[iHad].is_shth()) ? 2 : 1;
        stat *= (pos_ptn == 0) ? -1 : 1;
        int lab = (pos_ptn == 0) ? -1 : 1;
        int idH = HH_hadrons[iHad].id();
        double mH = HH_hadrons[iHad].mass();
        FourVector p(HH_hadrons[iHad].P());
        FourVector x(HH_hadrons[iHad].pos());
        hOut.push_back(
            std::make_shared<Hadron>(Hadron(lab, idH, stat, p, x, mH)));

        if (pos_ptn == 1) {  // used for scaling of negative hadrons
          energy_hadrons += HH_hadrons[iHad].e();
        }
      }
    }
    th_recofactor = tmp_threco;
    maxB_level = tmp_maxB_level;
    reco_hadrons_pythia = tmp_reco_hadrons_pythia;
    if (pythia_decays == "off") {
      pythia.readString("HadronLevel:Decay = off");
      pythia.init();
    }
  }
}

// scale the negative hadron energies/momenta such that the energy is
// conserved when subtracting negative hadrons from the positive ones
void HybridHadronization::scale_kinematics_negative_hadrons(
    hadron_collection& HH_hadrons, double shower_energy,
    double positive_hadrons_energy) {
  double negative_hadrons_energy_initial = 0;
  for (int i = 0; i < HH_hadrons.num(); i++) {
    negative_hadrons_energy_initial += HH_hadrons[i].e();
  }

  double e_scaling_factor = (positive_hadrons_energy - shower_energy) /
                            negative_hadrons_energy_initial;
  if (e_scaling_factor < 0.) {
    JSDEBUG << "e_scaling_factor negative, hadron energy below parton energy: "
            << positive_hadrons_energy - shower_energy << "GeV";
    return;
  }

  for (int i = 0; i < HH_hadrons.num(); i++) {
    double new_energy = HH_hadrons[i].e() * e_scaling_factor;
    if (new_energy > HH_hadrons[i].mass()) {
      HH_hadrons[i].e(new_energy);
      double new_abs_momentum =
          std::sqrt(new_energy * new_energy -
                    HH_hadrons[i].mass() * HH_hadrons[i].mass());
      double p_scaling_factor =
          new_abs_momentum / std::sqrt(HH_hadrons[i].px() * HH_hadrons[i].px() +
                                       HH_hadrons[i].py() * HH_hadrons[i].py() +
                                       HH_hadrons[i].pz() * HH_hadrons[i].pz());
      HH_hadrons[i].px(HH_hadrons[i].px() * p_scaling_factor);
      HH_hadrons[i].py(HH_hadrons[i].py() * p_scaling_factor);
      HH_hadrons[i].pz(HH_hadrons[i].pz() * p_scaling_factor);
    } else {
      HH_hadrons[i].e(HH_hadrons[i].mass());
      HH_hadrons[i].px(0.0);
      HH_hadrons[i].py(0.0);
      HH_hadrons[i].pz(0.0);
    }
  }
}

void HybridHadronization::convert_color_tags_to_int_type(
    vector<vector<shared_ptr<Parton>>>& shower) {
  // check which partons have anti-/color tags larger than the int limit
  std::vector<unsigned int> used_tags;
  // reduce the upper limit -> the maxtag could fail otherwise, if it becomes
  // larger than maximum int
  int max_allowed_tag = std::numeric_limits<int>::max() - 1e6;
  for (unsigned int ishower = 0; ishower < shower.size(); ishower++) {
    for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ipart++) {
      used_tags.push_back(shower.at(ishower).at(ipart)->color());
      used_tags.push_back(shower.at(ishower).at(ipart)->anti_color());
    }
  }
  // remove consecutive duplicates in the sorted vector
  std::sort(used_tags.begin(), used_tags.end());
  used_tags.resize(std::distance(
      used_tags.begin(), std::unique(used_tags.begin(), used_tags.end())));

  // count the number of needed tags
  std::vector<unsigned int> needed_tags_values;
  int number_needed_tags = 0;
  for (unsigned int tag : used_tags) {
    if (tag > max_allowed_tag) {
      needed_tags_values.push_back(tag);
      number_needed_tags++;
    }
  }

  // find the first free tag that is not yet used
  int lower_limit_tag = 100;
  std::vector<int> new_tags;
  while (new_tags.size() < number_needed_tags) {
    if (std::find(used_tags.begin(), used_tags.end(), lower_limit_tag) ==
        used_tags.end()) {
      new_tags.push_back(lower_limit_tag);
    }
    lower_limit_tag++;
  }

  // exchange the tags which are too large by the new ones
  for (int tag = 0; tag < number_needed_tags; tag++) {
    for (unsigned int ishower = 0; ishower < shower.size(); ishower++) {
      for (unsigned int ipart = 0; ipart < shower.at(ishower).size(); ipart++) {
        unsigned int current_tag = needed_tags_values.at(tag);
        if (shower.at(ishower).at(ipart)->color() == current_tag) {
          shower.at(ishower).at(ipart)->set_color(new_tags.at(tag));
        }
        if (shower.at(ishower).at(ipart)->anti_color() == current_tag) {
          shower.at(ishower).at(ipart)->set_anti_color(new_tags.at(tag));
        }
      }
    }
  }
}

// was used in the pygen code to create ordered strings (and assigned string
// ids) from color tags currently has junction functionality commented out until
// I know how to incorporate it within JETSCAPE junction functionality WILL be
// needed in the event that any given string configurations contain them!
void HybridHadronization::stringform() {
  int nstr = 1;

  for (int i = 0; i < HH_showerptns.num(); ++i) {
    if (HH_showerptns[i].string_id() != 0) {
      continue;
    }

    std::vector<int> colors;
    if (HH_showerptns[i].col() > 0) {
      colors.push_back(HH_showerptns[i].col());
    }
    if (HH_showerptns[i].acol() > 0) {
      colors.push_back(HH_showerptns[i].acol());
    }
    HH_showerptns[i].string_id(nstr);
    ++nstr;
    bool newcolor = true;
    while (newcolor) {
      newcolor = false;
      // checking all final particles not in a string for new color id
      // matches...
      int j = 0;
      while (j < HH_showerptns.num()) {
        if (HH_showerptns[j].string_id() != 0) {
          ++j;
          continue;
        }
        for (int k = 0; k < colors.size(); ++k) {
          if (colors[k] == HH_showerptns[j].col() &&
              (HH_showerptns[j].acol() > 0)) {
            colors.push_back(HH_showerptns[j].acol());
            newcolor = true;
            HH_showerptns[j].string_id(HH_showerptns[i].string_id());
            j = -1;
            break;
          } else if (colors[k] == HH_showerptns[j].acol() &&
                     (HH_showerptns[j].col() > 0)) {
            colors.push_back(HH_showerptns[j].col());
            newcolor = true;
            HH_showerptns[j].string_id(HH_showerptns[i].string_id());
            j = -1;
            break;
          } else if (colors[k] == HH_showerptns[j].col() ||
                     colors[k] == HH_showerptns[j].acol()) {
            HH_showerptns[j].string_id(HH_showerptns[i].string_id());
            break;
          }
        }
        ++j;
      }
      /* //this code WILL be needed in the event that any input string
      configurations contain junctions!
      //checking all the junctions in the event to see if there are any new
      color id matches... for(int
      iJun=0;iJun<pythia.event.sizeJunction();++iJun){ bool thisjuncolors =
      false; for(int j=0;j<3;++j){for(int
      k=0;k<colors.size();++k){if(colors[k]==pythia.event.colJunction(iJun,
      j)){thisjuncolors=true;}}} if(!thisjuncolors){continue;} for(int
      j=0;j<3;++j){ bool addcolor = true; for(int
      k=0;k<colors.size();++k){if(colors[k]==pythia.event.colJunction(iJun,j)){addcolor
      = false; break;}}
                      if(addcolor){colors.push_back(pythia.event.colJunction(iJun,j));
      newcolor=true;}
              }
      }
      */
    }
  }

  // now, finding all partons in string 'i'(out of nstr-1 total strings) and
  // reordering them for "string" i=0 final partons - these are noncolored
  // particles that should just be written out as-is (shouldn't actually be
  // present...)
  for (int i = 1; i < nstr; ++i) {
    std::vector<bool> is_used;
    std::vector<int> ptns_new;
    std::vector<int> ptns_strnow;
    for (int j = 0; j < HH_showerptns.num(); ++j) {
      if (HH_showerptns[j].string_id() == i) {
        ptns_strnow.push_back(j);
        is_used.push_back(false);
      }
    }
    // std::vector<int> usedJuns; for(int j=0; j<NumJunctions;
    // ++j){usedJuns.push_back(false);} int numJunsused = 0;

    // find 'a' quark in the event, if there are none, just grab a gluon.
    std::vector<int> stack;
    bool readcolor = true;
    int firstcol = 0;
    for (int j = 0; j < ptns_strnow.size(); ++j) {
      if (!(HH_showerptns[ptns_strnow[j]].id() == 21)) {
        if ((HH_showerptns[ptns_strnow[j]].id() > 0 &&
             HH_showerptns[ptns_strnow[j]].id() <= 6) ||
            (HH_showerptns[ptns_strnow[j]].id() < -6)) {
          stack.push_back(ptns_strnow[j]);
          is_used[j] = true;
          firstcol = HH_showerptns[stack[0]].acol();
          break;
        } else {
          stack.push_back(ptns_strnow[j]);
          is_used[j] = true;
          firstcol = HH_showerptns[stack[0]].col();
          readcolor = false;
          break;
        }
      }
      // gluon loop catch
      else if (j == ptns_strnow.size() - 1) {
        stack.push_back(ptns_strnow[0]);
        is_used[0] = true;
        firstcol = HH_showerptns[stack[0]].acol();
      }
    }
    ptns_new.push_back(stack[0]);

    while (stack.size() > 0) {
      // catching the end of either a string/junction leg or we've circled back
      // on a gluon loop.
      if (readcolor && (HH_showerptns[stack.back()].col() == firstcol)) {
        stack.pop_back();
        continue;
      } else if (!readcolor &&
                 (HH_showerptns[stack.back()].acol() == firstcol)) {
        stack.pop_back();
        continue;
      }

      // keeping track if we've added (a) new parton(s) to the stack
      bool found = false;

      for (int j = 0; j < ptns_strnow.size(); ++j) {
        if (readcolor && !is_used[j] &&
            (HH_showerptns[stack.back()].col() ==
             HH_showerptns[ptns_strnow[j]].acol())) {
          stack[stack.size() - 1] = ptns_strnow[j];
          is_used[j] = true;
          ptns_new.push_back(stack.back());  //--stack.size(); ++stack.size();
          found = true;
          break;
        } else if (!readcolor && !is_used[j] &&
                   (HH_showerptns[stack.back()].acol() ==
                    HH_showerptns[ptns_strnow[j]].col())) {
          stack[stack.size() - 1] = ptns_strnow[j];
          is_used[j] = true;
          ptns_new.push_back(stack.back());  //--stack.size(); ++stack.size();
          found = true;
          break;
        }
      }

      // this code WILL be needed in the event that any input string
      // configurations contain junctions!
      /*if(!found && ((pythia.event.sizeJunction() - numJunsused) > 0)){
                                int iJun(0), iCol(0);
                                if(readcolor){
                                        for(int
         iJ=0;iJ<pythia.event.sizeJunction();++iJ){ if(usedJuns[iJ]){continue;}
                                                for(int
         iC=0;iC<3;++iC){if(pythia.event.colJunction(iJ,iC)==HH_showerptns[stack.back()].col()
         ){iJun=iJ;iCol=iC;found=true;break;}}if(found){break;}
                                        }
                                }
                                else{
                                        for(int
         iJ=0;iJ<pythia.event.sizeJunction();++iJ){ if(usedJuns[iJ]){continue;}
                                                for(int
         iC=0;iC<3;++iC){if(pythia.event.colJunction(iJ,iC)==HH_showerptns[stack.back()].acol()){iJun=iJ;iCol=iC;found=true;break;}}if(found){break;}
                                        }
                                }

                                //find the next two partons here to add to the
         stack if(found){ usedJuns[iJun] = true; ++numJunsused; readcolor =
         !readcolor;

                                        int cols[2]={0,0};
                                        int offset=0; for(int
         j=0;j<3;++j){if(j==iCol){offset=1;continue;}
         cols[j-offset]=pythia.event.colJunction(iJun,j);}
                                        if(cols[0]>cols[1]){int temp=cols[0];
         cols[0]=cols[1]; cols[1]=temp;}

                                        if(readcolor){
                                                for(int
         j=0;j<ptns_strnow.size();++j){if(!is_used[j] &&
         (cols[0]==HH_showerptns[ptns_strnow[j]].acol())){ stack[stack.size()-1]
         = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back());
         //--stack.size(); ++stack.size(); break;
                                                }}
                                                for(int
         j=0;j<ptns_strnow.size();++j){if(!is_used[j] &&
         (cols[1]==HH_showerptns[ptns_strnow[j]].acol())){
                                                        stack.push_back(ptns_strnow[j]);
         is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size();
         ++stack.size(); break;
                                                }}
                                        }
                                        else{
                                                for(int
         j=0;j<ptns_strnow.size();++j){if(!is_used[j] &&
         (cols[0]==HH_showerptns[ptns_strnow[j]].col())){ stack[stack.size()-1]
         = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back());
         //--stack.size(); ++stack.size(); break;
                                                }}
                                                for(int
         j=0;j<ptns_strnow.size();++j){if(!is_used[j] &&
         (cols[1]==HH_showerptns[ptns_strnow[j]].col())){
                                                        stack.push_back(ptns_strnow[j]);
         is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size();
         ++stack.size(); break;
                                                }}
                                        }
                                }
                        }*/
      if (!found) {
        stack.pop_back();
      }
    }

    // ordering this string
    int endpt = 0;
    for (int j = 0; j < ptns_new.size(); ++j) {
      HH_showerptns[ptns_new[j]].pos_str(j);
      HH_showerptns[ptns_new[j]].endpt_id(
          (HH_showerptns[ptns_new[j]].id() == 21) ? 0 : ++endpt);
      HH_showerptns[ptns_new[j]].is_strendpt(
          (HH_showerptns[ptns_new[j]].id() == 21) ? false : true);
    }
  }

  // handling any remaining colored objects that have not had color tags set
  // these are all just thrown into a single string, which may be 'cut' up later
  // on in the string_prep function these will not be 'ordered' within the
  // string nicely...
  int numendpoint = 0;
  int num_str0 = -1;
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    if (HH_showerptns[i].string_id() != 0) {
      continue;
    }
    int iendpoint = 0;
    bool is_endpoint = false;
    if (HH_showerptns[i].id() != 21) {
      iendpoint = ++numendpoint;
      is_endpoint = true;
    }
    HH_showerptns[i].string_id(0);
    HH_showerptns[i].pos_str(++num_str0);
    HH_showerptns[i].endpt_id(iendpoint);
    HH_showerptns[i].is_strendpt(is_endpoint);
  }

  // reordering HH_showerptns - based on string, and position in the string
  // isn't strictly necessary to do this here, as the critical pos_str and
  // endpt_id variables have been set, but is 'nice'
  // std::stable_sort(&HH_showerptns[0],
  // (&HH_showerptns[HH_showerptns.num()-1])+1, strid_compare);
  std::stable_sort(&HH_showerptns[0],
                   (&HH_showerptns[HH_showerptns.num() - 1]) + 1,
                   [](const HHparton& parton1, const HHparton& parton2) {
                     return (parton1.string_id() < parton2.string_id());
                   });
  // now that the list is sorted based on string id, going to sort the partons
  // in each string based on the position of the partons in the string
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    int start, prev_pos, cur_pos, lastfix;
    lastfix = 0;
    if (i == HH_showerptns.num() - 1) {
      lastfix = 1;
    }
    cur_pos = HH_showerptns[i].string_id();
    if (i == 0) {
      prev_pos = HH_showerptns[0].string_id();
      start = 0;
    }
    if (cur_pos != prev_pos || i == HH_showerptns.num() - 1) {
      std::stable_sort(&HH_showerptns[start], &HH_showerptns[i] + lastfix,
                       [](const HHparton& parton1, const HHparton& parton2) {
                         return (parton1.pos_str() < parton2.pos_str());
                       });
      start = i;
      prev_pos = HH_showerptns[i].pos_str();
    }
  }
}

void HybridHadronization::recomb() {
  // parton list for treating thermal siblings for string repair
  parton_collection Extraparton;

  // declaring a few needed values (based on values in init)
  double hbarc2 = hbarc * hbarc;

  // should create a new parton collection here, one with only shower quarks
  // (inc. from gluon splitting) to consider for recombination - then afterwards
  // can reform gluons and output remnants
  parton_collection showerquarks;

  // clearing remnants, in case it hasn't been done before
  HH_remnants.clear();

  // constructing a list of all the strings in the event
  std::vector<int> list_strs;
  // adding the first string to the list
  list_strs.push_back(HH_showerptns[0].string_id());
  // looping over all the partons in the event, and writing each 'unique' string
  // to the list
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    bool str_match = false;
    for (int j = 0; j < list_strs.size(); ++j) {
      if (HH_showerptns[i].string_id() == list_strs[j]) {
        str_match = true;
      }
    }
    if (!str_match) {
      list_strs.push_back(HH_showerptns[i].string_id());
    }
  }

  //*********************************************************************************************************************
  //		Splitting gluons into q-qbar pairs for recombination
  //*********************************************************************************************************************

  // std::cout <<"Below is Color information of all particles in the String
  // (Col, Acol)"<<endl;
  std::vector<int*> ColInfo3;
  for (int i = 0; i < HH_showerptns.num(); i++) {
    int colinfo3[2] = {HH_showerptns[i].col(), HH_showerptns[i].acol()};
    ColInfo3.push_back(colinfo3);
    // std::cout <<" ( "<<HH_showerptns[i].string_id()<<" ,
    // "<<ColInfo3.at(i)[0]<<"   "<<ColInfo3.at(i)[1]<<"  ) ";
  }
  // std::cout <<endl;

  // while running PythiaBrickTest, error{same col tags and acol tag for the
  // different particles, somehow the cause would be in MATTER in determining
  // the color tag} detected, So, If there are problems from the initial
  // structure, correct it based on the color flow.
  int temptag = 1;
  for (int i = 0; i < HH_showerptns.num(); i++) {
    for (int j = i + 1; j < HH_showerptns.num(); j++) {
      if (HH_showerptns[i].col() != 0 &&
          HH_showerptns[i].col() ==
              HH_showerptns[j]
                  .col()) {  // same col tag for different particles detected.
        for (int k = 0; k < HH_showerptns.num(); k++) {
          if (HH_showerptns[i].col() == HH_showerptns[k].acol()) {
            vector<int> Info;  // location of i, location of k, temp color tag
                               // will be saved in this vector.
            Info.push_back(i);
            Info.push_back(k);
            Info.push_back(temptag);

            HH_showerptns[i].col(temptag);
            HH_showerptns[k].acol(
                temptag);  // apply the change for the col, acol pair to
                           // preserve the color flow with least impact to the
                           // final result.
            temptag++;
            Info.clear();
          }
        }
      }
    }
  }

  for (int i = 0; i < HH_showerptns.num(); i++) {
    for (int j = i + 1; j < HH_showerptns.num(); j++) {
      if (HH_showerptns[i].acol() != 0 &&
          HH_showerptns[i].acol() ==
              HH_showerptns[j]
                  .acol()) {  // same acol tag for different particles detected.
        for (int k = 0; k < HH_showerptns.num(); k++) {
          if (HH_showerptns[i].acol() == HH_showerptns[k].col()) {
            vector<int> Info;  // location of i, location of k, temp color tag
                               // will be saved in this vector.
            Info.push_back(i);
            Info.push_back(k);
            Info.push_back(temptag);

            HH_showerptns[i].acol(temptag);
            HH_showerptns[k].col(
                temptag);  // apply the change for the col, acol pair to
                           // preserve the color flow with least impact to the
                           // final result.
            temptag++;
            Info.clear();
          }
        }
      }
    }
  }

  // std::cout <<endl;
  // std::cout <<"Below is Color information of all particles in the Revised
  // String (Col, Acol)"<<endl;
  std::vector<int*> ColInfo2;
  for (int i = 0; i < HH_showerptns.num(); i++) {
    int colinfo2[2] = {HH_showerptns[i].col(), HH_showerptns[i].acol()};
    ColInfo2.push_back(colinfo2);
    // std::cout <<" ( "<<HH_showerptns[i].string_id()<<" ,
    // "<<ColInfo2.at(i)[0]<<"   "<<ColInfo2.at(i)[1]<<"  ) ";
  }
  // std::cout <<endl;

  set_initial_parton_masses(HH_showerptns);

  // starting with shower partons
  for (int i_pt = 0; i_pt < HH_showerptns.num(); ++i_pt) {
    // if parton is a quark, stick it into showerquarks - and set it's parent id
    // to the quark in the original shower
    if ((std::abs(HH_showerptns[i_pt].id()) <= 5) &&
        (HH_showerptns[i_pt].PY_origid() == 0)) {
      showerquarks.add(HH_showerptns.partons[i_pt]);
      showerquarks[showerquarks.num() - 1].par(i_pt);
    }
    // else parton is a gluon, decay into q-qbar and stick those into the
    // collection; set parent id appropriately for both
    else if ((HH_showerptns[i_pt].id() == 21) &&
             (HH_showerptns[i_pt].PY_origid() == 0)) {
      // setting is_decayedglu; not going to set used until one of it's quarks
      // has been used will set the status to -99 here, this will be changed to
      // -1 after first quark, and to 1 after second
      HH_showerptns[i_pt].is_decayedglu(true);
      HH_showerptns[i_pt].status(-99);

      // choosing a gluon mass - if implemented in the future, can (should) read
      // this from the gluon entry itself (or set if necessary) maybe discard
      // gluon if it is under some threshold of mass (eg < pion?) temporarily
      // saving previously set mass here - here's a good place to check if this
      // is even necessary?
      double temp_glumass = HH_showerptns[i_pt].mass();
      if (HH_showerptns[i_pt].mass() < 2. * xmq + 0.001) {
        HH_showerptns[i_pt].mass(2. * xmq + 0.001);
      }

      // gluon decay function reads in the gluon (and the overwritten random
      // mass), and writes the output q-qbar pair to qpair
      parton_collection qpair;
      gluon_decay(HH_showerptns[i_pt], qpair);

      // swapping back original gluon mass
      HH_showerptns[i_pt].mass(temp_glumass);

      // setting the parents of the q-qbar pair to the original gluon (and other
      // vars appropriately here)
      qpair[0].par(i_pt);
      qpair[1].par(i_pt);
      qpair[0].is_shower(true);
      qpair[1].is_shower(true);
      qpair[0].orig(HH_showerptns[i_pt].orig());
      qpair[1].orig(HH_showerptns[i_pt].orig());
      qpair[0].string_id(HH_showerptns[i_pt].string_id());
      qpair[1].string_id(HH_showerptns[i_pt].string_id());
      qpair[0].pos_str(HH_showerptns[i_pt].pos_str());
      qpair[1].pos_str(HH_showerptns[i_pt].pos_str());

      // adding these partons to the collection
      showerquarks.add(qpair);

      // setting up the sibling relations in showerquarks
      showerquarks[showerquarks.num() - 1].sibling(showerquarks.num() - 2);
      showerquarks[showerquarks.num() - 2].sibling(showerquarks.num() - 1);
    }
    // else it's not a quark(u,d,s,c,b) or gluon and we're skipping it.
    else {
      // JSINFO << "\n\nThere is a parton that is not a quark(u,d,s,c,b) or
      // gluon in the input shower.  Skipping parton in recombination.\n\n";
    }
  }
  //*********************************************************************************************************************
  //		Finished splitting gluons
  //*********************************************************************************************************************

  //*********************************************************************************************************************
  //		Beginning of Recombination routine
  //*********************************************************************************************************************

  // randomly permute an integer array with 'n' entries from 1 to n (actually
  // from 0 to n-1) the array will be used to access the i'th element of the
  // showerquarks collection done as showerquarks[perm0_sharray[i]] since we
  // have separate shower and thermal arrays, we will force the first quark to
  // always be from the shower
  //
  // construct perm1 array of size showerquarks.num() and perm2 array of size
  // showerquarks.num() + HH_thermal.num() to access an element from either
  // shower or thermal partons there will be showerquarks.num() positive entries
  // and HH_thermal.num() negative entries (going to exclude 0 in these, just -1
  // i's) thus, perm2 = { 3, -42, -33,...} will access i=2 from showerquarks,
  // then 41 from thermal, then 32 from thermal...
  //
  // lastly, need to bypass scenarios where the quark has been used, or we're
  // trying to use the same quark twice to determine if we're trying to use the
  // same quark twice, check if perm[i]=perm[j] if used=true OR p[i]=p[j]?, then
  // skip this attempt (continue works well for these cases?)

  // constructing permutation arrays; shower quarks are > 0, thermal are < 0
  // perm2 will always access element [std::abs(perm2[i]) - 1]
  int perm1[showerquarks.num()], perm2[showerquarks.num() + HH_thermal.num()];

  // option to either bias recombination, earlier particles attempt to recombine
  // first
  if (torder_reco) {
    // placeholders to sort by time
    std::vector<std::pair<double, int>> tosort1, tosort2;
    for (int i = 0; i < showerquarks.num(); ++i) {
      tosort1.push_back(std::make_pair(showerquarks[i].x_t(), i));
      tosort2.push_back(std::make_pair(showerquarks[i].x_t(), i + 1));
    }
    for (int i = 0; i < HH_thermal.num(); ++i) {
      tosort2.push_back({HH_thermal[i].x_t(), -i - 1});
    }

    // sorting
    std::stable_sort(std::begin(tosort1), std::end(tosort1));
    std::stable_sort(std::begin(tosort2), std::end(tosort2));

    // saving order into permutation arrays
    for (int i = 0; i < tosort1.size(); ++i) {
      perm1[i] = tosort1[i].second;
    }
    for (int i = 0; i < tosort2.size(); ++i) {
      perm2[i] = tosort2[i].second;
    }
  } else {
    // prepping permutation arrays to be shuffled
    for (int i = 0; i < showerquarks.num(); ++i) {
      perm1[i] = i;
    }
    for (int i = 0; i < HH_thermal.num(); ++i) {
      perm2[i] = (i - HH_thermal.num());
    }
    for (int i = HH_thermal.num(); i < showerquarks.num() + HH_thermal.num();
         ++i) {
      perm2[i] = (i - HH_thermal.num() + 1);
    }

    // permuting these arrays using Fisher-Yates algorithm
    //-- To shuffle an array a of n elements (indices 0..n-1):
    // for i from 0 to n−2 do
    //   j ← random integer such that i ≤ j < n
    //   exchange a[i] and a[j]
    for (int i = 0; i < showerquarks.num() - 1; ++i) {
      int ranelement = i + floor((showerquarks.num() - i) * ran());
      int temp = perm1[i];
      perm1[i] = perm1[ranelement];
      perm1[ranelement] = temp;
    }
    for (int i = 0; i < showerquarks.num() + HH_thermal.num() - 1; ++i) {
      int ranelement =
          i + floor((showerquarks.num() + HH_thermal.num() - i) * ran());
      int temp = perm2[i];
      perm2[i] = perm2[ranelement];
      perm2[ranelement] = temp;
    }
  }

  // std::cout <<"Below is Color information of all particles in the String
  // (Col, Acol)"<<endl;
  std::vector<int*> ColInfo;
  for (int i = 0; i < HH_showerptns.num(); i++) {
    int colinfo[2] = {HH_showerptns[i].col(), HH_showerptns[i].acol()};
    ColInfo.push_back(colinfo);
    // std::cout <<" ( "<<HH_showerptns[i].string_id()<<" ,
    // "<<ColInfo.at(i)[0]<<"   "<<ColInfo.at(i)[1]<<"  ) ";
  }
  // std::cout <<endl;

  // make all color tags to be indices for the matrix (so the indice represent
  // the color tag, don't need to classify anti or not)
  std::vector<vector<int>>
      IndiceForCol1;  // this vector will form vector of vectors(location in the
                      // string, anti or not, color tag)
  std::vector<int>
      IndiceForCol2;  // this vector is one component of vector above
  std::vector<int> IndiceForColFin;  // this vector is final vector that'll be
                                     // used for process

  for (int iIndice = 0; iIndice < HH_showerptns.num(); iIndice++) {
    IndiceForCol2.push_back(iIndice);
    IndiceForCol2.push_back(1);
    IndiceForCol2.push_back(HH_showerptns[iIndice].col());
    IndiceForCol1.push_back(IndiceForCol2);
    IndiceForCol2.clear();

    IndiceForCol2.push_back(iIndice);
    IndiceForCol2.push_back(-1);
    IndiceForCol2.push_back(HH_showerptns[iIndice].acol());
    IndiceForCol1.push_back(IndiceForCol2);
    IndiceForCol2.clear();
  }  // As a result, two vectors with 3 component keep being added to bigger
     // vector, now exclude same tag and zero
  // dignostic measure
  // std::cout <<endl<<" Below is list of all color tags " <<endl;
  // std::cout <<" ( ";
  // for (int icheck=0; icheck < IndiceForCol1.size(); icheck++) {
  // std::cout << IndiceForCol1.at(icheck).at(2) << " , ";
  //}
  // std::cout <<" ) " <<endl;

  vector<int> tempcol;
  // now put all the values in ColInfo into IndiceForColFin
  for (int icol = 0; icol < IndiceForCol1.size(); icol++) {
    tempcol.push_back(IndiceForCol1.at(icol).at(2));
  }
  if (tempcol.at(0) != 0) {
    IndiceForColFin.push_back(tempcol.at(0));
  }
  if (tempcol.at(1) != 0) {
    IndiceForColFin.push_back(tempcol.at(1));
  }

  for (int iclear1 = 0; iclear1 < tempcol.size(); iclear1++) {
    bool foundsame = false;
    for (int iclear2 = 0; iclear2 < IndiceForColFin.size(); iclear2++) {
      if (tempcol.at(iclear1) == IndiceForColFin.at(iclear2) ||
          tempcol.at(iclear1) == 0) {
        foundsame = true;
      }
    }
    if (!foundsame) {
      IndiceForColFin.push_back(tempcol.at(iclear1));
    }
    foundsame = false;
  }

  // last checking whether there is zero in the color tag list
  for (int i = 0; i < IndiceForColFin.size(); i++) {
    if (IndiceForColFin[i] == 0) {
      std::vector<int>::iterator i1 =
          std::find(IndiceForColFin.begin(), IndiceForColFin.end(), 0);
      int distance = std::distance(IndiceForColFin.begin(), i1);
      // std::cout <<endl<<" zero is located at "<<distance<<endl;
      IndiceForColFin.erase(i1);
    }
  }

  // dignostic syntax
  /*std::cout <<endl<<" Below is list of valid color tags " <<endl;
  std::cout <<" ( ";
  for (int icheck=0; icheck < IndiceForColFin.size(); icheck++) {
    std::cout << IndiceForColFin.at(icheck) << " , ";
  }
  std::cout <<" ) " <<endl;*/

  // when we includes the partons from LBT, color tags from them are both zero
  // in col and acol tags, Therefore, based on the maximum color tag from the
  // parton with col , acol tags, reassign the color tags for LBT partons. As a
  // result, one fake string should be formed based on these color Tags First,
  // check the maximum color tags in the vector of valid color tags
  int limit;
  int maxtag;
  maxtag = *max_element(
      IndiceForColFin.begin(),
      IndiceForColFin
          .end());  // this will be so important in dealing thermal partons,
                    // beacause they get color tag incremented from this max tag
  limit = maxtag;  // save initial maxtag for the future usage{It should be used
                   // for color reconnection matrix}
  // std::cout <<endl<<" max tag is "<<maxtag<<endl;

  // before determine the matrix, need to link color tag with the matrix
  // indices( by color tag and trace back to the particle, then determine the
  // probability) It's done above(IndiceForCol1) to define the probability,
  // parameters for distance should exist. possibly, that could be the par() or
  // string number the last
  std::vector<vector<double>>
      MesonrecoMatrix1;  // vector of double (prob) , ,,,// in a big scheme,
                         // these two vector form a Matrix
  std::vector<double> MesonrecoMatrix2;  // this is 1d vector of (probabilty)

  // first, Form Meson recombination Matrix
  for (int irow = 0; irow < IndiceForColFin.size(); irow++) {
    for (int icol = 0; icol < IndiceForColFin.size(); icol++) {
      double Mesonfactor;
      int distance;  // distance between the color tags in string(not c++
                     // function)
      int tag1 = IndiceForColFin.at(irow);
      int tag2 = IndiceForColFin.at(icol);
      std::vector<int>::iterator it1 =
          std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
      std::vector<int>::iterator it2 =
          std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                    tag2);  // set up for finding the location of co tag in
                            // vector that will be used to find distance
      int pos1 = std::distance(IndiceForColFin.begin(),
                               it1);  // location of col tag of irow
      int pos2 = std::distance(IndiceForColFin.begin(),
                               it2);  // location of col tag of icol
      // std::cout <<endl<<"position of "<<IndiceForColFin.at(irow)<<" is " <<
      // pos1 <<endl; std::cout <<endl<<"position of
      // "<<IndiceForColFin.at(icol)<<" is " << pos2 <<endl;

      distance = abs(pos1 - pos2);
      // std::cout <<endl<<" At ( "<<irow<<","<<icol<<" ) , distance is
      // "<<distance<<endl;
      if (distance == 0) {
        Mesonfactor = 1;
        // std::cout <<endl<<"Temp Meson Factor is "<<Mesonfactor<<endl;
        MesonrecoMatrix2.push_back(Mesonfactor);
      };
      if (distance > 0) {
        Mesonfactor = 0.111;
        // std::cout <<endl<<"Temp Meson Factor is "<<Mesonfactor<<endl;
        MesonrecoMatrix2.push_back(Mesonfactor);
      };
    }
    MesonrecoMatrix1.push_back(MesonrecoMatrix2);
    MesonrecoMatrix2.clear();
  }

  // now correct the component in MesonMatrix1 by searching through
  // HH_showerptns{checking gluon}
  for (int icheck1 = 0; icheck1 < HH_showerptns.num(); icheck1++) {
    if (HH_showerptns[icheck1].col() != 0 &&
        HH_showerptns[icheck1].acol() !=
            0) {  // finding gluon means that finding partons that can't form
                  // color singlet, but octet
      int tag1 = HH_showerptns[icheck1].col();
      int tag2 = HH_showerptns[icheck1].acol();
      std::vector<int>::iterator it1 =
          std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
      std::vector<int>::iterator it2 =
          std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                    tag2);  // set up for finding the location of co tag in
                            // vector that will be used to find distance
      int pos1 = std::distance(IndiceForColFin.begin(),
                               it1);  // location of col tag of irow
      int pos2 = std::distance(IndiceForColFin.begin(),
                               it2);  // location of col tag of icol
      MesonrecoMatrix1.at(pos1).at(pos2) = 0;
      MesonrecoMatrix1.at(pos2).at(pos1) = 0;
    }
  }
  /*
  //dignostic measure
  std::cout <<endl<<"Meson reco Matrix is same as below"<<endl;
  for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
      std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
    }
    std::cout <<endl<<endl;
  }
  */
  // MesonrecoMatrix is Formed so far, Now make that of Baryon!!

  std::vector<vector<vector<double>>>
      BaryonrecoMatrix1;  // vector of double (prob) , ,,,// in a big scheme,
                          // these three items form a Matrix of vector
  std::vector<vector<double>>
      BaryonrecoMatrix2;  // this is vector of 2d vector below (probabilty,
                          // required color charge for color neutrality), ( , )
                          // , ( , ) ...
  std::vector<double>
      BaryonrecoMatrix3;  // this is 2d vector of (probabilty, required color
                          // charge for color neutrality)

  for (int irow = 0; irow < IndiceForColFin.size(); irow++) {
    for (int icol = 0; icol < IndiceForColFin.size(); icol++) {
      double Baryonfactor;
      BaryonrecoMatrix3.push_back(
          1. / 27.);  // 1/27: probability to get singlet state
      BaryonrecoMatrix3.push_back(
          0.);  // now the probability factor and color(not specified yet, it'll
                // be done with baryon formation) for neutrality are saved.

      BaryonrecoMatrix2.push_back(BaryonrecoMatrix3);
      BaryonrecoMatrix3.clear();
    }
    BaryonrecoMatrix1.push_back(BaryonrecoMatrix2);
    BaryonrecoMatrix2.clear();
  }
  /*
  //dignostic measure
  std::cout <<endl<<"Baryon reco Matrix is same as below"<<endl;
  for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
      std::cout <<"  ( "<<BaryonrecoMatrix1.at(irow).at(icol).at(1)<<" )  ";
    }
    std::cout <<endl<<endl;
  }
  //std::cout <<endl<<endl;
  */

  // looping over all the quarks that we have...
  //'q1' loops over all quarks in the shower
  //'q2' loops over all quarks in the event
  //'q3' loops over all quarks in the event, starting from 'q2' and ending at
  // the last quark when 'q2' is at the last quark, we will not consider quark
  // 'q3' - can only make a meson at that point...

  parton_collection considering;
  int element[3];

  for (int q1 = 0; q1 < showerquarks.num(); ++q1) {
    // accessing first considered quark
    // set q1 variables here

    if (sh_recofactor < 0.0000000001) {
      continue;
    }  // turning off recombination in a computationally friendly way...

    // resetting madehadron flag, now that we're going to be considering a new
    // hadron
    bool madehadron = false;

    // taking the id from the permutation array, and turning it into quark array
    // element
    element[0] = perm1[q1];

    // skipping if current quark has been used or isn't a u,d,s,c,b quark or
    // antiquark
    if (showerquarks[element[0]].status() != 0 ||
        showerquarks[element[0]].is_used()) {
      continue;
    } else if (std::abs(showerquarks[element[0]].id()) > 5) {
      JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for "
                "recombination, THIS SHOULD NOT HAPPEN!";
      continue;
    }

    // assigning quark values to considering - and setting the status to -991
    // if at the end of the check, we make a hadron, we will set all -99*
    // entries to 1, else we'll set back to 0

    considering.add(showerquarks[element[0]]);
    showerquarks[element[0]].status(-991);

    for (int q2 = 0; q2 < showerquarks.num() + HH_thermal.num(); ++q2) {
      // set q2 variables here - if we can form a meson, then skip q3 loop
      // also skip q3 loop if q2 is at last quark

      double recofactor2 = 1. / 9.;

      // accessing the second considered quark
      // this will skip over non-quark entries in HH_thermal
      if (perm2[q2] > 0) { /*is shower quark*/
        element[1] = perm2[q2] - 1;
      } else { /*is thermal quark*/
        element[1] = perm2[q2] + 1;
      }

      // checking to see if quark is in shower or thermal
      // only need to check if quark is the same as the previous quark IFF it is
      // in the shower only need to check if quark is from the same gluon if it
      // is from a gluon decay in the shower... want to check if we've messed up
      // and are accessing a gluon if it is in original shower/thermal need to
      // check if used for all cases...
      if (perm2[q2] > 0) {
        // skipping if current quark has been used
        if (showerquarks[element[1]].status() != 0 ||
            showerquarks[element[1]].is_used()) {
          continue;
        }
        // skipping if current quark is the same as q1
        else if (element[0] == element[1]) {
          continue;
        }
        // skipping if the current quark is from the same gluon as q1
        else if ((showerquarks[element[0]].par() != -1) &&
                 (showerquarks[element[0]].par() ==
                  showerquarks[element[1]].par())) {
          continue;
        }
        // skipping if the current quark is not in the same string as q1 (both
        // are shower partons) else if(showerquarks[element[1]].string_id() !=
        // showerquarks[element[0]].string_id()){continue;} skipping if current
        // quark is not a u,d,s,c,b quark
        else if (std::abs(showerquarks[element[1]].id()) > 5) {
          JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for "
                    "recombination, THIS SHOULD NOT HAPPEN!";
          continue;
        }

        considering.add(showerquarks[element[1]]);
        showerquarks[element[1]].status(-992);
      } else if (perm2[q2] < 0) {
        // skipping if current quark has been used
        if (HH_thermal[-element[1]].status() != 0 ||
            HH_thermal[-element[1]].is_used()) {
          continue;
        }
        // skipping if current quark is not a u,d,s,c,b quark
        else if (std::abs(HH_thermal[-element[1]].id()) > 5) {
          JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for "
                    "recombination, THIS SHOULD NOT HAPPEN!";
          continue;
        }
        // turning off recombination for thermal partons in a computationally
        // friendly way...
        else if (th_recofactor < 0.001) {
          continue;
        }

        // checking distance cut ONLY if this is a thermal parton - skip if
        // dist2 > dist2cut
        FourVector pos_ptn1 = considering[0].pos();
        FourVector pos_ptn2 = HH_thermal[-element[1]].pos();
        double dt = pos_ptn1.t() - pos_ptn2.t();
        if (dt > 0.) {
          double dt_E =
              dt / HH_thermal[-element[1]].e();  // P/E * dT = dist = P*(dT/E)
          pos_ptn2.Set(pos_ptn2.x() + HH_thermal[-element[1]].px() * dt_E,
                       pos_ptn2.y() + HH_thermal[-element[1]].py() * dt_E,
                       pos_ptn2.z() + HH_thermal[-element[1]].pz() * dt_E, 0.);
        } else {
          double dt_E = -dt / considering[0].e();
          pos_ptn1.Set(pos_ptn1.x() + considering[0].px() * dt_E,
                       pos_ptn1.y() + considering[0].py() * dt_E,
                       pos_ptn1.z() + considering[0].pz() * dt_E, 0.);
        }
        if (dif2(pos_ptn1, pos_ptn2) > dist2cut) {
          continue;
        }

        considering.add(HH_thermal[-element[1]]);
        HH_thermal[-element[1]].status(-992);
      } else {
        JSWARN << "SOMETHING WENT HORRIBLY WRONG - DO NOT KNOW WHERE CURRENT "
                  "QUARK CAME FROM?!";
      }

      // now that we have two 'acceptable' quarks to check hadron formation
      // against, there is no reason to bother checking if we can make a baryon
      // if we have a q-qbar at this point... will skip third loop in this case
      // - otherwise we will check if we can make a baryon...
      if ((considering[0].id() * considering[1].id() > 0) &&
          (q2 < showerquarks.num() + HH_thermal.num() - 1) &&
          (maxB_level > -1)) {
        for (int q3 = q2 + 1; q3 < showerquarks.num() + HH_thermal.num();
             ++q3) {
          double recofactor3 = 2. / 27.;

          // removing all but the first two entries in the considering
          // collection... this should have been done before, but have this here
          // just in case - remove if this doesn't ever trigger
          while (considering.num() > 2) {
            considering.partons.pop_back();
          }

          // accessing the third considered quark
          // this will skip over non-quark entries in HH_thermal
          if (perm2[q3] > 0) { /*is shower quark*/
            element[2] = perm2[q3] - 1;
          } else { /*is thermal quark*/
            element[2] = perm2[q3] + 1;
          }

          // now that we have q3, we need to check if it is valid:
          // q3 needs to be checked if used (all cases)
          // q3 needs to be checked if it is a erroneously accessed parton (not
          // u,d,s,c,b quark) q3 needs to be checked against q1 to see if it is
          // the same, or from the same gluon q3 does not need to be checked
          // against q2 to see if it is the same q3 does need to be checked to
          // see if it is from the same gluon as q2 q3 needs to be checked to
          // make sure that it can form a baryon with q1 and q2 (checking
          // against either is ok) check:  used, sameq1, sameg_q1, sameg_q2,
          // isglu
          if (perm2[q3] > 0) {
            // skipping if the current quark cannot make a baryon with other two
            // quarks
            if (showerquarks[element[2]].id() * considering[0].id() < 0) {
              continue;
            }
            // skipping if current quark has been used
            if (showerquarks[element[2]].status() != 0 ||
                showerquarks[element[2]].is_used()) {
              continue;
            }
            // skipping if current quark is the same as q1
            else if (element[0] == element[2] || element[1] == element[2]) {
              continue;
            }
            // skipping if the current quark is from the same gluon as q1
            else if ((showerquarks[element[0]].par() != -1) &&
                     (showerquarks[element[0]].par() ==
                      showerquarks[element[2]].par())) {
              continue;
            }
            // skipping if the current quark is from the same gluon as q2
            else if ((perm2[q2] > 0) &&
                     (showerquarks[element[1]].par() != -1) &&
                     (showerquarks[element[1]].par() ==
                      showerquarks[element[2]].par())) {
              continue;
            }
            // skipping if the current quark is not in the same string as q1
            //(q2 MUST be in the same string as q1 if it's in the shower, and
            // doesn't need to be checked if thermal) else
            // if(showerquarks[element[2]].string_id() !=
            // showerquarks[element[0]].string_id()){continue;} skipping if
            // current quark is not a u,d,s,c,b quark
            else if (std::abs(showerquarks[element[2]].id()) > 5) {
              JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for "
                        "recombination, THIS SHOULD NOT HAPPEN!";
              continue;
            }

            considering.add(showerquarks[element[2]]);
            showerquarks[element[2]].status(-993);

            // TODO: Here we need to establish the probability about baryon
            // formation. first, need to find 3sets of col tag pair are included
            // junction info at the same time. if all 3 does, prob is 1/27, if 1
            // does, find other component's tag with element[2](3rd one)'s in
            // Meson Matrix! second, scan all junction lists before that, we
            // need to evaluate temporary junction, because some junctions could
            // be eliminated, So need to make tempjunction list. then we need to
            // proocedure to erase the factor in the vectors! Now evaluating the
            // Tempjunctions element first, before this, declare bool variables
            // for the next proocedure
            bool element1 = false;
            bool element2 = false;
            bool element3 = false;
            int juncnum1 = 999999999;
            int juncnum2 = 999999999;
            int juncnum3 = 999999999;

            int standard =
                IndiceForColFin.size();  // to apply the "if" statement, set the
                                         // condition for later syntax
            int tagformatrix =
                999999999;  // if all three tags are in same junction, we need
                            // to set large number not to confuse, this will be
                            // filtered by later syntax
            int loc1 = standard;  // location in the MesonrecoMatrix
            std::vector<int>::iterator I1;
            std::vector<int>::iterator I2;
            int loc2 = standard;

            if (considering[0].id() * considering[1].id() *
                    considering[2].id() >
                0) {
              // baryon case, and this will check whether two particles are in
              // same junction and the other is not. and evaluate(by
              // mesonrecomatrix) the other particle's color tag with the the
              // color tag of 3rd particle collected to be baryon.
              for (int ijunc = 0; ijunc < Tempjunctions.size(); ijunc++) {
                if ((considering[0].col() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[0].col() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[0].col() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element1 = true;
                  juncnum1 = ijunc;
                }
                if ((considering[1].col() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[1].col() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[1].col() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element2 = true;
                  juncnum2 = ijunc;
                }
                if ((considering[2].col() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[2].col() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[2].col() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element3 = true;
                  juncnum3 = ijunc;
                }
              }
              if ((juncnum1 == juncnum2) && (juncnum1 != juncnum3) &&
                  juncnum1 != 999999999 && considering[2].col() != 0) {
                tagformatrix = Tempjunctions.at(juncnum1).at(3).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[2].col());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag1 is
                // "<<HH_showerptns[showerquarks[element[2]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };
              if ((juncnum2 == juncnum3) && (juncnum2 != juncnum1) &&
                  juncnum2 != 999999999 && considering[0].col() != 0) {
                tagformatrix = Tempjunctions.at(juncnum2).at(1).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[0].col());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag2 is
                // "<<HH_showerptns[showerquarks[element[0]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };
              if ((juncnum1 == juncnum3) && (juncnum1 != juncnum2) &&
                  juncnum3 != 999999999 && considering[1].col() != 0) {
                tagformatrix = Tempjunctions.at(juncnum3).at(2).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[1].col());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag3 is
                // "<<HH_showerptns[showerquarks[element[1]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };

              I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                             tagformatrix);
              loc1 = std::distance(IndiceForColFin.begin(), I1);
              // std::cout <<endl<<"chosen orginal color tag is "<<tagformatrix;
              // std::cout <<endl<<"and corresponding indice in the matrix is
              // "<<loc1<<endl;// now we find the locations of the color tags.

              // std::cout <<endl<<"loc2 is given as "<<loc2<<endl;
              //  now, search the mesonrecomatrix
              if (tagformatrix < 999999999 && loc1 < standard &&
                  loc2 <
                      standard) {  // Not all the tags are in same junction so
                                   // loc1,2 has always same or smaller value
                                   // than the size of INdiceForColFin vector
                if (MesonrecoMatrix1.at(loc1).at(loc2) == 1) {
                  recofactor3 = 1;  // one specific tag, which is located in
                                    // different location from other two strings
                                    // from a junction, is neighbored with the
                                    // tag in the not used tag in junction
                } else {
                  recofactor3 = 1. / 9.;
                }  // general case
              } else {
                recofactor3 = 1. / 27.;
              }  // general case

              if ((juncnum1 == juncnum2) && (juncnum2 == juncnum3) &&
                  (tagformatrix != 999999999)) {
                recofactor3 = 1.;
              }
            } else if (considering[0].id() * considering[1].id() *
                           considering[2].id() <
                       0) {
              for (int ijunc = 0; ijunc < Tempjunctions.size(); ijunc++) {
                if ((considering[0].acol() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[0].acol() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[0].acol() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element1 = true;
                  juncnum1 = ijunc;
                }
                if ((considering[1].acol() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[1].acol() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[1].acol() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element2 = true;
                  juncnum2 = ijunc;
                }
                if ((considering[2].acol() ==
                     Tempjunctions.at(ijunc).at(1).at(1)) ||
                    (considering[2].acol() ==
                     Tempjunctions.at(ijunc).at(2).at(1)) ||
                    (considering[2].acol() ==
                     Tempjunctions.at(ijunc).at(3).at(1))) {
                  element3 = true;
                  juncnum3 = ijunc;
                }
              }
              if ((juncnum1 == juncnum2) && (juncnum1 != juncnum3) &&
                  juncnum1 != 999999999) {
                tagformatrix = Tempjunctions.at(juncnum1).at(3).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[2].acol());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag1 is
                // "<<HH_showerptns[showerquarks[element[2]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };
              if ((juncnum2 == juncnum3) && (juncnum2 != juncnum1) &&
                  juncnum2 != 999999999) {
                tagformatrix = Tempjunctions.at(juncnum2).at(1).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[0].acol());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag2 is
                // "<<HH_showerptns[showerquarks[element[0]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };
              if ((juncnum1 == juncnum3) && (juncnum1 != juncnum2) &&
                  juncnum3 != 999999999) {
                tagformatrix = Tempjunctions.at(juncnum3).at(2).at(1);
                I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                               considering[1].acol());
                loc2 = std::distance(IndiceForColFin.begin(), I2);
                // std::cout <<endl<<"chosen color tag3 is
                // "<<HH_showerptns[showerquarks[element[1]].par()].col();
                // std::cout <<endl<<"and corresponding indice in the matrix is
                // "<<loc2<<endl;
              };

              I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(),
                             tagformatrix);
              loc1 = std::distance(IndiceForColFin.begin(), I1);
              // std::cout <<endl<<"chosen orginal color tag is "<<tagformatrix;
              // std::cout <<endl<<"and corresponding indice in the matrix is
              // "<<loc1<<endl;// now we find the locations of the color tags.

              // std::cout <<endl<<"loc2 is given as "<<loc2<<endl;
              //  now, search the mesonrecomatrix
              if (tagformatrix < 999999999 && loc1 < standard &&
                  loc2 <
                      standard) {  // Not all the tags are in same junction so
                                   // loc1,2 has always same or smaller value
                                   // than the size of INdiceForColFin vector
                if (MesonrecoMatrix1.at(loc1).at(loc2) == 1) {
                  recofactor3 = 1;  // one specific tag, which is located in
                                    // different location from other two strings
                                    // from a junction, is neighbored with the
                                    // tag in the not used tag in junction
                } else {
                  recofactor3 = 1. / 9.;
                }  // general case
              } else {
                recofactor3 = 1. / 27.;
              }  // general case

              if ((juncnum1 == juncnum2) && (juncnum2 == juncnum3) &&
                  (tagformatrix != 999999999)) {
                recofactor3 = 1.;
              }
            }
            // recofactor3 = sh_recofactor*recofactor2;
          } else if (perm2[q3] < 0) {
            // skipping if the current quark cannot make a baryon with other two
            // quarks
            if (HH_thermal[-element[2]].id() * considering[0].id() < 0) {
              continue;
            }
            // skipping if current quark has been used
            if (HH_thermal[-element[2]].status() != 0 ||
                HH_thermal[-element[2]].is_used()) {
              continue;
            }
            // skipping if the current quark is from the same gluon as q2
            else if ((perm2[q2] < 0) && (HH_thermal[-element[1]].par() != -1) &&
                     (HH_thermal[-element[1]].par() ==
                      HH_thermal[-element[2]].par())) {
              continue;
            }
            // skipping if current quark is not a u,d,s,c,b quark
            else if (std::abs(HH_thermal[-element[2]].id()) > 5) {
              JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for "
                        "recombination, THIS SHOULD NOT HAPPEN!";
              continue;
            }
            // turning off recombination for thermal partons in a
            // computationally friendly way...
            else if (th_recofactor < 0.001) {
              continue;
            }

            // checking distance cut ONLY if this is a thermal parton - skip if
            // dist2 > dist2cut
            FourVector pos_ptn1 = considering[0].pos();
            FourVector pos_ptn2 = considering[1].pos();
            FourVector pos_ptn3 = HH_thermal[-element[2]].pos();
            if ((pos_ptn1.t() > pos_ptn2.t()) &&
                (pos_ptn1.t() > pos_ptn3.t())) {
              double dt_E2 = (pos_ptn1.t() - pos_ptn2.t()) / considering[1].e();
              double dt_E3 =
                  (pos_ptn1.t() - pos_ptn3.t()) / HH_thermal[-element[2]].e();
              pos_ptn2.Set(pos_ptn2.x() + considering[1].px() * dt_E2,
                           pos_ptn2.y() + considering[1].py() * dt_E2,
                           pos_ptn2.z() + considering[1].pz() * dt_E2, 0.);
              pos_ptn3.Set(pos_ptn3.x() + HH_thermal[-element[2]].px() * dt_E3,
                           pos_ptn3.y() + HH_thermal[-element[2]].py() * dt_E3,
                           pos_ptn3.z() + HH_thermal[-element[2]].pz() * dt_E3,
                           0.);
            } else if ((pos_ptn2.t() > pos_ptn1.t()) &&
                       (pos_ptn2.t() > pos_ptn3.t())) {
              double dt_E1 = (pos_ptn2.t() - pos_ptn1.t()) / considering[0].e();
              double dt_E3 =
                  (pos_ptn2.t() - pos_ptn3.t()) / HH_thermal[-element[2]].e();
              pos_ptn1.Set(pos_ptn1.x() + considering[0].px() * dt_E1,
                           pos_ptn1.y() + considering[0].py() * dt_E1,
                           pos_ptn1.z() + considering[0].pz() * dt_E1, 0.);
              pos_ptn3.Set(pos_ptn3.x() + HH_thermal[-element[2]].px() * dt_E3,
                           pos_ptn3.y() + HH_thermal[-element[2]].py() * dt_E3,
                           pos_ptn3.z() + HH_thermal[-element[2]].pz() * dt_E3,
                           0.);
            } else {
              double dt_E1 = (pos_ptn3.t() - pos_ptn1.t()) / considering[0].e();
              double dt_E2 = (pos_ptn3.t() - pos_ptn2.t()) / considering[1].e();
              pos_ptn1.Set(pos_ptn1.x() + considering[0].px() * dt_E1,
                           pos_ptn1.y() + considering[0].py() * dt_E1,
                           pos_ptn1.z() + considering[0].pz() * dt_E1, 0.);
              pos_ptn2.Set(pos_ptn2.x() + considering[1].px() * dt_E2,
                           pos_ptn2.y() + considering[1].py() * dt_E2,
                           pos_ptn2.z() + considering[1].pz() * dt_E2, 0.);
            }
            if ((dif2(pos_ptn3, pos_ptn1) > dist2cut) ||
                (dif2(pos_ptn3, pos_ptn2) > dist2cut) ||
                (dif2(pos_ptn1, pos_ptn2) > dist2cut)) {
              continue;
            }

            considering.add(HH_thermal[-element[2]]);
            HH_thermal[-element[2]].status(-993);
          } else {
            JSWARN << "SOMETHING WENT HORRIBLY WRONG - DO NOT KNOW WHERE "
                      "CURRENT QUARK CAME FROM?!";
          }

          // now that we *could* form a baryon, now we check if we actually do
          // form one baryon momentum
          FourVector Pbaryon;
          Pbaryon.Set(
              considering[0].px() + considering[1].px() + considering[2].px(),
              considering[0].py() + considering[1].py() + considering[2].py(),
              considering[0].pz() + considering[1].pz() + considering[2].pz(),
              0.);

          // baryon(CM) velocity
          FourVector betaB;  // really p[i]/e below
          betaB.Set(Pbaryon.x() / (considering[0].e() + considering[1].e() +
                                   considering[2].e()),
                    Pbaryon.y() / (considering[0].e() + considering[1].e() +
                                   considering[2].e()),
                    Pbaryon.z() / (considering[0].e() + considering[1].e() +
                                   considering[2].e()),
                    0.);
          betaB.Set(
              betaB.x(), betaB.y(), betaB.z(),
              1. / (sqrt(1. - (betaB.x() * betaB.x() + betaB.y() * betaB.y() +
                               betaB.z() * betaB.z()))));

          // boosting into CM frame
          FourVector pos_BCM[3], p_BCM[3];
          pos_BCM[0] = considering[0].boost_pos(betaB);
          pos_BCM[1] = considering[1].boost_pos(betaB);
          pos_BCM[2] = considering[2].boost_pos(betaB);
          p_BCM[0] = considering[0].boost_P(betaB);
          p_BCM[1] = considering[1].boost_P(betaB);
          p_BCM[2] = considering[2].boost_P(betaB);

          // velocities in CM frame
          FourVector v_BCM[3];
          v_BCM[0].Set(p_BCM[0].x() / p_BCM[0].t(), p_BCM[0].y() / p_BCM[0].t(),
                       p_BCM[0].z() / p_BCM[0].t(),
                       0.);  // these are really p[i]/e
          v_BCM[1].Set(p_BCM[1].x() / p_BCM[1].t(), p_BCM[1].y() / p_BCM[1].t(),
                       p_BCM[1].z() / p_BCM[1].t(), 0.);
          v_BCM[2].Set(p_BCM[2].x() / p_BCM[2].t(), p_BCM[2].y() / p_BCM[2].t(),
                       p_BCM[2].z() / p_BCM[2].t(), 0.);

          // propagating quarks until time of youngest quark
          double curtime = std::max(std::max(pos_BCM[0].t(), pos_BCM[1].t()),
                                    pos_BCM[2].t());
          FourVector cur_pos[3];
          cur_pos[0].Set(
              pos_BCM[0].x() + v_BCM[0].x() * (curtime - pos_BCM[0].t()),
              pos_BCM[0].y() + v_BCM[0].y() * (curtime - pos_BCM[0].t()),
              pos_BCM[0].z() + v_BCM[0].z() * (curtime - pos_BCM[0].t()),
              curtime);
          cur_pos[1].Set(
              pos_BCM[1].x() + v_BCM[1].x() * (curtime - pos_BCM[1].t()),
              pos_BCM[1].y() + v_BCM[1].y() * (curtime - pos_BCM[1].t()),
              pos_BCM[1].z() + v_BCM[1].z() * (curtime - pos_BCM[1].t()),
              curtime);
          cur_pos[2].Set(
              pos_BCM[2].x() + v_BCM[2].x() * (curtime - pos_BCM[2].t()),
              pos_BCM[2].y() + v_BCM[2].y() * (curtime - pos_BCM[2].t()),
              pos_BCM[2].z() + v_BCM[2].z() * (curtime - pos_BCM[2].t()),
              curtime);

          // finding position of CM at curtime
          FourVector pos_CM;
          pos_CM.Set((cur_pos[0].x() * considering[0].mass() +
                      cur_pos[1].x() * considering[1].mass() +
                      cur_pos[2].x() * considering[2].mass()) /
                         (considering[0].mass() + considering[1].mass() +
                          considering[2].mass()),
                     (cur_pos[0].y() * considering[0].mass() +
                      cur_pos[1].y() * considering[1].mass() +
                      cur_pos[2].y() * considering[2].mass()) /
                         (considering[0].mass() + considering[1].mass() +
                          considering[2].mass()),
                     (cur_pos[0].z() * considering[0].mass() +
                      cur_pos[1].z() * considering[1].mass() +
                      cur_pos[2].z() * considering[2].mass()) /
                         (considering[0].mass() + considering[1].mass() +
                          considering[2].mass()),
                     curtime);

          // finding position of baryon in lab frame
          betaB.Set(-betaB.x(), -betaB.y(), -betaB.z(), betaB.t());
          FourVector pos_lab = HHboost(betaB, pos_CM);

          // finding relative momenta of partons in CM frame
          FourVector k_rel_square[2];
          k_rel_square[0].Set(
              (considering[1].mass() * p_BCM[0].x() -
               considering[0].mass() * p_BCM[1].x()) /
                  (considering[0].mass() + considering[1].mass()),
              (considering[1].mass() * p_BCM[0].y() -
               considering[0].mass() * p_BCM[1].y()) /
                  (considering[0].mass() + considering[1].mass()),
              (considering[1].mass() * p_BCM[0].z() -
               considering[0].mass() * p_BCM[1].z()) /
                  (considering[0].mass() + considering[1].mass()),
              0.);
          k_rel_square[1].Set(
              (considering[2].mass() * (p_BCM[0].x() + p_BCM[1].x()) -
               (considering[0].mass() + considering[1].mass()) * p_BCM[2].x()) /
                  (considering[0].mass() + considering[1].mass() +
                   considering[2].mass()),
              (considering[2].mass() * (p_BCM[0].y() + p_BCM[1].y()) -
               (considering[0].mass() + considering[1].mass()) * p_BCM[2].y()) /
                  (considering[0].mass() + considering[1].mass() +
                   considering[2].mass()),
              (considering[2].mass() * (p_BCM[0].z() + p_BCM[1].z()) -
               (considering[0].mass() + considering[1].mass()) * p_BCM[2].z()) /
                  (considering[0].mass() + considering[1].mass() +
                   considering[2].mass()),
              0.);

          // finding relative positions of partons in CM frame
          FourVector pos_rel_square[2];
          pos_rel_square[0].Set((cur_pos[0].x() - cur_pos[1].x()),
                                (cur_pos[0].y() - cur_pos[1].y()),
                                (cur_pos[0].z() - cur_pos[1].z()), 0.);
          pos_rel_square[1].Set(
              ((cur_pos[0].x() * considering[0].mass() +
                cur_pos[1].x() * considering[1].mass()) /
                   (considering[0].mass() + considering[1].mass()) -
               cur_pos[2].x()),
              ((cur_pos[0].y() * considering[0].mass() +
                cur_pos[1].y() * considering[1].mass()) /
                   (considering[0].mass() + considering[1].mass()) -
               cur_pos[2].y()),
              ((cur_pos[0].z() * considering[0].mass() +
                cur_pos[1].z() * considering[1].mass()) /
                   (considering[0].mass() + considering[1].mass()) -
               cur_pos[2].z()),
              0.);

          double SigRB2 = SigNucR2;
          double SigLB2 = SigNucL2;
          int sortid[3] = {std::abs(considering[0].id()),
                           std::abs(considering[1].id()),
                           std::abs(considering[2].id())};
          std::stable_sort(std::begin(sortid), std::end(sortid),
                           std::greater<int>());

          // for particles we don't want to form, setting recofactor3 to 0
          if (sortid[0] == 3) {
            if (sortid[1] == 3) {
              if (sortid[2] == 3) {
                SigRB2 = SigOmgR2;
                SigLB2 = SigOmgL2;
              } else {
                SigRB2 = SigXiR2;
                SigLB2 = SigXiL2;
              }
            } else {
              SigRB2 = SigSigR2;
              SigLB2 = SigSigL2;
            }
          } else if (sortid[0] == 4) {
            if (sortid[1] == 4) {
              if (sortid[2] == 4) {
                SigRB2 = SigOcccR2;
                SigLB2 = SigOcccL2;
              } else if (sortid[2] == 3) {
                SigRB2 = SigOccR2;
                SigLB2 = SigOccL2;
              } else {
                SigRB2 = SigXiccR2;
                SigLB2 = SigXiccL2;
              }
            } else if (sortid[1] == 3) {
              if (sortid[2] == 3) {
                SigRB2 = SigOcR2;
                SigLB2 = SigOcL2;
              } else {
                SigRB2 = SigXicR2;
                SigLB2 = SigXicL2;
              }
            } else {
              SigRB2 = SigSigcR2;
              SigLB2 = SigSigcL2;
            }
          } else if (sortid[0] == 5) {
            if (sortid[1] == 5) {
              if (sortid[2] == 5) {
                SigRB2 = SigObbbR2;
                SigLB2 = SigObbbL2;
                recofactor3 = 0.;
              } else if (sortid[2] == 4) {
                SigRB2 = SigObbcR2;
                SigLB2 = SigObbcL2;
                recofactor3 = 0.;
              } else if (sortid[2] == 3) {
                SigRB2 = SigObbR2;
                SigLB2 = SigObbL2;
                recofactor3 = 0.;
              } else {
                SigRB2 = SigXibbR2;
                SigLB2 = SigXibbL2;
                recofactor3 = 0.;
              }
            } else if (sortid[1] == 4) {
              if (sortid[2] == 4) {
                SigRB2 = SigObccR2;
                SigLB2 = SigObccL2;
                recofactor3 = 0.;
              } else if (sortid[2] == 3) {
                SigRB2 = SigObcR2;
                SigLB2 = SigObcL2;
                recofactor3 = 0.;
              } else {
                SigRB2 = SigXibcR2;
                SigLB2 = SigXibcL2;
                recofactor3 = 0.;
              }
            } else if (sortid[1] == 3) {
              if (sortid[2] == 3) {
                SigRB2 = SigObR2;
                SigLB2 = SigObL2;
              } else {
                SigRB2 = SigXibR2;
                SigLB2 = SigXibL2;
              }
            } else {
              SigRB2 = SigSigbR2;
              SigLB2 = SigSigbL2;
            }
          }

          // precalc's for Wigner Wavefunction
          // 0:x, 1:y, 2:z ::: urho:(rel. between partons 1,2), ulamb:(rel.
          // between partons (1,2),3)
          double urho[3], ulamb[3];
          urho[0] =
              0.5 *
              (pos_rel_square[0].x() * pos_rel_square[0].x() / SigRB2 +
               k_rel_square[0].x() * k_rel_square[0].x() * SigRB2 / hbarc2);
          urho[1] =
              0.5 *
              (pos_rel_square[0].y() * pos_rel_square[0].y() / SigRB2 +
               k_rel_square[0].y() * k_rel_square[0].y() * SigRB2 / hbarc2);
          urho[2] =
              0.5 *
              (pos_rel_square[0].z() * pos_rel_square[0].z() / SigRB2 +
               k_rel_square[0].z() * k_rel_square[0].z() * SigRB2 / hbarc2);
          ulamb[0] =
              0.5 *
              (pos_rel_square[1].x() * pos_rel_square[1].x() / SigLB2 +
               k_rel_square[1].x() * k_rel_square[1].x() * SigLB2 / hbarc2);
          ulamb[1] =
              0.5 *
              (pos_rel_square[1].y() * pos_rel_square[1].y() / SigLB2 +
               k_rel_square[1].y() * k_rel_square[1].y() * SigLB2 / hbarc2);
          ulamb[2] =
              0.5 *
              (pos_rel_square[1].z() * pos_rel_square[1].z() / SigLB2 +
               k_rel_square[1].z() * k_rel_square[1].z() * SigLB2 / hbarc2);

          // 1D GS Wig. wavefunction
          double wig0[2][3];
          wig0[0][0] = std::exp(-urho[0]);
          wig0[0][1] = std::exp(-urho[1]);
          wig0[0][2] = std::exp(-urho[2]);
          wig0[1][0] = std::exp(-ulamb[0]);
          wig0[1][1] = std::exp(-ulamb[1]);
          wig0[1][2] = std::exp(-ulamb[2]);
          // 3D GS Wig. wavefunction
          double WigB[2];
          WigB[0] = wig0[0][0] * wig0[0][1] * wig0[0][2] * wig0[1][0] *
                    wig0[1][1] * wig0[1][2];

          // summing up 3D Wig. wavefunctions over nlev excited states
          WigB[1] = 0.;  // sumWigB;

          double wigE[2][3];
          for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
              wigE[i][j] = 0.;
            }
          }

          wigE[0][0] = wig0[0][0];
          for (int iRx = 0; iRx <= maxB_level; ++iRx) {
            wigE[0][1] = wig0[0][1];
            for (int iRy = 0; iRy <= maxB_level - iRx; ++iRy) {
              wigE[0][2] = wig0[0][2];
              for (int iRz = 0; iRz <= maxB_level - iRx - iRy; ++iRz) {
                wigE[1][0] = wig0[1][0];
                for (int iLx = 0; iLx <= maxB_level - iRx - iRy - iRz; ++iLx) {
                  wigE[1][1] = wig0[1][1];
                  for (int iLy = 0; iLy <= maxB_level - iRx - iRy - iRz - iLx;
                       ++iLy) {
                    wigE[1][2] = wig0[1][2];
                    for (int iLz = 0;
                         iLz <= maxB_level - iRx - iRy - iRz - iLx - iLy;
                         ++iLz) {
                      WigB[1] += wigE[0][0] * wigE[0][1] * wigE[0][2] *
                                 wigE[1][0] * wigE[1][1] * wigE[1][2];
                      wigE[1][2] *= ulamb[2] / ((double(iLz)) + 1.);
                    }
                    wigE[1][1] *= ulamb[1] / ((double(iLy)) + 1.);
                  }
                  wigE[1][0] *= ulamb[0] / ((double(iLx)) + 1.);
                }
                wigE[0][2] *= urho[2] / ((double(iRz)) + 1.);
              }
              wigE[0][1] *= urho[1] / ((double(iRy)) + 1.);
            }
            wigE[0][0] *= urho[0] / ((double(iRx)) + 1.);
          }

          if (maxB_level == -1) {
            WigB[0] = 0.;
            WigB[1] = 0.;
          }

          // Checking if baryon is formed (either ground or excited state)
          double rndbaryon = ran();
          double mult =
              (considering[1].is_thermal() || considering[2].is_thermal())
                  ? th_recofactor
                  : sh_recofactor;
          if (WigB[1] * recofactor3 * mult >= rndbaryon) {
            int junction_with_thermal_parton = 0;

            /*std::cout << "Baryons" << std::endl;
            std::cout << considering[0].id() << "," << considering[0].col() <<
            "," << considering[0].acol() << std::endl; std::cout <<
            considering[1].id() << "," << considering[1].col() << "," <<
            considering[1].acol() << std::endl;
            std::cout << considering[2].id() << "," << considering[2].col() <<
            "," << considering[2].acol() << std::endl;*/
            if (considering[0].id() * considering[1].id() *
                    considering[2].id() >
                0) {
              // baryon is to be formed, so temporary anti junction is defined
              if (considering[0].col() > 0 && considering[1].col() > 0 &&
                  considering[2].col() > 0) {
                IdColInfo1.push_back(-1);  // antijunction tag(-1)
                IdColInfo1.push_back(junction_with_thermal_parton);
                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    considering[0].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(considering[1].col());
                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(considering[2].col());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();

                // Since Baryon is formed, color tag for neutrality should be
                // added, casting int into double is needed!!(ex (double)
                // intvalue
                int coltag1 = considering[0].col();
                int coltag2 = considering[1].col();
                int coltag3 = considering[2].col();
                if (coltag1 > 0 && coltag2 > 0 && coltag3 > 0 &&
                    coltag1 <= limit && coltag2 <= limit && coltag3 <= limit) {
                  double tag1 = (double)
                      coltag1;  // they are casted to be inserted into the
                                // matrix(since it's vector of double
                  double tag2 = (double)coltag2;
                  double tag3 = (double)coltag3;
                  std::vector<int>::iterator I1 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                  std::vector<int>::iterator I2 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                  std::vector<int>::iterator I3 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag3);
                  int loc1 = std::distance(IndiceForColFin.begin(), I1);
                  int loc2 = std::distance(IndiceForColFin.begin(), I2);
                  int loc3 = std::distance(
                      IndiceForColFin.begin(),
                      I3);  // set up for find matrix indices corresponding to
                            // the color tags (we just found the corresponding
                            // indice in BaryonrecoMatrix1 with col tags )

                  BaryonrecoMatrix1.at(loc1).at(loc2).at(1) = tag3;
                  BaryonrecoMatrix1.at(loc2).at(loc1).at(1) = tag3;
                  BaryonrecoMatrix1.at(loc2).at(loc3).at(1) = tag1;
                  BaryonrecoMatrix1.at(loc3).at(loc2).at(1) = tag1;
                  BaryonrecoMatrix1.at(loc1).at(loc3).at(1) = tag2;
                  BaryonrecoMatrix1.at(loc3).at(loc1).at(1) =
                      tag2;  // now the color tag info for color neutrality is
                             // saved in Matrix, also we need to consider the
                             // impact from this to Meson Formation

                  // MesonrecoMatrix1 is modified by below
                  MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                  MesonrecoMatrix1.at(loc2).at(loc1) = 0;
                  MesonrecoMatrix1.at(loc2).at(loc3) = 0;
                  MesonrecoMatrix1.at(loc3).at(loc2) = 0;
                  MesonrecoMatrix1.at(loc1).at(loc3) = 0;
                  MesonrecoMatrix1.at(loc3).at(loc1) =
                      0;  // since three color tags of baryon are different from
                          // each other, so that these tags can't form meson
                          // with each other.
                }
                /*
                //dignostic measure
                std::cout <<endl<<"Meson reco Matrix revised by Baryon formation
                is same as below"<<endl; for(int irow=0; irow <
                IndiceForColFin.size(); irow++){ for(int icol=0; icol <
                IndiceForColFin.size(); icol++){ std::cout <<"
                "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                  }
                  std::cout <<endl<<endl;
                }
                //dignostic measure
                std::cout <<endl<<"Revised Baryon reco Matrix is same as
                below"<<endl; for(int irow=0; irow < IndiceForColFin.size();
                irow++){ for(int icol=0; icol < IndiceForColFin.size(); icol++){
                    std::cout <<"  (
                "<<BaryonrecoMatrix1.at(irow).at(icol).at(1)<<" )  ";
                  }
                  std::cout <<endl<<endl;
                }
                */
              }
              if (considering[0].col() > 0 && considering[1].col() == 0 &&
                  considering[2].col() > 0) {  // MAT + Therm/LBT + MAT case
                IdColInfo1.push_back(-1);      // antijunction tag(-1)

                maxtag++;
                int loc = findcloserepl(considering[1], perm2[q2], true, true,
                                        showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].acol(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                // TODO: give thermal partons "ANTI COLOR TAGS" to form anti
                // junction and conserve baryon number.
                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    considering[0].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(maxtag);
                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(considering[2].col());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].col() == 0 && considering[1].col() > 0 &&
                  considering[2].col() > 0) {  // LBT + MAT + MAT case
                IdColInfo1.push_back(-1);      // antijunction tag(-1)

                maxtag++;
                int loc = findcloserepl(considering[0], element[0] + 1, true,
                                        true, showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].acol(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    maxtag);  // {-1, anticolor tag} will be at 2nd, 3rd, 4th
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(considering[1].col());
                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(considering[2].col());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].col() > 0 && considering[1].col() > 0 &&
                  considering[2].col() == 0) {  // MAT + MAT + Therm/LBT case
                IdColInfo1.push_back(-1);       // antijunction tag(-1)

                maxtag++;
                int loc = findcloserepl(considering[2], perm2[q3], true, true,
                                        showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].acol(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    considering[0].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(considering[1].col());
                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(maxtag);

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with one zero
                 // component.
              if (considering[0].col() > 0 && considering[1].col() == 0 &&
                  considering[2].col() ==
                      0) {                 // MAT + Therm/LBT + Therm/LBT case
                IdColInfo1.push_back(-1);  // antijunction tag(-1)
                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    considering[0].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th

                maxtag++;
                int loc1 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].acol(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].acol(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(maxtag);
                // TODO: give thermal partons "ANTI COLOR TAGS" to form anti
                // junction and conserve baryon number.

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].col() == 0 && considering[1].col() == 0 &&
                  considering[2].col() >
                      0) {                 // Therm/LBT + Therm/LBT + MAT case
                IdColInfo1.push_back(-1);  // antijunction tag(-1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].acol(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].acol(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(maxtag);
                // TODO: give thermal partons "ANTI COLOR TAGS" to form anti
                // junction and conserve baryon number.

                IdColInfo2.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo2.push_back(
                    considering[2].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].col() == 0 && considering[1].col() > 0 &&
                  considering[2].col() ==
                      0) {                 // Therm/LBT + Therm/LBT + MAT case
                IdColInfo1.push_back(-1);  // antijunction tag(-1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].acol(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo2.push_back(-1);
                IdColInfo2.push_back(maxtag);

                IdColInfo3.push_back(
                    -1);  // means anticolor tag(negative color charge)
                IdColInfo3.push_back(
                    considering[1].col());  // {-1, anticolor tag} will be at
                                            // 2nd, 3rd, 4th

                maxtag++;
                int loc2 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].acol(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(maxtag);

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with two zeros
              if (considering[0].col() == 0 && considering[1].col() == 0 &&
                  considering[2].col() ==
                      0) {                 // Therm/LBT + Therm/LBT + Therm/LBT
                IdColInfo1.push_back(-1);  // antijunction tag(-1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].acol(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo2.push_back(-1);
                IdColInfo2.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].acol(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo3.push_back(-1);
                IdColInfo3.push_back(maxtag);

                maxtag++;
                int loc3 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc3 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.acol(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc3 > 0) {
                  showerquarks[loc3 - 1].acol(maxtag);
                } else if (loc3 < 0) {
                  HH_thermal[-loc3 - 1].acol(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(-1);
                IdColInfo4.push_back(maxtag);

                // TODO: give thermal partons "ANTI COLOR TAGS" to form anti
                // junction and conserve baryon number.
                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);
                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry
                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with three zeros
            } else if (considering[0].id() * considering[1].id() *
                           considering[2].id() <
                       0) {  // anti baryon would be made so temporary junction
                             // is defined
              // std::cout <<endl<<"Anti Baryon is made from "<<" (
              // "<<considering[0].acol()<<" , "<<considering[1].acol()<<" ,
              // "<<considering[2].acol()<<" ) "<<endl; std::cout <<endl<<"Anti
              // Baryon is made from "<<" (
              // "<<HH_showerptns[showerquarks[element[0]].par()].acol()<<" ,
              // "<<HH_showerptns[showerquarks[element[1]].par()].acol()<<" ,
              // "<<HH_showerptns[showerquarks[element[2]].par()].acol()<<" )
              // "<<endl;

              if (considering[0].acol() > 0 && considering[1].acol() > 0 &&
                  considering[2].acol() > 0) {
                IdColInfo1.push_back(1);  // junction tag(1)
                IdColInfo1.push_back(
                    junction_with_thermal_parton);  // zero(just room for the
                                                    // other usage)   : {1, 0}
                                                    // at 1st
                IdColInfo2.push_back(
                    1);  // means color tag(positive color charge)   : { 1,
                         // color tag } at 2nd, 3rd, 4th in the vector
                IdColInfo2.push_back(considering[0].acol());
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(considering[1].acol());
                IdColInfo4.push_back(1);
                IdColInfo4.push_back(considering[2].acol());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();

                // Since Baryon is formed, color tag for neutrality should be
                // added, casting int into double is needed!!(ex (double)
                // intvalue
                int coltag1 = considering[0].acol();
                int coltag2 = considering[1].acol();
                int coltag3 = considering[2].acol();
                if (coltag1 > 0 && coltag2 > 0 && coltag3 > 0 &&
                    coltag1 <= limit && coltag2 <= limit && coltag3 <= limit) {
                  double tag1 = (double)
                      coltag1;  // they are casted to be inserted into the
                                // matrix(since it's vector of double
                  double tag2 = (double)coltag2;
                  double tag3 = (double)coltag3;
                  std::vector<int>::iterator I1 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                  std::vector<int>::iterator I2 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                  std::vector<int>::iterator I3 = std::find(
                      IndiceForColFin.begin(), IndiceForColFin.end(), coltag3);
                  int loc1 = std::distance(IndiceForColFin.begin(), I1);
                  int loc2 = std::distance(IndiceForColFin.begin(), I2);
                  int loc3 = std::distance(
                      IndiceForColFin.begin(),
                      I3);  // set up for find matrix indices corresponding to
                            // the color tags (we just found the corresponding
                            // indice in BaryonrecoMatrix1 with col tags )

                  BaryonrecoMatrix1.at(loc1).at(loc2).at(1) = tag3;
                  BaryonrecoMatrix1.at(loc2).at(loc1).at(1) = tag3;
                  BaryonrecoMatrix1.at(loc2).at(loc3).at(1) = tag1;
                  BaryonrecoMatrix1.at(loc3).at(loc2).at(1) = tag1;
                  BaryonrecoMatrix1.at(loc1).at(loc3).at(1) = tag2;
                  BaryonrecoMatrix1.at(loc3).at(loc1).at(1) =
                      tag2;  // now the color tag info for color neutrality is
                             // saved in Matrix, also we need to consider the
                             // impact from this to Meson Formation

                  // MesonrecoMatrix1 is modified by below
                  MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                  MesonrecoMatrix1.at(loc2).at(loc1) = 0;
                  MesonrecoMatrix1.at(loc2).at(loc3) = 0;
                  MesonrecoMatrix1.at(loc3).at(loc2) = 0;
                  MesonrecoMatrix1.at(loc1).at(loc3) = 0;
                  MesonrecoMatrix1.at(loc3).at(loc1) =
                      0;  // since three color tags of baryon are different from
                          // each other, so that these tags can't form meson
                          // with each other.
                }
              }
              if (considering[0].acol() > 0 && considering[1].acol() == 0 &&
                  considering[2].acol() > 0) {  // MAT + Therm/LBT + MAT case
                IdColInfo1.push_back(1);        // junction tag(1)

                maxtag++;
                int loc = findcloserepl(considering[1], perm2[q2], true, true,
                                        showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].col(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo2.push_back(1);  // color tag
                IdColInfo2.push_back(
                    considering[0]
                        .acol());  // {1, color tag} will be at 2nd, 3rd, 4th
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(maxtag);
                IdColInfo4.push_back(1);
                IdColInfo4.push_back(considering[2].acol());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].acol() == 0 && considering[1].acol() > 0 &&
                  considering[2].acol() > 0) {  // LBT + MAT + MAT case
                IdColInfo1.push_back(1);        // junction tag(1)

                maxtag++;
                int loc = findcloserepl(considering[0], element[0] + 1, true,
                                        true, showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].col(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo2.push_back(1);  // means color tag
                IdColInfo2.push_back(
                    maxtag);  // {1, color tag} will be at 2nd, 3rd, 4th
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(considering[1].acol());
                IdColInfo4.push_back(1);
                IdColInfo4.push_back(considering[2].acol());

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].acol() > 0 && considering[1].acol() > 0 &&
                  considering[2].acol() == 0) {  // MAT + MAT + Therm/LBT case
                IdColInfo1.push_back(1);         // junction tag(1)

                maxtag++;
                int loc = findcloserepl(considering[2], perm2[q3], true, true,
                                        showerquarks, HH_thermal);
                if (loc == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc > 0) {
                  showerquarks[loc - 1].col(maxtag);
                } else if (loc < 0) {
                  HH_thermal[-loc - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo2.push_back(1);  // means color tag
                IdColInfo2.push_back(
                    considering[0]
                        .acol());  // {1, color tag} will be at 2nd, 3rd, 4th
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(considering[1].acol());
                IdColInfo4.push_back(1);
                IdColInfo4.push_back(maxtag);

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with one zero
                 // component.
              if (considering[0].acol() > 0 && considering[1].acol() == 0 &&
                  considering[2].acol() ==
                      0) {                // MAT + Therm/LBT + Therm/LBT case
                IdColInfo1.push_back(1);  // junction tag(1)
                IdColInfo2.push_back(1);  // means color tag
                IdColInfo2.push_back(
                    considering[0]
                        .acol());  // {1, color tag} will be at 2nd, 3rd, 4th

                maxtag++;
                int loc1 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].col(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].col(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(1);
                IdColInfo4.push_back(maxtag);
                // TODO: give thermal partons "ANTI COLOR TAGS" to form anti
                // junction and conserve baryon number.

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].acol() == 0 && considering[1].acol() == 0 &&
                  considering[2].acol() >
                      0) {                // Therm/LBT + Therm/LBT + MAT case
                IdColInfo1.push_back(1);  // junction tag(1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].col(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo2.push_back(1);
                IdColInfo2.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].col(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo3.push_back(1);
                IdColInfo3.push_back(maxtag);

                IdColInfo4.push_back(1);  // means color tag
                IdColInfo4.push_back(
                    considering[2]
                        .acol());  // {1, color tag} will be at 2nd, 3rd, 4th

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }
              if (considering[0].acol() == 0 && considering[1].acol() > 0 &&
                  considering[2].acol() ==
                      0) {                // Therm/LBT + Therm/LBT + MAT case
                IdColInfo1.push_back(1);  // junction tag(1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].col(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo2.push_back(1);
                IdColInfo2.push_back(maxtag);

                IdColInfo3.push_back(1);  // means color tag
                IdColInfo3.push_back(
                    considering[1]
                        .acol());  // {1, color tag} will be at 2nd, 3rd, 4th

                maxtag++;
                int loc2 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].col(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(1);
                IdColInfo4.push_back(maxtag);

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);

                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry

                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with two zeros
              if (considering[0].acol() == 0 && considering[1].acol() == 0 &&
                  considering[2].acol() ==
                      0) {                // Therm/LBT + Therm/LBT + Therm/LBT
                IdColInfo1.push_back(1);  // junction tag(1)

                maxtag++;
                int loc1 = findcloserepl(considering[0], element[0] + 1, true,
                                         true, showerquarks, HH_thermal);
                if (loc1 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[0];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc1 > 0) {
                  showerquarks[loc1 - 1].col(maxtag);
                } else if (loc1 < 0) {
                  HH_thermal[-loc1 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo2.push_back(1);
                IdColInfo2.push_back(maxtag);

                maxtag++;
                int loc2 = findcloserepl(considering[1], perm2[q2], true, true,
                                         showerquarks, HH_thermal);
                if (loc2 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[1];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc2 > 0) {
                  showerquarks[loc2 - 1].col(maxtag);
                } else if (loc2 < 0) {
                  HH_thermal[-loc2 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }
                IdColInfo3.push_back(1);
                IdColInfo3.push_back(maxtag);

                maxtag++;
                int loc3 = findcloserepl(considering[2], perm2[q3], true, true,
                                         showerquarks, HH_thermal);
                if (loc3 == 999999999) {
                  // std::cout <<endl<<"Warning : extra parton used for string
                  // repair!!"<<endl;
                  HHparton fakep = considering[2];
                  fakep.id(-fakep.id());
                  fakep.col(maxtag);
                  fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                  fakep.px(0.);
                  fakep.py(0.);
                  fakep.pz(0.);
                  fakep.e(fakep.mass());
                  fakep.is_fakeparton(true);
                  Extraparton.add(fakep);
                  junction_with_thermal_parton = 1;
                  // somewhere, we need to make for loop to toss all partons to
                  // remnants list.
                } else if (loc3 > 0) {
                  showerquarks[loc3 - 1].col(maxtag);
                } else if (loc3 < 0) {
                  HH_thermal[-loc3 - 1].col(maxtag);
                  junction_with_thermal_parton = 1;
                }

                IdColInfo1.push_back(junction_with_thermal_parton);

                IdColInfo4.push_back(1);
                IdColInfo4.push_back(maxtag);

                JunctionInfo.push_back(IdColInfo1);
                JunctionInfo.push_back(IdColInfo2);
                JunctionInfo.push_back(IdColInfo3);
                JunctionInfo.push_back(IdColInfo4);
                Tempjunctions.push_back(
                    JunctionInfo);  // information of tempjunction is saved so
                                    // et clear subordinate vector for next
                                    // entry
                IdColInfo1.clear();
                IdColInfo2.clear();
                IdColInfo3.clear();
                IdColInfo4.clear();
                JunctionInfo.clear();
              }  // so far, we are finished with the list with three zeros
                 // dignostic measure
              /*
              std::cout <<endl<<"Meson reco Matrix revised by Baryon formation
              is same as below"<<endl; for(int irow=0; irow <
              IndiceForColFin.size(); irow++){ for(int icol=0; icol <
              IndiceForColFin.size(); icol++){ std::cout <<"
              "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                  }
              std::cout <<endl<<endl;
              }
              */
            }

            // now we're forming the hadron
            HHhadron formedhadron;
            // setting the hadron values: is a recombined hadron, mass, and
            // parents
            formedhadron.is_recohad(true);
            formedhadron.mass(p_BCM[0].t() + p_BCM[1].t() + p_BCM[2].t());
            formedhadron.add_par(showerquarks[element[0]].par());
            if (perm2[q2] > 0) {
              formedhadron.add_par(showerquarks[element[1]].par());
            } else {
              formedhadron.add_par(element[1]);
            }
            if (perm2[q3] > 0) {
              formedhadron.add_par(showerquarks[element[2]].par());
            } else {
              formedhadron.add_par(element[2]);
            }

            // now setting if baryon is in excited state
            if (WigB[0] * recofactor3 * mult < rndbaryon) {
              formedhadron.is_excited(true);
            }

            // setting if there are any thermal partons used to make the hadron
            // (with the '0' thermal parent as -99999 so that it doesn't
            // conflict with '0' shower parton)
            if (perm2[q2] > 0 && perm2[q3] > 0) { /*is sh-sh-sh*/
              formedhadron.is_shsh(true);
            } else if (perm2[q2] > 0 && perm2[q3] < 0) { /*is sh-sh-th*/
              formedhadron.is_shth(true);
              if (element[2] == 0) {
                formedhadron.parents[2] = -99999;
              }
            } else if (perm2[q2] < 0 && perm2[q3] > 0) { /*is sh-th-sh*/
              formedhadron.is_shth(true);
              if (element[1] == 0) {
                formedhadron.parents[1] = -99999;
              }
            } else if (perm2[q2] < 0 && perm2[q3] < 0) { /*is sh-th-th*/
              formedhadron.is_shth(true);
              if (element[1] == 0) {
                formedhadron.parents[1] = -99999;
              }
              if (element[2] == 0) {
                formedhadron.parents[2] = -99999;
              }
            }

            // setting hadron position and momentum vectors
            Pbaryon.Set(
                Pbaryon.x(), Pbaryon.y(), Pbaryon.z(),
                sqrt(Pbaryon.x() * Pbaryon.x() + Pbaryon.y() * Pbaryon.y() +
                     Pbaryon.z() * Pbaryon.z() +
                     formedhadron.mass() * formedhadron.mass()));
            formedhadron.pos(pos_lab);
            formedhadron.P(Pbaryon);

            // setting hadron color tags (for tracing colors of the constituent
            // partons) [**] will need to update to reflect color tags given to
            // random (thermal/lbt) partons
            formedhadron.add_col((considering[0].col() > 0)
                                     ? considering[0].col()
                                     : considering[0].acol());
            formedhadron.add_col((considering[1].col() > 0)
                                     ? considering[1].col()
                                     : considering[1].acol());
            formedhadron.add_col((considering[2].col() > 0)
                                     ? considering[2].col()
                                     : considering[2].acol());

            // need to choose *what* hadron we've formed... base this on the
            // parton id's, mass, & if excited might want to do this
            // differently? void f'n(partoncollection, formedhadron)?
            set_baryon_id(considering, formedhadron);

            // need to add the hadron to the collection
            HH_hadrons.add(formedhadron);

            // now that we've formed the hadron, need to set ALL (3) the
            // 'considering' flags to used
            showerquarks[element[0]].status(1);
            showerquarks[element[0]].is_used(true);
            if (perm2[q2] > 0) {
              showerquarks[element[1]].status(1);
              showerquarks[element[1]].is_used(true);
            } else {
              HH_thermal[-element[1]].status(1);
              HH_thermal[-element[1]].is_used(true);
            }
            if (perm2[q3] > 0) {
              showerquarks[element[2]].status(1);
              showerquarks[element[2]].is_used(true);
            } else {
              HH_thermal[-element[2]].status(1);
              HH_thermal[-element[2]].is_used(true);
            }

            madehadron = true;
            considering.clear();
            break;
          } else {
            // since we've not formed a baryon on this try, need to revert the
            // third quark 'considering' used flag
            if (perm2[q3] > 0) {
              showerquarks[element[2]].status(0);
            } else {
              HH_thermal[-element[2]].status(0);
            }

            // and remove the third entry in considering
            considering.partons.pop_back();
          }
        }
      } else if (considering[0].id() * considering[1].id() < 0) {
        // Key point is determing recofactor2
        if (considering[0].id() > 0 && considering[1].id() < 0) {
          int tag0 = considering[0].col();
          int tag1 = considering[1].acol();
          if (tag0 > 0 && tag1 > 0 && tag0 <= limit && tag1 <= limit) {
            std::vector<int>::iterator L1 =
                std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag0);
            std::vector<int>::iterator L2 =
                std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
            int indexMatrix1 = std::distance(IndiceForColFin.begin(), L1);
            int indexMatrix2 = std::distance(IndiceForColFin.begin(), L2);

            recofactor2 = (MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2));
            // std::cout <<endl<<" the recofactor for meson is
            // "<<(MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2))<<endl;
          } else {
            recofactor2 = 1. / 9.;
          }
        } else if (considering[1].id() > 0 && considering[0].id() < 0) {
          int tag0 = considering[1].col();
          int tag1 = considering[0].acol();
          if (tag0 > 0 && tag1 > 0 && tag0 <= limit && tag1 <= limit) {
            std::vector<int>::iterator L1 =
                std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag0);
            std::vector<int>::iterator L2 =
                std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
            int indexMatrix1 = std::distance(IndiceForColFin.begin(), L1);
            int indexMatrix2 = std::distance(IndiceForColFin.begin(), L2);

            recofactor2 = (MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2));
          } else {
            recofactor2 = 1. / 9.;
          }
          // std::cout <<endl<<" the recofactor for meson is
          // "<<(MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2))<<endl;
        }
        // now that we *could* form a meson, now we check if we actually do form
        // one meson momentum
        FourVector Pmeson;
        Pmeson.Set(considering[0].px() + considering[1].px(),
                   considering[0].py() + considering[1].py(),
                   considering[0].pz() + considering[1].pz(), 0.);

        // meson(CM) velocity
        FourVector betaM;
        betaM.Set(Pmeson.x() / (considering[0].e() + considering[1].e()),
                  Pmeson.y() / (considering[0].e() + considering[1].e()),
                  Pmeson.z() / (considering[0].e() + considering[1].e()), 0.);
        betaM.Set(
            betaM.x(), betaM.y(), betaM.z(),
            1. / (sqrt(1. - (betaM.x() * betaM.x() + betaM.y() * betaM.y() +
                             betaM.z() * betaM.z()))));

        // boosting into CM frame
        FourVector pos_MCM[2], p_MCM[2];
        pos_MCM[0] = considering[0].boost_pos(betaM);
        pos_MCM[1] = considering[1].boost_pos(betaM);
        p_MCM[0] = considering[0].boost_P(betaM);
        p_MCM[1] = considering[1].boost_P(betaM);

        // velocities in CM frame
        FourVector v_MCM[2];
        v_MCM[0].Set(p_MCM[0].x() / p_MCM[0].t(), p_MCM[0].y() / p_MCM[0].t(),
                     p_MCM[0].z() / p_MCM[0].t(), 0.);
        v_MCM[1].Set(p_MCM[1].x() / p_MCM[1].t(), p_MCM[1].y() / p_MCM[1].t(),
                     p_MCM[1].z() / p_MCM[1].t(), 0.);

        // propagating quarks until time of youngest quark
        // is just max(pos_MCM[0].t(), pos_MCM[1].t());
        double curtime =
            (pos_MCM[0].t() > pos_MCM[1].t()) ? pos_MCM[0].t() : pos_MCM[1].t();
        FourVector cur_pos[2];
        cur_pos[0].Set(
            pos_MCM[0].x() + v_MCM[0].x() * (curtime - pos_MCM[0].t()),
            pos_MCM[0].y() + v_MCM[0].y() * (curtime - pos_MCM[0].t()),
            pos_MCM[0].z() + v_MCM[0].z() * (curtime - pos_MCM[0].t()),
            curtime);
        cur_pos[1].Set(
            pos_MCM[1].x() + v_MCM[1].x() * (curtime - pos_MCM[1].t()),
            pos_MCM[1].y() + v_MCM[1].y() * (curtime - pos_MCM[1].t()),
            pos_MCM[1].z() + v_MCM[1].z() * (curtime - pos_MCM[1].t()),
            curtime);

        // finding position of CM at curtime
        FourVector pos_CM;
        pos_CM.Set((cur_pos[0].x() * considering[0].mass() +
                    cur_pos[1].x() * considering[1].mass()) /
                       (considering[0].mass() + considering[1].mass()),
                   (cur_pos[0].y() * considering[0].mass() +
                    cur_pos[1].y() * considering[1].mass()) /
                       (considering[0].mass() + considering[1].mass()),
                   (cur_pos[0].z() * considering[0].mass() +
                    cur_pos[1].z() * considering[1].mass()) /
                       (considering[0].mass() + considering[1].mass()),
                   curtime);

        // finding position of meson in lab frame
        betaM.Set(-betaM.x(), -betaM.y(), -betaM.z(), betaM.t());
        FourVector pos_lab = HHboost(betaM, pos_CM);

        // finding the squares of the relative momenta of partons in CM frame
        FourVector k_rel_square;
        double sum_mass_square =
            (considering[0].mass() + considering[1].mass()) *
            (considering[0].mass() + considering[1].mass());
        k_rel_square.Set(std::pow(considering[1].mass() * p_MCM[0].x() -
                                      considering[0].mass() * p_MCM[1].x(),
                                  2.) /
                             sum_mass_square,
                         std::pow(considering[1].mass() * p_MCM[0].y() -
                                      considering[0].mass() * p_MCM[1].y(),
                                  2.) /
                             sum_mass_square,
                         std::pow(considering[1].mass() * p_MCM[0].z() -
                                      considering[0].mass() * p_MCM[1].z(),
                                  2.) /
                             sum_mass_square,
                         0.);
        k_rel_square.Set(
            k_rel_square.x(), k_rel_square.y(), k_rel_square.z(),
            k_rel_square.x() + k_rel_square.y() + k_rel_square.z());

        // finding the squares of relative positions of partons in CM frame
        FourVector pos_rel_square;
        pos_rel_square.Set((cur_pos[0].x() - cur_pos[1].x()) *
                               (cur_pos[0].x() - cur_pos[1].x()),
                           (cur_pos[0].y() - cur_pos[1].y()) *
                               (cur_pos[0].y() - cur_pos[1].y()),
                           (cur_pos[0].z() - cur_pos[1].z()) *
                               (cur_pos[0].z() - cur_pos[1].z()),
                           0.);
        pos_rel_square.Set(
            pos_rel_square.x(), pos_rel_square.y(), pos_rel_square.z(),
            pos_rel_square.x() + pos_rel_square.y() + pos_rel_square.z());

        // setting appropriate sigma...
        double SigM2 = SigPi2;
        int sortid[2] = {0, 0};
        if (std::abs(considering[0].id()) >= std::abs(considering[1].id())) {
          sortid[0] = std::abs(considering[0].id());
          sortid[1] = std::abs(considering[1].id());
        } else {
          sortid[0] = std::abs(considering[1].id());
          sortid[1] = std::abs(considering[0].id());
        }

        if (sortid[0] == 3) {
          if (sortid[1] == 3) {
            SigM2 = SigPhi2;
          } else {
            SigM2 = SigK2;
          }
        } else if (sortid[0] == 4) {
          if (sortid[1] == 4) {
            SigM2 = SigJpi2;
          } else if (sortid[1] == 3) {
            SigM2 = SigDs2;
          } else {
            SigM2 = SigD2;
          }
        } else if (sortid[0] == 5) {
          if (sortid[1] == 5) {
            SigM2 = SigUps2;
          } else if (sortid[1] == 4) {
            SigM2 = SigBc2;
          } else if (sortid[1] == 3) {
            SigM2 = SigB2;
          } else {
            SigM2 = SigB2;
          }
        }

        double u[4];
        // This is the squared distance in phase space weighted with the widths
        u[1] = 0.5 *
               (pos_rel_square.x() / SigM2 + k_rel_square.x() * SigM2 / hbarc2);
        u[2] = 0.5 *
               (pos_rel_square.y() / SigM2 + k_rel_square.y() * SigM2 / hbarc2);
        u[3] = 0.5 *
               (pos_rel_square.z() / SigM2 + k_rel_square.z() * SigM2 / hbarc2);
        u[0] = u[1] + u[2] + u[3];

        // Ground state wave function
        double WigM = std::exp(-u[0]);

        // Computing s ~ L^2 for the system
        double rdotr = pos_rel_square.t();
        double pdotp = k_rel_square.t();
        double pdotr =
            std::sqrt(pos_rel_square.x()) * std::sqrt(k_rel_square.x()) +
            std::sqrt(pos_rel_square.y()) * std::sqrt(k_rel_square.y()) +
            std::sqrt(pos_rel_square.z()) * std::sqrt(k_rel_square.z());

        double s = 1 / hbarc2 * (pdotp * rdotr - pdotr * pdotr);

        // Random number for recombination dice roll
        double rndmeson = ran();

        // Initialize quantum numbers and recombination probability
        int angular_qnum = -1;  // l
        int radial_qnum = -1;   // k
        double mult1 =
            considering[1].is_thermal() ? th_recofactor : sh_recofactor;
        double total_prob = WigM * recofactor2 * mult1;

        if (total_prob >= rndmeson) {
          angular_qnum = 0;
          radial_qnum = 0;
        } else {
          total_prob += WigM * u[0];
          if (total_prob >= rndmeson && maxM_level > 0) {
            angular_qnum = 1;
            radial_qnum = 0;
          } else {
            total_prob += (1. / 2.) * WigM *
                          ((2. / 3.) * std::pow(u[0], 2) + (1. / 3.) * s);
            if (total_prob >= rndmeson && maxM_level > 1) {
              angular_qnum = 2;
              radial_qnum = 0;
            } else {
              total_prob += (1. / 2.) * WigM *
                            ((1. / 3.) * std::pow(u[0], 2) - (1. / 3.) * s);
              if (total_prob >= rndmeson && maxM_level > 1) {
                angular_qnum = 0;
                radial_qnum = 1;
              } else {
                total_prob +=
                    (1. / 6.) * WigM *
                    ((2. / 5.) * std::pow(u[0], 3) + (3. / 5.) * u[0] * s);
                if (total_prob >= rndmeson && maxM_level > 2) {
                  angular_qnum = 1;
                  radial_qnum = 1;
                } else {
                  total_prob +=
                      (1. / 6.) * WigM *
                      ((3. / 5.) * std::pow(u[0], 3) - (3. / 5.) * u[0] * s);
                  if (total_prob >= rndmeson && maxM_level > 2) {
                    angular_qnum = 3;
                    radial_qnum = 0;
                  } else {
                    total_prob += (1. / 120.) * WigM *
                                  std::pow((std::pow(u[0], 2) - s), 2);
                    if (total_prob >= rndmeson && maxM_level > 3) {
                      angular_qnum = 0;
                      radial_qnum = 2;
                    } else {
                      total_prob += (1. / 24.) * WigM *
                                    ((4. / 7.) * std::pow(u[0], 4) -
                                     (2. / 7.) * std::pow(u[0], 2) * s -
                                     (2. / 7.) * std::pow(s, 2));
                      if (total_prob >= rndmeson && maxM_level > 3) {
                        angular_qnum = 2;
                        radial_qnum = 1;
                      } else {
                        total_prob += (1. / 24.) * WigM *
                                      ((8. / 35.) * std::pow(u[0], 4) +
                                       (24. / 35.) * std::pow(u[0], 2) * s +
                                       (3. / 35.) * std::pow(s, 2));
                        if (total_prob >= rndmeson && maxM_level > 3) {
                          angular_qnum = 4;
                          radial_qnum = 0;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }

        // Habemus Recombinationem!
        if (angular_qnum >= 0) {
          /*std::cout << "Mesons" << std::endl;
          std::cout << considering[0].id() << "," << considering[0].col() << ","
          << considering[0].acol() << std::endl; std::cout <<
          considering[1].id() << "," << considering[1].col() << "," <<
          considering[1].acol() << std::endl;*/
          if (considering[0].id() > 0 &&
              considering[1].id() <
                  0) {  // case of first parton is q and second is q-bar
            // std::cout <<endl<<"chosen partons are "<<
            // int(considering[0].id()) << " and " << int(considering[1].id())
            // << endl <<endl; std::cout <<endl<<"and their color tag is "<<
            // int(considering[0].col()) << " and " <<
            // int(considering[1].acol()) << endl <<endl; std::cout <<"color
            // correction implemented as Follows: "<< considering[0].col() <<" =
            // " << considering[1].acol()<<endl;

            // Overwrite non-dominant colortags in tempjunctions with dominant
            // ones q-qbar case
            if (considering[0].col() != 0 && considering[1].acol() != 0) {
              for (int ijunc = 0; ijunc < Tempjunctions.size(); ijunc++) {
                if (Tempjunctions.at(ijunc).at(1).at(1) ==
                    considering[1].acol()) {
                  Tempjunctions.at(ijunc).at(1).pop_back();
                  Tempjunctions.at(ijunc).at(1).push_back(considering[0].col());
                }
                if (Tempjunctions.at(ijunc).at(2).at(1) ==
                    considering[1].acol()) {
                  Tempjunctions.at(ijunc).at(2).pop_back();
                  Tempjunctions.at(ijunc).at(2).push_back(considering[0].col());
                }
                if (Tempjunctions.at(ijunc).at(3).at(1) ==
                    considering[1].acol()) {
                  Tempjunctions.at(ijunc).at(3).pop_back();
                  Tempjunctions.at(ijunc).at(3).push_back(considering[0].col());
                }
              }
            }

            // before changing the tags, revise the MesonrecoMatrix elements!
            int coltag1 = considering[0].col();
            int coltag2 = considering[1].acol();
            if (coltag1 > 0 && coltag2 > 0 && coltag1 <= limit &&
                coltag2 <= limit) {
              std::vector<int>::iterator I1 = std::find(
                  IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
              std::vector<int>::iterator I2 = std::find(
                  IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
              int loc1 = std::distance(IndiceForColFin.begin(), I1);
              int loc2 = std::distance(IndiceForColFin.begin(),
                                       I2);  // set up for find matrix indices
                                             // corresponding to the color tags

              MesonrecoMatrix1.at(loc1).at(loc2) = 1;
              MesonrecoMatrix1.at(loc2).at(loc1) = 1;  // Matrix revised

              for (int iM = 0; iM < MesonrecoMatrix1[0].size(); ++iM) {
                if (MesonrecoMatrix1[loc1][iM] == 1) {
                  MesonrecoMatrix1[loc2][iM] = 1;
                }
                if (MesonrecoMatrix1[loc2][iM] == 1) {
                  MesonrecoMatrix1[loc1][iM] = 1;
                }
                if (MesonrecoMatrix1[iM][loc1] == 1) {
                  MesonrecoMatrix1[iM][loc2] = 1;
                }
                if (MesonrecoMatrix1[iM][loc2] == 1) {
                  MesonrecoMatrix1[iM][loc1] = 1;
                }
                if (MesonrecoMatrix1[loc1][iM] == 0) {
                  MesonrecoMatrix1[loc2][iM] = 0;
                }
                if (MesonrecoMatrix1[loc2][iM] == 0) {
                  MesonrecoMatrix1[loc1][iM] = 0;
                }
                if (MesonrecoMatrix1[iM][loc1] == 0) {
                  MesonrecoMatrix1[iM][loc2] = 0;
                }
                if (MesonrecoMatrix1[iM][loc2] == 0) {
                  MesonrecoMatrix1[iM][loc1] = 0;
                }
              }
              /*
              std::cout <<endl<<"Revised Matrix is same as below"<<endl;
              for(int irow=0; irow < IndiceForColFin.size(); irow++){
                for(int icol=0; icol < IndiceForColFin.size(); icol++){
                  std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                }
                std::cout <<endl<<endl;
              }
              */
            }

            // possible case  MAT+LBT, MAT+THERM,
            // treatment 1 : based on distance.
            if (considering[0].col() > 0 &&
                considering[1].acol() >
                    0) {  // MAT/lbt or therm with color tags + MAT/lbt or therm
                          // with color tags
              if (perm2[q2] > 0) {
                HH_showerptns[showerquarks[element[1]].par()].acol(
                    considering[0]
                        .col());  // now color tags from both partons are same
                // also need to set the remaining color tag in showerquarks (if
                // present)
                for (int ishq = 0; ishq < showerquarks.num(); ++ishq) {
                  if (!showerquarks[ishq].is_used() &&
                      showerquarks[ishq].col() == considering[1].acol()) {
                    showerquarks[ishq].col(considering[0].col()); /*break;*/
                  }
                }
                for (int ishq = 0; ishq < HH_showerptns.num(); ++ishq) {
                  if (HH_showerptns[ishq].col() == considering[1].acol()) {
                    HH_showerptns[ishq].col(considering[0].col()); /*break;*/
                  }
                }
              } else {
                HH_thermal[-element[1]].acol(considering[0].col());
              }
            } else if (considering[0].col() > 0) {  // MAT + LBT/THERM
              int loc =
                  findcloserepl(considering[1], perm2[q2], true, true,
                                showerquarks, HH_thermal);  // functon to find
              if (loc == 999999999) {
                // std::cout <<endl<<"Warning : extra parton used for string
                // repair!!"<<endl;
                HHparton fakep = considering[1];
                fakep.id(-fakep.id());
                fakep.col(considering[0].col());
                fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                fakep.px(0.);
                fakep.py(0.);
                fakep.pz(0.);
                fakep.e(fakep.mass());
                fakep.is_fakeparton(true);
                Extraparton.add(fakep);
                // somewhere, we need to make for loop to toss all partons to
                // remnants list.
              } else if (loc < 0) {
                HH_thermal[-loc - 1].col(considering[0].col());
              } else {
                showerquarks[loc - 1].col(considering[0].col());
              }
            } else if (considering[1].acol() > 0) {  // LBT + MAT
              int loc = findcloserepl(considering[0], perm1[q1] + 1, true, true,
                                      showerquarks, HH_thermal);
              if (loc == 999999999) {
                // std::cout <<endl<<"Warning : extra parton used for string
                // repair!!"<<endl;
                HHparton fakep = considering[0];
                fakep.id(-fakep.id());
                fakep.acol(considering[1].acol());
                fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                fakep.px(0.);
                fakep.py(0.);
                fakep.pz(0.);
                fakep.e(fakep.mass());
                fakep.is_fakeparton(true);
                Extraparton.add(fakep);
                // somewhere, we need to make for loop to toss all partons to
                // remnants list.
              } else if (loc < 0) {
                HH_thermal[-loc - 1].acol(considering[1].acol());
              } else {
                showerquarks[loc - 1].acol(considering[1].acol());
              }
            }
          } else if (considering[0].id() < 0 &&
                     considering[1].id() >
                         0) {  // case of first parton is q-bar and second is q
            // std::cout <<endl<<"chosen partons are "<<
            // int(considering[0].id()) << " and " << int(considering[1].id())
            // << endl <<endl; std::cout <<endl<<"and their color tag is "<<
            // int(considering[0].acol()) << " and " <<
            // int(considering[1].col()) << endl <<endl; std::cout <<"color
            // correction implemented as Follows: "<< considering[0].acol()<<" =
            // "<<considering[1].col()<<endl;

            // Overwrite non-dominant colortags in tempjunctions with dominant
            // ones qbar-q case
            if (considering[0].acol() != 0 && considering[1].col() != 0) {
              for (int ijunc = 0; ijunc < Tempjunctions.size(); ijunc++) {
                if (Tempjunctions.at(ijunc).at(1).at(1) ==
                    considering[1].col()) {
                  Tempjunctions.at(ijunc).at(1).pop_back();
                  Tempjunctions.at(ijunc).at(1).push_back(
                      considering[0].acol());
                }
                if (Tempjunctions.at(ijunc).at(2).at(1) ==
                    considering[1].col()) {
                  Tempjunctions.at(ijunc).at(2).pop_back();
                  Tempjunctions.at(ijunc).at(2).push_back(
                      considering[0].acol());
                }
                if (Tempjunctions.at(ijunc).at(3).at(1) ==
                    considering[1].col()) {
                  Tempjunctions.at(ijunc).at(3).pop_back();
                  Tempjunctions.at(ijunc).at(3).push_back(
                      considering[0].acol());
                }
              }
            }

            // before changing the tags, revise the MesonrecoMatrix elements!
            int coltag1 = considering[0].acol();
            int coltag2 = considering[1].col();
            if (coltag1 > 0 && coltag2 > 0 && coltag1 <= limit &&
                coltag2 <= limit) {
              std::vector<int>::iterator I1 = std::find(
                  IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
              std::vector<int>::iterator I2 = std::find(
                  IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
              int loc1 = std::distance(IndiceForColFin.begin(), I1);
              int loc2 = std::distance(IndiceForColFin.begin(),
                                       I2);  // set up for find matrix indices
                                             // corresponding to the color tags

              MesonrecoMatrix1.at(loc1).at(loc2) = 1;
              MesonrecoMatrix1.at(loc2).at(loc1) = 1;  // Matrix revised

              for (int iM = 0; iM < MesonrecoMatrix1[0].size(); ++iM) {
                if (MesonrecoMatrix1[loc1][iM] == 1) {
                  MesonrecoMatrix1[loc2][iM] = 1;
                }
                if (MesonrecoMatrix1[loc2][iM] == 1) {
                  MesonrecoMatrix1[loc1][iM] = 1;
                }
                if (MesonrecoMatrix1[iM][loc1] == 1) {
                  MesonrecoMatrix1[iM][loc2] = 1;
                }
                if (MesonrecoMatrix1[iM][loc2] == 1) {
                  MesonrecoMatrix1[iM][loc1] = 1;
                }
                if (MesonrecoMatrix1[loc1][iM] == 0) {
                  MesonrecoMatrix1[loc2][iM] = 0;
                }
                if (MesonrecoMatrix1[loc2][iM] == 0) {
                  MesonrecoMatrix1[loc1][iM] = 0;
                }
                if (MesonrecoMatrix1[iM][loc1] == 0) {
                  MesonrecoMatrix1[iM][loc2] = 0;
                }
                if (MesonrecoMatrix1[iM][loc2] == 0) {
                  MesonrecoMatrix1[iM][loc1] = 0;
                }
              }
            }
            /*
            std::cout <<endl<<"Revised Matrix is same as below"<<endl;
            for(int irow=0; irow < IndiceForColFin.size(); irow++){
              for(int icol=0; icol < IndiceForColFin.size(); icol++){
                std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
              }
              std::cout <<endl<<endl;
            }*/
            if (considering[0].acol() > 0 &&
                considering[1].col() > 0) {  // MAT + MAT
              if (perm2[q2] > 0) {
                HH_showerptns[showerquarks[element[1]].par()].col(
                    considering[0]
                        .acol());  // now color tags from both partons are same
                // also need to set the remaining color tag in showerquarks (if
                // present)
                for (int ishq = 0; ishq < showerquarks.num(); ++ishq) {
                  if (!showerquarks[ishq].is_used() &&
                      showerquarks[ishq].acol() == considering[1].col()) {
                    showerquarks[ishq].acol(considering[0].acol()); /*break;*/
                  }
                }
                for (int ishq = 0; ishq < HH_showerptns.num(); ++ishq) {
                  if (HH_showerptns[ishq].acol() == considering[1].col()) {
                    HH_showerptns[ishq].acol(considering[0].acol()); /*break;*/
                  }
                }
              } else {
                HH_thermal[-element[1]].col(considering[0].acol());
              }
            } else if (considering[0].acol() > 0) {  // MAT + LBT/THERM
              int loc = findcloserepl(considering[1], perm2[q2], true, true,
                                      showerquarks, HH_thermal);
              if (loc == 999999999) {
                // std::cout <<endl<<"Warning : extra parton used for string
                // repair!!"<<endl;
                HHparton fakep = considering[1];
                fakep.id(-fakep.id());
                fakep.acol(considering[0].acol());
                fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                fakep.px(0.);
                fakep.py(0.);
                fakep.pz(0.);
                fakep.e(fakep.mass());
                fakep.is_fakeparton(true);
                Extraparton.add(fakep);
                // somewhere, we need to make for loop to toss all partons to
                // remnants list.
              } else if (loc > 0) {
                showerquarks[loc - 1].acol(considering[0].acol());
              } else if (loc < 0) {
                HH_thermal[-loc - 1].acol(considering[0].acol());
              }
            } else if (considering[1].col() > 0) {  // LBT + MAT
              int loc = findcloserepl(considering[0], perm1[q1] + 1, true, true,
                                      showerquarks, HH_thermal);
              if (loc == 999999999) {
                // std::cout <<endl<<"Warning : extra parton used for string
                // repair!!"<<endl;
                HHparton fakep = considering[0];
                fakep.id(-fakep.id());
                fakep.col(considering[1].col());
                fakep.mass(pythia.particleData.m0(std::abs(fakep.id())));
                fakep.px(0.);
                fakep.py(0.);
                fakep.pz(0.);
                fakep.e(fakep.mass());
                fakep.is_fakeparton(true);
                Extraparton.add(fakep);
                // somewhere, we need to make for loop to toss all partons to
                // remnants list.
              } else if (loc > 0) {
                showerquarks[loc - 1].col(considering[1].col());
              } else if (loc < 0) {
                HH_thermal[-loc - 1].col(considering[1].col());
              }
            }
            // now color tags from both partons are same
          }

          // now we're forming the hadron
          HHhadron formedhadron;
          // setting the hadron values: is a recombined hadron, mass, and
          // parents - setting the par3 flag to 999999 to denote the lack of a
          // 3rd parent parton
          formedhadron.is_recohad(true);
          formedhadron.mass(p_MCM[0].t() + p_MCM[1].t());
          // formedhadron.par1 = element[0]; formedhadron.par2 = element[1];
          // formedhadron.par3 = 999999;
          formedhadron.add_par(showerquarks[element[0]].par());
          if (perm2[q2] > 0) {
            formedhadron.add_par(showerquarks[element[1]].par());
          } else {
            formedhadron.add_par(element[1]);
          }

          // now setting if meson is in excited state
          if (angular_qnum + radial_qnum > 0) {
            formedhadron.is_excited(true);
          }

          // setting if there are any thermal partons used to make the hadron:
          // (setting the '0' thermal parent to -99999 so that it doesn't
          // conflict with '0' shower parton)
          if (perm2[q2] > 0) { /*is sh-sh*/
            formedhadron.is_shsh(true);
          } else if (perm2[q2] < 0) { /*is sh-th*/
            formedhadron.is_shth(true);
            if (element[1] == 0) {
              formedhadron.parents[1] = -99999;
            }
          }

          // setting hadron position and momentum vectors
          Pmeson.Set(Pmeson.x(), Pmeson.y(), Pmeson.z(),
                     sqrt(Pmeson.x() * Pmeson.x() + Pmeson.y() * Pmeson.y() +
                          Pmeson.z() * Pmeson.z() +
                          formedhadron.mass() * formedhadron.mass()));
          formedhadron.pos(pos_lab);
          formedhadron.P(Pmeson);

          // setting hadron color tags (for tracing colors of the constituent
          // partons) will need to update to reflect color tags given to random
          // (thermal/lbt) partons [**]
          if (considering[0].id() > 0) {
            formedhadron.add_col(considering[0].col());
            formedhadron.add_col(considering[1].acol());
          } else {
            formedhadron.add_col(considering[1].col());
            formedhadron.add_col(considering[0].acol());
          }

          // need to choose *what* hadron we've formed... base this on the
          // parton id's, mass, & if excited
          set_meson_id(considering, formedhadron, angular_qnum, radial_qnum);

          // need to add the hadron to the collection
          HH_hadrons.add(formedhadron);

          // now that we've formed the hadron, need to set ALL (both) the
          // 'considering' flags to used
          showerquarks[element[0]].status(1);
          showerquarks[element[0]].is_used(true);
          if (perm2[q2] > 0) {
            showerquarks[element[1]].status(1);
            showerquarks[element[1]].is_used(true);
          } else {
            HH_thermal[-element[1]].status(1);
            HH_thermal[-element[1]].is_used(true);
          }

          // now that we've formed the hadron, break to first loop here!
          madehadron = true;
          considering.clear();
          break;
        }
      }

      // if we've formed a hadron - break to first loop
      if (madehadron) {
        break;
      }

      // since we CAN'T form a baryon on this try, need to revert the second
      // quark 'considering' used flag
      if (perm2[q2] > 0) {
        showerquarks[element[1]].status(0);
      } else {
        HH_thermal[-element[1]].status(0);
      }

      // and remove the second entry in considering (putting in if statement in
      // case we tried and failed to make a hadron;
      considering.partons.pop_back();
    }
    // if we've formed a hadron - continue to next parton in first loop...
    if (madehadron) {
      continue;
    }

    // since we've not formed a hadron with the first quark, need to revert the
    // first quark 'considering' used flag only need to reset these for shower
    // quarks as the first quark cannot be a thermal quark and remove the
    // first(only) entry in considering
    showerquarks[element[0]].status(0);
    considering.partons.pop_back();
  }

  // all possibilities have been considered, all used quark flags for
  // showerquarks are set appropriately time for cleanup

  // set the used quarks in the original shower to reflect that they were used
  // in reco module and that they were actually used set the fully used gluons
  // in the original shower to reflect that they were used in reco module and
  // that they were actually used set the partially used gluons in the original
  // shower to reflect that they were used in reco module and that they were
  // actually used give used quarks and fully used gluons a status of '1'; give
  // partially used gluons a status of '-1' stick all unused quarks and
  // 'completely' unused gluons into remnants (make sure to set the parents to
  // the original shower partons appropriately) write updated string information
  // into shower, so that it is set properly in remnants for the thermal array,
  // all the used flags should already be set (and have a status of '1') if this
  // is used in a loop (as the original version should be doing) - those will
  // have to be reset before reuse otherwise, this is perfect for medium
  // feedback

  // using quarks in showerquarks to set the flags appropriately for partons in
  // shower
  for (int i = 0; i < showerquarks.num(); ++i) {
    // if we have a quark in the original shower
    if ((std::abs(HH_showerptns[showerquarks[i].par()].id()) <= 5) &&
        (showerquarks[i].is_used())) {
      HH_showerptns[showerquarks[i].par()].is_used(true);
      HH_showerptns[showerquarks[i].par()].status(1);
      HH_showerptns[showerquarks[i].par()].used_reco(true);
    }
    // if this quark is from a split gluon in the original shower
    else if (std::abs(HH_showerptns[showerquarks[i].par()].id()) == 21 &&
             showerquarks[i].is_used()) {
      // if this is the first used quark in a splitting, set the parent gluon to
      // used, and status of -1
      if (HH_showerptns[showerquarks[i].par()].status() == -99) {
        HH_showerptns[showerquarks[i].par()].status(-1);
        HH_showerptns[showerquarks[i].par()].is_used(true);
        HH_showerptns[showerquarks[i].par()].used_reco(true);
      }
      // if this is the second (last) used quark in a splitting, set the status
      // to 1
      else if (HH_showerptns[showerquarks[i].par()].status() == -1) {
        HH_showerptns[showerquarks[i].par()].status(1);
      }
      // remove this check if it never throws.
      else {
        JSWARN << "SOMETHING HAS GONE VERY WRONG WITH REFORMING GLUON IN POS: "
               << showerquarks[i].par();
        int val;
        showerquarks[i].par(0);
      }
    }
  }
  // need to run back through shower; if there are any gluons that didn't get
  // used at all in the shower - restore them (and output to remnants below)
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    if (HH_showerptns[i].status() == -99) {
      HH_showerptns[i].is_decayedglu(false);
      HH_showerptns[i].status(0);
    }
  }

  // need to update string information in shower from showerquarks
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    if (!HH_showerptns[i].is_used()) {
      for (int j = 0; j < showerquarks.num(); ++j) {
        if (showerquarks[j].par() == i && !showerquarks[j].is_used()) {
          HH_showerptns[i].string_id(showerquarks[j].string_id());
          HH_showerptns[i].is_strendpt(showerquarks[j].is_strendpt());
          HH_showerptns[i].pos_str(showerquarks[j].pos_str());
          HH_showerptns[i].endpt_id(showerquarks[j].endpt_id());
          break;
        }
      }
    }
  }
  // now all the partons in shower have had flags appropriately set; the
  // 'partially' used gluons have a status of -1

  // TODO: here we need new function to find thermal sibling for gluon
  // sol 1 : gluon loops.
  // sol 2 : decay gluon find pairs.{concern : energy conservation violated}
  // sol 3 : find two theraml siblings for gluon. declare fincloserepl twice.
  // first pick quark and antiquarks to the 2nd.--> this is chosen
  for (int i = 0; i < HH_showerptns.num(); i++) {
    if (HH_showerptns[i].is_used()) {
      continue;
    }
    if (HH_showerptns[i].id() == 21 && HH_showerptns[i].col() == 0 &&
        HH_showerptns[i].acol() == 0) {
      int sel_out[2] = {0, 0};
      findcloserepl_glu(HH_showerptns[i], i + 1, true, true, HH_showerptns,
                        HH_thermal, sel_out);
      if (sel_out[0] == 999999999 || sel_out[1] == 999999999) {
        HHparton fakeg = HH_showerptns[i];
        fakeg.id(21);
        fakeg.acol(++maxtag);
        fakeg.col(++maxtag);
        fakeg.mass(1e-6);
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakeg.px(fake_pT * cos(fake_phi));
        fakeg.py(fake_pT * sin(fake_phi));
        fakeg.pz(0.);
        fakeg.e(std::sqrt(fakeg.px() * fakeg.px() + fakeg.py() * fakeg.py() +
                          fakeg.pz() * fakeg.pz() +
                          fakeg.mass() * fakeg.mass()));
        fakeg.orig(-1);
        fakeg.is_remnant(true);
        fakeg.is_fakeparton(true);
        Extraparton.add(fakeg);
        HH_showerptns[i].col(fakeg.acol());
        HH_showerptns[i].acol(fakeg.col());
      } else if (sel_out[0] > 0) {
        HH_showerptns[sel_out[0] - 1].col(++maxtag);
        HH_showerptns[i].acol(maxtag);
      } else if (sel_out[0] < 0) {
        HH_thermal[-sel_out[0] - 1].col(++maxtag);
        HH_showerptns[i].acol(maxtag);
      }

      if (sel_out[1] > 0 && sel_out[1] != 999999999 &&
          sel_out[0] != 999999999) {
        HH_showerptns[sel_out[1] - 1].acol(++maxtag);
        HH_showerptns[i].col(maxtag);
      } else if (sel_out[1] < 0 && sel_out[1] != 999999999 &&
                 sel_out[0] != 999999999) {
        HH_thermal[-sel_out[1] - 1].acol(++maxtag);
        HH_showerptns[i].col(maxtag);
      }
    } else if (HH_showerptns[i].id() == 21 &&
               (HH_showerptns[i].col() == 0 ||
                HH_showerptns[i].acol() ==
                    0)) {  // gluon that needs a single quark to fix
      if (HH_showerptns[i].col() == 0) {
        HH_showerptns[i].id(1);  // need to temporarily pretend gluon is a
                                 // quark, swap back id after...
        int loc = findcloserepl(HH_showerptns[i], i + 1, true, true,
                                HH_showerptns, HH_thermal);
        HH_showerptns[i].id(21);
        if (loc == 999999999) {
          double fid = (ran() > 0.5) ? -1 : -2;
          HHparton fakep = HH_showerptns[i];
          fakep.id(fid);
          fakep.acol(++maxtag);
          double dir = (ran() < 0.5) ? 1. : -1.;
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          Extraparton.add(fakep);
        } else if (loc > 0) {
          HH_showerptns[loc - 1].acol(++maxtag);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].acol(++maxtag);
        }

        HH_showerptns[i].col(maxtag);
      } else if (HH_showerptns[i].acol() == 0) {
        HH_showerptns[i].id(-1);  // need to temporarily pretend gluon is an
                                  // antiquark, swap back id after...
        int loc = findcloserepl(HH_showerptns[i], i + 1, true, true,
                                HH_showerptns, HH_thermal);
        HH_showerptns[i].id(21);
        if (loc == 999999999) {
          double fid = (ran() > 0.5) ? 1 : 2;
          HHparton fakep = HH_showerptns[i];
          fakep.id(fid);
          fakep.col(++maxtag);
          double dir = (ran() < 0.5) ? 1. : -1.;
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          Extraparton.add(fakep);
        } else if (loc > 0) {
          HH_showerptns[loc - 1].col(++maxtag);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].col(++maxtag);
        }

        HH_showerptns[i].acol(maxtag);
      }
    } else if (HH_showerptns[i].id() > 0 && HH_showerptns[i].col() == 0) {
      int loc = findcloserepl(HH_showerptns[i], i + 1, true, true,
                              HH_showerptns, HH_thermal);
      if (loc == 999999999) {
        double fid = (ran() > 0.5) ? -1 : -2;
        HHparton fakep = HH_showerptns[i];
        fakep.id(fid);
        fakep.acol(++maxtag);
        double dir = (ran() < 0.5) ? 1. : -1.;
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        Extraparton.add(fakep);
      } else if (loc > 0) {
        HH_showerptns[loc - 1].acol(++maxtag);
      } else if (loc < 0) {
        HH_thermal[-loc - 1].acol(++maxtag);
      }

      HH_showerptns[i].col(maxtag);
    } else if (HH_showerptns[i].id() < 0 && HH_showerptns[i].acol() == 0) {
      int loc = findcloserepl(HH_showerptns[i], i + 1, true, true,
                              HH_showerptns, HH_thermal);
      if (loc == 999999999) {
        double fid = (ran() > 0.5) ? 1 : 2;
        HHparton fakep = HH_showerptns[i];
        fakep.id(fid);
        fakep.col(++maxtag);
        double dir = (ran() < 0.5) ? 1. : -1.;
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        Extraparton.add(fakep);
      } else if (loc > 0) {
        HH_showerptns[loc - 1].col(++maxtag);
      } else if (loc < 0) {
        HH_thermal[-loc - 1].col(++maxtag);
      }

      HH_showerptns[i].acol(maxtag);
    }
  }
  // sticking all unused partons into remnants; keeping order intact (a
  // partially used gluon is replaced with it's unused quark)
  for (int i = 0; i < HH_showerptns.num(); ++i) {
    // if unused parton, write into remnants
    if (HH_showerptns[i].status() == 0) {
      HH_remnants.add(HH_showerptns[i]);
      HH_remnants[HH_remnants.num() - 1].par(i);
      HH_showerptns[i].is_remnant(true);
    }
    // if 'partially' used gluon, write unused daughter quark into remnants
    else if (HH_showerptns[i].status() == -1) {
      // finding the unused quark for this gluon and adding it to remnants (have
      // to loop over as we only keep track of parents, not daughters)
      for (int j = 0; j < showerquarks.num(); ++j) {
        if (showerquarks[j].par() == i && !showerquarks[j].is_used()) {
          if (showerquarks[j].col() == 0 &&
              showerquarks[j].id() > 0) {  // quark with zero col tag is left,
            int loc = findcloserepl(showerquarks[j], j + 1, true, true,
                                    HH_showerptns, HH_thermal);
            if (loc == 999999999) {
              double fid = (ran() > 0.5) ? -1 : -2;
              HHparton fakep = showerquarks[j];
              fakep.id(fid);
              fakep.acol(++maxtag);
              double dir = (ran() < 0.5) ? 1. : -1.;
              double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
              double fake_phi = 2. * 3.14159265358979 * ran();
              fakep.px(fake_pT * cos(fake_phi));
              fakep.py(fake_pT * sin(fake_phi));
              if (number_p_fake == 0 || number_p_fake == 2) {
                fakep.pz(-p_fake);
              } else if (number_p_fake == 1 || number_p_fake == 3) {
                fakep.pz(p_fake);
              } else {
                fakep.pz(0.);
              }
              number_p_fake++;
              fakep.e(std::sqrt(
                  fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                  fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
              fakep.orig(-1);
              fakep.is_remnant(true);
              fakep.is_fakeparton(true);
              Extraparton.add(fakep);
            } else if (loc > 0) {
              HH_showerptns[loc - 1].acol(++maxtag);
            } else if (loc < 0) {
              HH_thermal[-loc - 1].acol(++maxtag);
            }

            showerquarks[j].col(maxtag);
          } else if (showerquarks[j].acol() == 0 &&
                     showerquarks[j].id() <
                         0) {  // antiquark with zero col tag is left,
            int loc = findcloserepl(showerquarks[j], j + 1, true, true,
                                    HH_showerptns, HH_thermal);
            if (loc == 999999999) {
              double fid = (ran() > 0.5) ? 1 : 2;
              HHparton fakep = showerquarks[j];
              fakep.id(fid);
              fakep.col(++maxtag);
              double dir = (ran() < 0.5) ? 1. : -1.;
              double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
              double fake_phi = 2. * 3.14159265358979 * ran();
              fakep.px(fake_pT * cos(fake_phi));
              fakep.py(fake_pT * sin(fake_phi));
              if (number_p_fake == 0 || number_p_fake == 2) {
                fakep.pz(-p_fake);
              } else if (number_p_fake == 1 || number_p_fake == 3) {
                fakep.pz(p_fake);
              } else {
                fakep.pz(0.);
              }
              number_p_fake++;
              fakep.e(std::sqrt(
                  fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                  fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
              fakep.orig(-1);
              fakep.is_remnant(true);
              fakep.is_fakeparton(true);
              Extraparton.add(fakep);
            } else if (loc > 0) {
              HH_showerptns[loc - 1].col(++maxtag);
            } else if (loc < 0) {
              HH_thermal[-loc - 1].col(++maxtag);
            }

            showerquarks[j].acol(maxtag);
          }

          HH_remnants.add(showerquarks[j]);
          break;
        }
      }
      HH_showerptns[i].is_remnant(true);
    }
  }

  // appending the thermal partons used in the string repair functionality into
  // remnants
  for (int i = 0; i < HH_thermal.num(); ++i) {
    // if this thermal parton is ... then add it to the remnants collection
    if (HH_thermal[i].is_used()) {
      HH_thermal[i].status(1);
      HH_thermal[i].used_reco(true);
      continue;
    }
    if (HH_thermal[i].col() > 0 || HH_thermal[i].acol() > 0) {
      // std::cout << "Thermal parton added to remnants: " << HH_thermal[i].id()
      // << "," << HH_thermal[i].col() << "," << HH_thermal[i].acol() <<
      // std::endl;
      HH_remnants.add(HH_thermal[i]);
      HH_remnants[HH_remnants.num() - 1].par(-i - 1);
      HH_thermal[i].is_remnant(true);
      HH_thermal[i].is_used(true);
    }
  }

  for (int i = 0; i < Extraparton.num(); ++i) {
    HH_remnants.add(Extraparton[i]);
    Extraparton[i].is_used(
        true);  // use the extrapartons during hadronizations of negative
                // partons (background subtraction)

    if (Extraparton[i].id() == 21) {
      HHparton fakep1 = Extraparton[i];
      fakep1.id(1);
      fakep1.acol(0);
      fakep1.px(fakep1.px() / 2.);
      fakep1.py(fakep1.py() / 2.);
      fakep1.pz(fakep1.pz() / 2.);
      fakep1.mass(xmq);
      fakep1.e(std::sqrt(fakep1.px() * fakep1.px() + fakep1.py() * fakep1.py() +
                         fakep1.pz() * fakep1.pz() +
                         fakep1.mass() * fakep1.mass()));

      HHparton fakep2 = Extraparton[i];
      fakep2.id(-1);
      fakep2.col(0);
      fakep2.px(fakep2.px() / 2.);
      fakep2.py(fakep2.py() / 2.);
      fakep2.pz(fakep2.pz() / 2.);
      fakep2.mass(xmq);
      fakep2.e(std::sqrt(fakep2.px() * fakep2.px() + fakep2.py() * fakep2.py() +
                         fakep2.pz() * fakep2.pz() +
                         fakep2.mass() * fakep2.mass()));

      HH_recomb_extrapartons.add(fakep1);
      HH_recomb_extrapartons.add(fakep2);
    } else {
      HH_recomb_extrapartons.add(Extraparton[i]);
    }
  }

  // hadrons have been recombined, and output to hadron collection
  // remnants have been collected, and output to remnant collection
  // shower partons have all been updated appropriately
  // thermal partons have all been updated appropriately - make sure that the
  // thermal partons are reset before the recomb module is called again...

  // end of recombination routine
}

// sets id of formed baryon based on quark content, mass of quark system, and if
// the baryon formed into an excited state
void HybridHadronization::set_baryon_id(parton_collection& qrks,
                                        HHhadron& had) {
  // assigning quark_ids in descending order to construct baryon id
  int id[3] = {qrks[0].id(), qrks[1].id(), qrks[2].id()};
  std::sort(id, id + 3, [](int a, int b) { return std::abs(a) > std::abs(b); });

  /*//http://pdg.lbl.gov/2017/listings/contents_listings.html
        double mdelta, msigma, mlambda, mxi, mdelE, msigmaE, mxiE, momega;
        mdelta  = 1.232;
        msigma  = 1.190;
        msigmaE = 1.382;
        mlambda = 1.115;
        mxi     = 1.315;
        mxiE    = 1.530;
        momega  = 1.672;
  //	^^^^^ MIGHT NOT NEED MASSES BUT KEEP FOR NOW ^^^^^^
  */
  // id and mass list:
  // ground state: (111, 222, 333 configurations prohibited...)
  //	2212, 2112 - .938  p,n
  //	3222, 3212, 3112 - 1.190   sigma;  3122 - 1.115 lambda
  //	3322, 3312 - 1.315  xi
  //
  // excited:
  //	2224, 2214, 2114, 1114 - 1.232  delta
  //	3224, 3214, 3114 - 1.190   sigma;
  //	3324, 3314 - 1.315  xi
  //	3334 - 1.672 omega
  //	^^^^^ MIGHT NOT NEED BUT KEEP FOR NOW ^^^^^^

  if (id[0] == id[1] && id[0] == id[2]) {
    had.id(1000 * std::abs(id[0]) + 100 * std::abs(id[1]) +
           10 * std::abs(id[2]) + 4);
  }  // J=3/2 only
  else if (id[0] == id[1] || id[0] == id[2] || id[1] == id[2]) {
    if (ran() > 0.333) {
      had.id(1000 * std::abs(id[0]) + 100 * std::abs(id[1]) +
             10 * std::abs(id[2]) + 4);
    }  // J=3/2
    else {
      had.id(1000 * std::abs(id[0]) + 100 * std::abs(id[1]) +
             10 * std::abs(id[2]) + 2);
    }
  } else {
    double prb = ran();
    if (prb > 0.333) {
      had.id(1000 * std::abs(id[0]) + 100 * std::abs(id[1]) +
             10 * std::abs(id[2]) + 4);
    }  // J=3/2
    else if (prb > 0.166) {
      had.id(1000 * std::abs(id[0]) + 100 * std::abs(id[1]) +
             10 * std::abs(id[2]) + 2);
    }  // J=1/2 higher mass
    else {
      had.id(1000 * std::abs(id[0]) + 10 * std::abs(id[1]) +
             100 * std::abs(id[2]) + 2);
    }  // J=1/2 lower mass
  }    // note the swap of quark index in last line

  // This would be how to make excited N, Delta and Lambda. RIGHT NOW DISABLED
  // IN FIRST IF STATEMENT
  // there are no excited baryon codes of this form; there are only spin excited
  // states in PYTHIA's ParticleData class
  // if(false && had.is_excited){
  //	if(had.id == 2212 || had.id == 2112 || had.id == 2214 || had.id == 2114
  //|| had.id == 2224 || had.id ==1114 || had.id == 3122){ 		had.id
  //+= 100000;  // This needs to be adjusted to make the actual excited baryon
  // codes
  //	}
  // }

  had.id(had.id() * (id[0] < 0 ? -1 : 1));
  return;
}

// sets id of formed meson based on quark content, mass of quark system, and if
// the meson formed into an excited state
void HybridHadronization::set_meson_id(parton_collection& qrks, HHhadron& had,
                                       int l, int k) {
  // assigning quark_ids in descending order to construct meson id
  int id[2] = {qrks[0].id(), qrks[1].id()};
  if (std::abs(qrks[1].id()) > std::abs(qrks[0].id())) {
    std::swap(id[0], id[1]);
  }

  //  Don't need the following any more
  //	double mass_pi0, mass_eta, mass_omegam, mass_etap, mass_phi; //mass_rho;
  //	mass_pi0    = 0.1349770;
  //	mass_eta    = 0.547862;
  //	mass_rho    = 0.7690;
  //	mass_omegam = 0.78265;
  //	mass_etap   = 0.95778;
  //	mass_phi    = 1.019460;
  //// MASSES ABOVE NOT REALLY NEEDED IN SIMPLE QUARK MODEL BASED APPROACH

  int baseid = 0;

  // Determine spin quantum number statistically
  int spin_qnum = (ran() > 0.25) ? 1 : 0;

  // Determine total angular momentum statistically
  int j_qnum = spin_qnum;
  if (l > 0 && spin_qnum == 1) {
    double random = ran();

    if ((2. * l - 1.) / (3. * (2. * l + 1.)) >= random) {
      j_qnum = l - 1;
    } else if (1. / 3. + (2. * l - 1.) / (3. * (2. * l + 1.)) >= random) {
      j_qnum = l;
    } else {
      j_qnum = l + 1;
    }
  } else if (spin_qnum == 0) {
    j_qnum = l;
  }

  // All meson quantum numbers are now known

  int basesign = 1;
  // if isospin I3=0
  if (id[0] == -id[1]) {
    if (ran() > 0.25) {  // spin triplet
      if (std::abs(id[0]) == 5) {
        baseid = 550;
      }  // Upsilon
      if (std::abs(id[0]) == 4) {
        baseid = 440;
      }  // J/psi
      if (std::abs(id[0]) == 3) {
        baseid = 330;
      }  // phi
      if ((std::abs(id[0]) == 1) || (std::abs(id[0]) == 2)) {
        if (ran() > 0.5) {
          baseid = 220;
        } else {
          baseid = 110;
        }
      }       // omega and rho
    } else {  // spin singlet
      if (std::abs(id[0]) == 5) {
        baseid = 550;
      }  // etaB
      if (std::abs(id[0]) == 4) {
        baseid = 440;
      }  // etac
      if (std::abs(id[0]) == 3) {
        if (ran() > 0.666) {
          baseid = 330;
        } else {
          baseid = 220;
        }
      }  // eta' and eta
      if (std::abs(id[0]) < 3) {
        double prb = ran();
        if (prb > 0.5) {
          baseid = 110;
        }  // pi0
        else if (prb > 0.333) {
          baseid = 220;
        }  // eta
        else {
          baseid = 330;
        }  // eta'
      }
    }
  }
  // if isospin I3 not 0
  else {
    baseid = 100 * std::abs(id[0]) + 10 * std::abs(id[1]);
    if (id[0] % 2 == 0) {
      basesign = 2 * std::signbit(-id[0]) - 1;
    } else {
      basesign = 2 * std::signbit(id[0]) - 1;
    }
  }

  // Put everything together: first digit (from the right)
  if (j_qnum < 5) {
    baseid += 2 * j_qnum + 1;
  } else {
    baseid += 8;
  }

  // 5th digit
  if (l > 0 && spin_qnum == 0) {
    baseid += 10000;
  } else if (l > 0 && spin_qnum == 1 && l == j_qnum) {
    baseid += 20000;
  } else if (l > 1 && spin_qnum == 1 && l == j_qnum + 1) {
    baseid += 30000;
  } else if (l == 1 && spin_qnum == 1 && l == j_qnum + 1) {
    baseid += 10000;
  }

  // 6th digit
  if (k < 10) {
    baseid += k * 100000;
  } else {
    baseid += 900000;
  }

  // If we don't want Goldstone bosons convert them to vector mesons here
  if (!goldstonereco) {
    if (baseid == 211) {
      baseid = 213;
    }  // fixed pi+- -> rho+-
    if (baseid == 311) {
      baseid = 313;
    }  // fixed for K+- -> K*+-
    if (baseid == 321) {
      baseid = 323;
    }  // fixed for K0 -> K*0
    if (baseid == 111) {
      baseid = 113;
    }  // fixed for pi0 -> rho0
  }

  baseid *= basesign;

  had.id(baseid);
  return;
}

// gluon decay function
void HybridHadronization::gluon_decay(HHparton& glu, parton_collection& qrks) {
  HHparton q1, q2;
  // these are placeholders - might want to instead use values directly from
  // partons...
  double qmass, glu_e;

  // if set to already be not on-shell, but not initially set!
  // glu.mass = sqrt(glu.e()*glu.e() - glu.px()*glu.px() - glu.py()*glu.py() -
  // glu.pz()*glu.pz()); glu_e = glu.e();
  glu_e = sqrt(glu.mass() * glu.mass() + glu.px() * glu.px() +
               glu.py() * glu.py() + glu.pz() * glu.pz());

  // choosing qqbar ids (u, d, or s)
  // assuming that xms >= xmq (bad things *could* happen if not...)
  if (glu.mass() > 2. * xms) {
    //******** ratio = Gamma(g->ssbar)/Gamma(g->uubar, ddbar) ******
    double ratio = 0.5 *
                   sqrt((glu.mass() * glu.mass() - 4. * xms * xms) /
                        (glu.mass() * glu.mass() - 4. * xmq * xmq)) *
                   ((glu.mass() * glu.mass() + 2. * xms * xms) /
                    (glu.mass() * glu.mass() + 2. * xmq * xmq));
    double prob = ran();
    if (prob <= ratio / (1. + ratio)) {
      qmass = xms;
      q1.id(3);
      q2.id(-3);
      q1.mass(xms);
      q2.mass(xms);
    } else if ((prob > ratio / (1. + ratio)) &&
               (prob <= (0.5 + ratio) / (1. + ratio))) {
      qmass = xmq;
      q1.id(1);
      q2.id(-1);
      q1.mass(xmq);
      q2.mass(xmq);
    } else { /*if (prob > (0.5+ratio)/(1.+ratio))*/
      qmass = xmq;
      q1.id(2);
      q2.id(-2);
      q1.mass(xmq);
      q2.mass(xmq);
    }
  } else {
    double prob = ran();
    if (prob <= 0.5) {
      qmass = xmq;
      q1.id(1);
      q2.id(-1);
      q1.mass(xmq);
      q2.mass(xmq);
    } else {
      qmass = xmq;
      q1.id(2);
      q2.id(-2);
      q1.mass(xmq);
      q2.mass(xmq);
    }
  }

  // gluon velocity
  FourVector Betag;
  Betag.Set(-glu.px() / glu_e, -glu.py() / glu_e, -glu.pz() / glu_e, 0.);
  double sum2 =
      Betag.x() * Betag.x() + Betag.y() * Betag.y() + Betag.z() * Betag.z();
  if (sum2 < 1.) {
    Betag.Set(Betag.x(), Betag.y(), Betag.z(), 1. / sqrt(1. - sum2));
  }

  // setting the q-qbar momenta, starting in gluon rest frame
  FourVector Pq_CM, Pq1, Pq2;
  double pq = (glu.mass() > 2. * qmass)
                  ? sqrt(glu.mass() * glu.mass() / 4. - qmass * qmass)
                  : 0.;
  double theta = acos(1. - 2. * ran());
  double phi = 2 * pi * ran();
  Pq_CM.Set(pq * sin(theta) * cos(phi), pq * sin(theta) * sin(phi),
            pq * cos(theta), sqrt(qmass * qmass + pq * pq));
  Pq1 = HHboost(Betag, Pq_CM);
  Pq_CM.Set(-Pq_CM.x(), -Pq_CM.y(), -Pq_CM.z(), Pq_CM.t());
  Pq2 = HHboost(Betag, Pq_CM);
  q1.P(Pq1);
  q2.P(Pq2);
  q1.col(glu.col());
  q2.acol(glu.acol());

  // propagate quarks
  FourVector position1;
  FourVector position2;
  position1.Set(glu.x() + (Pq1.x() / qmass) * part_prop,
                glu.y() + (Pq1.y() / qmass) * part_prop,
                glu.z() + (Pq1.z() / qmass) * part_prop, glu.x_t() + part_prop);
  position2.Set(glu.x() + (Pq2.x() / qmass) * part_prop,
                glu.y() + (Pq2.y() / qmass) * part_prop,
                glu.z() + (Pq2.z() / qmass) * part_prop, glu.x_t() + part_prop);
  q1.pos(position1);
  q2.pos(position2);

  qrks.add(q1);
  qrks.add(q2);
}

// finding a thermal sibling for a thermal parton in therm
int HybridHadronization::findthermalsibling(int ithm,
                                            parton_collection& therm) {
  if (!therm[therm[ithm].sibling()].is_used() &&
      (therm[therm[ithm].sibling()].string_id() < 0) &&
      (ithm != therm[ithm].sibling())) {
    return therm[ithm].sibling();
  }
  int qrk_close = -1;
  double dist2min = 999999999999.;
  for (int i = 0; i < therm.num(); ++i) {
    if ((therm[ithm].id() * therm[i].id() > 0) || therm[i].is_used()) {
      continue;
    }
    double distnow = therm[ithm].posDif2(therm[i]) +
                     (therm[ithm].x_t() - therm[i].x_t()) *
                         (therm[ithm].x_t() - therm[i].x_t());
    if (distnow < dist2min) {
      qrk_close = i;
      dist2min = distnow;
    }
  }
  if (qrk_close == -1) {
    qrk_close = ithm;
  }
  return qrk_close;
}

int HybridHadronization::findcloserepl(HHparton ptn, int iptn, bool lbt,
                                       bool thm, parton_collection& sh_lbt,
                                       parton_collection& therm) {
  // should not happen
  if (iptn == 0 || ptn.id() == 21) {
    throw std::runtime_error(
        "Parton index is incorrect (should not be 0 or gluon)");
  }

  // if the parton is thermal, and we're only looking at thermal partons, and
  // the sibling was already found & not used, then return that.
  if ((iptn < 0) && thm && !lbt && !therm[ptn.sibling()].is_used() &&
      (((therm[ptn.sibling()].id() > 0) && (therm[ptn.sibling()].col() != 0)) ||
       ((therm[ptn.sibling()].id() < 0) &&
        (therm[ptn.sibling()].acol() != 0))) &&
      (iptn + 1 != ptn.sibling())) {
    return ptn.sibling();
  }

  // initializing vars
  int qrk_close = 999999999;
  double dist2min = 999999999999.;
  // checking thermal partons for closest parton
  if (thm) {
    for (int i = 0; i < therm.num(); ++i) {
      // if the parton's color tag is not 0, then skip (is already repairing a
      // string)
      if (((therm[i].id() > 0) && (therm[i].col() != 0)) ||
          ((therm[i].id() < 0) && (therm[i].acol() != 0))) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if ((ptn.id() * therm[i].id() > 0) || therm[i].is_used() ||
          ((iptn < 0) && (-i == iptn + 1))) {
        continue;
      }
      double distnow = ptn.posDif2(therm[i]) + (ptn.x_t() - therm[i].x_t()) *
                                                   (ptn.x_t() - therm[i].x_t());
      if (distnow < dist2min) {
        qrk_close = -i - 1;
        dist2min = distnow;
      }
    }
  }
  // checking lbt partons for closest parton
  if (lbt) {
    for (int i = 0; i < sh_lbt.num(); ++i) {
      // if gluon, skip
      if (sh_lbt[i].id() == 21) {
        continue;
      }
      // if the parton's color tag is not 0, then skip (not an lbt/martini?
      // parton)
      if (((sh_lbt[i].id() > 0) && (sh_lbt[i].col() != 0)) ||
          ((sh_lbt[i].id() < 0) && (sh_lbt[i].acol() != 0))) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if ((ptn.id() * sh_lbt[i].id() > 0) || sh_lbt[i].is_used() ||
          ((iptn > 0) && (i == iptn - 1))) {
        continue;
      }
      double distnow =
          ptn.posDif2(sh_lbt[i]) +
          (ptn.x_t() - sh_lbt[i].x_t()) * (ptn.x_t() - sh_lbt[i].x_t());
      if (distnow < dist2min) {
        qrk_close = i + 1;
        dist2min = distnow;
      }
    }
  }

  // positive return values indicate an lbt parton in the shower was found (the
  // (i-1)th parton) negative return values indicate a thermal parton was found
  // (the -(i+1)th parton)
  return qrk_close;
}

// version to handle finding a q-qbar pair for lbt gluons
void HybridHadronization::findcloserepl_glu(HHparton ptn, int iptn, bool lbt,
                                            bool thm, parton_collection& sh_lbt,
                                            parton_collection& therm,
                                            int sel_out[]) {
  // should not happen
  if (iptn == 0 || ptn.id() != 21) {
    throw std::runtime_error(
        "Parton index is incorrect (should not be 0, or anything other than "
        "gluon)");
  }
  if (ptn.is_thermal()) {
    throw std::runtime_error("Parton is thermal (should not be so)");
  }

  // initializing vars
  int qrk_close = 999999999;
  double dist2min = 999999999999.;
  // checking thermal partons for closest quark
  if (thm) {
    for (int i = 0; i < therm.num(); ++i) {
      // if the parton's color tag is not 0, then skip (is already repairing a
      // string), or skip if antiquark
      if (((therm[i].id() > 0) && (therm[i].col() != 0)) ||
          (therm[i].id() < 0)) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if (therm[i].is_used() || ((iptn < 0) && (-i == iptn + 1))) {
        continue;
      }
      double distnow = ptn.posDif2(therm[i]) + (ptn.x_t() - therm[i].x_t()) *
                                                   (ptn.x_t() - therm[i].x_t());
      if (distnow < dist2min) {
        qrk_close = -i - 1;
        dist2min = distnow;
      }
    }
  }
  // checking lbt partons for closest parton
  if (lbt) {
    for (int i = 0; i < sh_lbt.num(); ++i) {
      // if gluon, skip
      if (sh_lbt[i].id() == 21) {
        continue;
      }
      // if the parton's color tag is not 0, then skip (not an lbt/martini?
      // parton), or skip if antiquark
      if (((sh_lbt[i].id() > 0) && (sh_lbt[i].col() != 0)) ||
          (sh_lbt[i].id() < 0)) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if (sh_lbt[i].is_used() || ((iptn > 0) && (i == iptn - 1))) {
        continue;
      }
      double distnow =
          ptn.posDif2(sh_lbt[i]) +
          (ptn.x_t() - sh_lbt[i].x_t()) * (ptn.x_t() - sh_lbt[i].x_t());
      if (distnow < dist2min) {
        qrk_close = i + 1;
        dist2min = distnow;
      }
    }
  }

  // positive return values indicate an lbt parton in the shower was found (the
  // (i-1)th parton) negative return values indicate a thermal parton was found
  // (the -(i+1)th parton)
  sel_out[0] = qrk_close;

  qrk_close = 999999999;
  dist2min = 999999999999.;
  // checking thermal partons for closest antiquark
  if (thm) {
    for (int i = 0; i < therm.num(); ++i) {
      // if the parton's color tag is not 0, then skip (is already repairing a
      // string), or skip if quark
      if ((therm[i].id() > 0) ||
          ((therm[i].id() < 0) && (therm[i].acol() != 0))) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if (therm[i].is_used() || ((iptn < 0) && (-i == iptn + 1))) {
        continue;
      }
      double distnow = ptn.posDif2(therm[i]) + (ptn.x_t() - therm[i].x_t()) *
                                                   (ptn.x_t() - therm[i].x_t());
      if (distnow < dist2min) {
        qrk_close = -i - 1;
        dist2min = distnow;
      }
    }
  }
  // checking lbt partons for closest parton
  if (lbt) {
    for (int i = 0; i < sh_lbt.num(); ++i) {
      // if gluon, skip
      if (sh_lbt[i].id() == 21) {
        continue;
      }
      // if the parton's color tag is not 0, then skip (not an lbt/martini?
      // parton), or skip if quark
      if ((sh_lbt[i].id() > 0) ||
          ((sh_lbt[i].id() < 0) && (sh_lbt[i].acol() != 0))) {
        continue;
      }
      // if the parton is not a partner (eg. both quarks/antiquarks), or is
      // used, or is the same as the input parton, skip it.
      if (sh_lbt[i].is_used() || ((iptn > 0) && (i == iptn - 1))) {
        continue;
      }
      double distnow =
          ptn.posDif2(sh_lbt[i]) +
          (ptn.x_t() - sh_lbt[i].x_t()) * (ptn.x_t() - sh_lbt[i].x_t());
      if (distnow < dist2min) {
        qrk_close = i + 1;
        dist2min = distnow;
      }
    }
  }

  // positive return values indicate an lbt parton in the shower was found (the
  // (i-1)th parton) negative return values indicate a thermal parton was found
  // (the -(i+1)th parton)
  sel_out[1] = qrk_close;
}

// prepares remnant partons/strings for PYTHIA string hadronization
// sorts strings, ensures strings are in 'valid' configurations, assigns
// color/anticolor tags
// TODO: this might be where to use thermal partons to enforce color neutrality
void HybridHadronization::stringprep(parton_collection& SP_remnants,
                                     parton_collection& SP_prepremn,
                                     bool cutstr) {
  // dignostic measure
  /*std::cout <<endl<<"below is all the colors in the Tempjunction lists"<<endl;
  for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
    std::cout <<" { "<<Tempjunctions.at(ijunc).at(0).at(0)<<" ,
  "<<Tempjunctions.at(ijunc).at(1).at(1)<<" ,
  "<<Tempjunctions.at(ijunc).at(2).at(1)<<" ,
  "<<Tempjunctions.at(ijunc).at(3).at(1)<<" } "<<endl; std::cout <<"col/acol: ("
  << Tempjunctions.at(ijunc).at(1).at(0) << "," <<
  Tempjunctions.at(ijunc).at(2).at(0) << "," <<
  Tempjunctions.at(ijunc).at(3).at(0) << ")" << std::endl;
  }*/
  // Declare essential vectors for string repair
  vector<vector<vector<HHparton>>> JuncStructure;
  vector<vector<HHparton>>
      JuncLegs;  // vector of all junction ( Junction Num, Leg1, Leg2, Leg3)
  vector<HHparton> Leg1;  // partons in the legs for juncion formation
  vector<HHparton> Leg2;
  vector<HHparton> Leg3;
  vector<int> Legconsidering;  // address of partons in SP_remnants
  vector<vector<vector<int>>>
      IMStructure1;  // Intermediate structure for saving the indices in
                     // Tempjunction and Corresponding Leg Number, this will be
                     // used for cutting string for appending info to PYTHIA
  vector<vector<int>> IMStructure2;  // this will contain Tempjunction Indice
                                     // and -+1 and leg numbers
  // ex. 2nd lef in 3rd temp anti junction is same as 1st leg 2nd tempjunction.
  // 3,-1,0,0   3,1,2,0 would be saved in the vectorIMS1
  vector<int> IMStructure3;  // 4 element vector of Tempjunction order, kind and
                             // leg address in JuncStructure vector
  vector<vector<HHparton>>
      Recombearly1;  // these partons are in the junction with three or two
                     // shared legs with others , and they will be recombined
                     // into baryon to cut? or arrange the string for PYTHIA to
                     // understand input come from this code
  vector<HHparton> Recombearly2;
  vector<vector<vector<HHparton>>>
      Dijunction1;  // these partons are in the junction with two legs shared
                    // legs!
  vector<vector<HHparton>> Dijunction2;
  // Vector of Indices for Dijunction1 vector to be ordered for invoking Pythia
  vector<vector<int>> DijunctionInfo1;
  vector<int> DijunctionInfo2;  //  { -1: antiJ, +1 : J , 0 : shared leg}
  // these tags should be assigned to corresponding legs in Dijunction1 vector!
  // so DijunctionInfo2 has five elements{since, dijunction structure has five
  // legs!}

  vector<vector<vector<HHparton>>>
      Singlejunction1;  // easiest case! just single string
  vector<vector<HHparton>> Singlejunction2;
  vector<vector<HHparton>>
      Tailoredstring1;  // when there are the junction with two shared legs, it
                        // also should be recombined into baryon and one string
                        // will remain after, which would be saved in this
                        // vector
  vector<HHparton> Tailoredstring2;
  vector<int>
      realjuncindice;  // vector to save the indice of Tempjunction when real
                       // junction is formed by three initiating particle

  vector<HHparton> finalstring;  // final space for all of the remnant particles
                                 // being corrected by fake parton addition

  /*JSINFO << "SP_remnants before stringprep:";
  for(int irem=0; irem < SP_remnants.num(); ++irem) {
    std::cout << SP_remnants[irem].id() << "," << SP_remnants[irem].col() << ","
  << SP_remnants[irem].acol() << std::endl;
  }*/

  // Tempjunctions missing partons? Add thermal or fake
  // find the maximum (anti-)color tag in the remnants list and the
  // tempjunctions
  int maxtag = 0;
  for (int irem = 0; irem < SP_remnants.num(); ++irem) {
    if (SP_remnants[irem].col() > maxtag) {
      maxtag = SP_remnants[irem].col();
    }
    if (SP_remnants[irem].acol() > maxtag) {
      maxtag = SP_remnants[irem].acol();
    }
  }
  for (int ijunc = 0; ijunc < Tempjunctions.size(); ++ijunc) {
    if (Tempjunctions.at(ijunc).at(1).at(1) > maxtag) {
      maxtag = Tempjunctions.at(ijunc).at(1).at(1);
    }
    if (Tempjunctions.at(ijunc).at(2).at(1) > maxtag) {
      maxtag = Tempjunctions.at(ijunc).at(2).at(1);
    }
    if (Tempjunctions.at(ijunc).at(3).at(1) > maxtag) {
      maxtag = Tempjunctions.at(ijunc).at(3).at(1);
    }
  }

  HHparton tempparton1;
  HHparton tempparton2;
  HHparton tempparton3;
  for (int ijunc = 0; ijunc < Tempjunctions.size(); ++ijunc) {
    int i1 = 0;
    int i2 = 0;
    int i3 = 0;

    // First try to match to remant parton color/anticolor tags and write
    // matching parton index into tag
    for (int irem = 0; irem < SP_remnants.num(); ++irem) {
      if (Tempjunctions.at(ijunc).at(0).at(0) == -1) {
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(1).at(1) &&
            Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          i1 = irem + 1;
        }
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(2).at(1) &&
            Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          i2 = irem + 1;
        }
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(3).at(1) &&
            Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          i3 = irem + 1;
        }
      } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1) {
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(1).at(1) &&
            Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          i1 = irem + 1;
        }
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(2).at(1) &&
            Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          i2 = irem + 1;
        }
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(3).at(1) &&
            Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          i3 = irem + 1;
        }
      }
    }

    // Now try to find partners in another temp junction, store junction and leg
    // number as a negative tag
    for (int ijunc2 = ijunc + 1; ijunc2 < Tempjunctions.size(); ++ijunc2) {
      if ((Tempjunctions.at(ijunc).at(0).at(0) == -1 &&
           Tempjunctions.at(ijunc2).at(0).at(0) == 1) ||
          (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
           Tempjunctions.at(ijunc2).at(0).at(0) == -1)) {
        for (int ileg = 1; ileg <= 3; ++ileg) {
          if (Tempjunctions.at(ijunc2).at(ileg).at(1) ==
                  Tempjunctions.at(ijunc).at(1).at(1) &&
              Tempjunctions.at(ijunc).at(1).at(1) != 0) {
            i1 = -ijunc2 * 10 - ileg;
          }
          if (Tempjunctions.at(ijunc2).at(ileg).at(1) ==
                  Tempjunctions.at(ijunc).at(2).at(1) &&
              Tempjunctions.at(ijunc).at(2).at(1) != 0) {
            i2 = -ijunc2 * 10 - ileg;
          }
          if (Tempjunctions.at(ijunc2).at(ileg).at(1) ==
                  Tempjunctions.at(ijunc).at(3).at(1) &&
              Tempjunctions.at(ijunc).at(3).at(1) != 0) {
            i3 = -ijunc2 * 10 - ileg;
          }
        }
      } else {
        bool warning = false;
        for (int ileg = 1; ileg <= 3; ++ileg) {
          for (int jleg = 1; jleg <= 3; ++jleg) {
            if (Tempjunctions.at(ijunc).at(ileg).at(1) ==
                Tempjunctions.at(ijunc2).at(jleg).at(1)) {
              warning = true;
            }
          }
        }
        if (warning) {
          JSWARN << "There is a junction pair which is not "
                    "junction-antijunction, but junction-junction or "
                    "antijunction-antijunction. This should not happen!";
        }
      }
    }

    // define temporary parton to be used in findclosereplica
    if (i1 <= 0) {
      if (i2 > 0) {
        tempparton1 = SP_remnants[i2 - 1];
        if (std::abs(tempparton1.id()) <
            6) {  // no change of colors here, they are not important in
                  // findclosereplica
          tempparton1.id(-tempparton1.id());
        }
      } else if (i3 > 0) {
        tempparton1 = SP_remnants[i3 - 1];
        if (std::abs(tempparton1.id()) < 6) {
          tempparton1.id(-tempparton1.id());
        }
      } else {  // should be very rare: junction with 3 ghost legs
        tempparton1 = SP_remnants[0];
        if (Tempjunctions.at(ijunc).at(0).at(0) == -1 &&
            std::abs(tempparton1.id()) < 6 && tempparton1.id() < 0) {
          tempparton1.id(-tempparton1.id());
        } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                   std::abs(tempparton1.id()) < 6 && tempparton1.id() > 0) {
          tempparton1.id(-tempparton1.id());
        }
      }
    }
    if (i2 <= 0) {
      if (i1 > 0) {
        tempparton2 = SP_remnants[i1 - 1];
        if (std::abs(tempparton2.id()) <
            6) {  // no change of colors here, they are not important in
                  // findclosereplica
          tempparton2.id(-tempparton2.id());
        }
      } else if (i3 > 0) {
        tempparton2 = SP_remnants[i3 - 1];
        if (std::abs(tempparton2.id()) < 6) {
          tempparton2.id(-tempparton2.id());
        }
      } else {  // should be very rare: junction with 3 ghost legs
        tempparton2 = SP_remnants[0];
        if (Tempjunctions.at(ijunc).at(0).at(0) == -1 &&
            std::abs(tempparton2.id()) < 6 && tempparton2.id() < 0) {
          tempparton2.id(-tempparton2.id());
        } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                   std::abs(tempparton2.id()) < 6 && tempparton2.id() > 0) {
          tempparton2.id(-tempparton2.id());
        }
      }
    }
    if (i3 <= 0) {
      if (i1 > 0) {
        tempparton3 = SP_remnants[i1 - 1];
        if (std::abs(tempparton3.id()) <
            6) {  // no change of colors here, they are not important in
                  // findclosereplica
          tempparton3.id(-tempparton3.id());
        }
      } else if (i2 > 0) {
        tempparton3 = SP_remnants[i2 - 1];
        if (std::abs(tempparton3.id()) < 6) {
          tempparton3.id(-tempparton3.id());
        }
      } else {  // should be very rare: junction with 3 ghost legs
        tempparton3 = SP_remnants[0];
        if (Tempjunctions.at(ijunc).at(0).at(0) == -1 &&
            std::abs(tempparton3.id()) < 6 && tempparton3.id() < 0) {
          tempparton3.id(-tempparton3.id());
        } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                   std::abs(tempparton3.id()) < 6 && tempparton3.id() > 0) {
          tempparton3.id(-tempparton3.id());
        }
      }
    }

    if (i1 == 0 || i2 == 0 || i3 == 0) {
      Tempjunctions.at(ijunc).at(0).at(1) =
          1;  // at least one thermal parton in (anti-)junction
    }

    // Now find partons to add to missing junction legs; either thermal or fake
    if (i1 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int loc = 0;
      if (tempparton1.id() == 21) {
        tempparton1.id(1);
        loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton1.id(21);
      } else {
        loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3);
        HHparton fakep = tempparton1;
        fakep.id(fid);
        fakep.set_color(0);
        if (Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          fakep.set_anti_color(Tempjunctions.at(ijunc).at(1).at(1));
        } else {
          fakep.set_anti_color(++maxtag);
          Tempjunctions.at(ijunc).at(1).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          HH_thermal[-loc - 1].acol(Tempjunctions.at(ijunc).at(1).at(1));
        } else {
          HH_thermal[-loc - 1].acol(++maxtag);
          Tempjunctions.at(ijunc).at(1).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].col(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    if (i2 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int loc = 0;
      if (tempparton2.id() == 21) {
        tempparton2.id(1);
        loc = findcloserepl(tempparton2, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton2.id(21);
      } else {
        loc = findcloserepl(tempparton2, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3);
        HHparton fakep = tempparton2;
        fakep.id(fid);
        fakep.set_color(0);
        if (Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          fakep.set_anti_color(Tempjunctions.at(ijunc).at(2).at(1));
        } else {
          fakep.set_anti_color(++maxtag);
          Tempjunctions.at(ijunc).at(2).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          HH_thermal[-loc - 1].acol(Tempjunctions.at(ijunc).at(2).at(1));
        } else {
          HH_thermal[-loc - 1].acol(++maxtag);
          Tempjunctions.at(ijunc).at(2).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].col(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    if (i3 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int loc = 0;
      if (tempparton3.id() == 21) {
        tempparton3.id(1);
        loc = findcloserepl(tempparton3, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton3.id(21);
      } else {
        loc = findcloserepl(tempparton3, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3);
        HHparton fakep = tempparton3;
        fakep.id(fid);
        fakep.set_color(0);
        if (Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          fakep.set_anti_color(Tempjunctions.at(ijunc).at(3).at(1));
        } else {
          fakep.set_anti_color(++maxtag);
          Tempjunctions.at(ijunc).at(3).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          HH_thermal[-loc - 1].acol(Tempjunctions.at(ijunc).at(3).at(1));
        } else {
          HH_thermal[-loc - 1].acol(++maxtag);
          Tempjunctions.at(ijunc).at(3).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].col(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    if (i1 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int loc = 0;
      if (tempparton1.id() == 21) {
        tempparton1.id(-1);
        loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton1.id(21);
      } else {
        loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3);
        HHparton fakep = tempparton1;
        fakep.id(fid);
        fakep.set_anti_color(0);
        if (Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          fakep.set_color(Tempjunctions.at(ijunc).at(1).at(1));
        } else {
          fakep.set_color(++maxtag);
          Tempjunctions.at(ijunc).at(1).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(1).at(1) != 0) {
          HH_thermal[-loc - 1].col(Tempjunctions.at(ijunc).at(1).at(1));
        } else {
          HH_thermal[-loc - 1].col(++maxtag);
          Tempjunctions.at(ijunc).at(1).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].acol(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    if (i2 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int loc = 0;
      if (tempparton2.id() == 21) {
        tempparton2.id(-1);
        loc = findcloserepl(tempparton2, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton2.id(21);
      } else {
        loc = findcloserepl(tempparton2, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3);
        HHparton fakep = tempparton2;
        fakep.id(fid);
        fakep.set_anti_color(0);
        if (Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          fakep.set_color(Tempjunctions.at(ijunc).at(2).at(1));
        } else {
          fakep.set_color(++maxtag);
          Tempjunctions.at(ijunc).at(2).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(2).at(1) != 0) {
          HH_thermal[-loc - 1].col(Tempjunctions.at(ijunc).at(2).at(1));
        } else {
          HH_thermal[-loc - 1].col(++maxtag);
          Tempjunctions.at(ijunc).at(2).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].acol(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    if (i3 == 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int loc = 0;
      if (tempparton3.id() == 21) {
        tempparton3.id(-1);
        loc = findcloserepl(tempparton3, 1, false, true, HH_showerptns,
                            HH_thermal);
        tempparton3.id(21);
      } else {
        loc = findcloserepl(tempparton3, 1, false, true, HH_showerptns,
                            HH_thermal);
      }
      if (loc == 999999999 || loc > 0) {
        int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3);
        HHparton fakep = tempparton3;
        fakep.id(fid);
        fakep.set_anti_color(0);
        if (Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          fakep.set_color(Tempjunctions.at(ijunc).at(3).at(1));
        } else {
          fakep.set_color(++maxtag);
          Tempjunctions.at(ijunc).at(3).at(1) = maxtag;
        }
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        if (Tempjunctions.at(ijunc).at(3).at(1) != 0) {
          HH_thermal[-loc - 1].col(Tempjunctions.at(ijunc).at(3).at(1));
        } else {
          HH_thermal[-loc - 1].col(++maxtag);
          Tempjunctions.at(ijunc).at(3).at(1) = maxtag;
        }
        HH_thermal[-loc - 1].acol(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }

    // Add fake gluons (no beam energies) between junctions with empty
    // connecting legs [handling could be improved later]
    if (i1 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int jleg = ((-i1) % 10);
      int jjunc = ((-i1) - jleg) / 10;
      HHparton fakep = tempparton1;
      fakep.id(21);
      fakep.set_color(++maxtag);
      fakep.set_anti_color(Tempjunctions.at(ijunc).at(1).at(1));
      fakep.mass(1e-2);  // 10 MeV "accuracy"
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }

    if (i2 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int jleg = ((-i2) % 10);
      int jjunc = ((-i2) - jleg) / 10;
      HHparton fakep = tempparton2;
      fakep.id(21);
      fakep.set_color(++maxtag);
      fakep.set_anti_color(Tempjunctions.at(ijunc).at(2).at(1));
      fakep.mass(1e-2);
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }

    if (i3 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == -1) {
      int jleg = ((-i3) % 10);
      int jjunc = ((-i3) - jleg) / 10;
      HHparton fakep = tempparton3;
      fakep.id(21);
      fakep.set_color(++maxtag);
      fakep.set_anti_color(Tempjunctions.at(ijunc).at(3).at(1));
      fakep.mass(1e-2);
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }

    if (i1 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int jleg = ((-i1) % 10);
      int jjunc = ((-i1) - jleg) / 10;
      HHparton fakep = tempparton1;
      fakep.id(21);
      fakep.set_anti_color(++maxtag);
      fakep.set_color(Tempjunctions.at(ijunc).at(1).at(1));
      fakep.mass(1e-2);
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }

    if (i2 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int jleg = ((-i2) % 10);
      int jjunc = ((-i2) - jleg) / 10;
      HHparton fakep = tempparton2;
      fakep.id(21);
      fakep.set_anti_color(++maxtag);
      fakep.set_color(Tempjunctions.at(ijunc).at(2).at(1));
      fakep.mass(1e-2);
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }

    if (i3 < 0 && Tempjunctions.at(ijunc).at(0).at(0) == 1) {
      int jleg = ((-i3) % 10);
      int jjunc = ((-i3) - jleg) / 10;
      HHparton fakep = tempparton3;
      fakep.id(21);
      fakep.set_anti_color(++maxtag);
      fakep.set_color(Tempjunctions.at(ijunc).at(3).at(1));
      fakep.mass(1e-2);
      Tempjunctions.at(jjunc).at(jleg).at(1) = maxtag;
      double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
      double fake_phi = 2. * 3.14159265358979 * ran();
      fakep.px(fake_pT * cos(fake_phi));
      fakep.py(fake_pT * sin(fake_phi));
      fakep.pz(0.);
      fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                        fakep.pz() * fakep.pz() + fakep.mass() * fakep.mass()));
      fakep.orig(-1);
      fakep.is_remnant(true);
      fakep.is_fakeparton(true);
      SP_remnants.add(fakep);
    }
  }

  // repair all strings (add missing partons)
  for (int ptn1 = 0; ptn1 < SP_remnants.num(); ptn1++) {
    bool colors_match = false;
    bool anti_colors_match = false;
    int ptn1_col = SP_remnants[ptn1].col();
    int ptn1_acol = SP_remnants[ptn1].acol();
    for (int ptn2 = 0; ptn2 < SP_remnants.num(); ptn2++) {
      int ptn2_acol = SP_remnants[ptn2].acol();
      int ptn2_col = SP_remnants[ptn2].col();
      if (ptn1_col == ptn2_acol) {
        colors_match = true;
      }
      if (ptn1_acol == ptn2_col) {
        anti_colors_match = true;
      }
    }
    if (!colors_match || !anti_colors_match) {
      for (int ijunc = 0; ijunc < Tempjunctions.size(); ijunc++) {
        if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
            (ptn1_col == Tempjunctions.at(ijunc).at(1).at(1) ||
             ptn1_col == Tempjunctions.at(ijunc).at(2).at(1) ||
             ptn1_col == Tempjunctions.at(ijunc).at(3).at(1))) {
          colors_match = true;
        }
        if (Tempjunctions.at(ijunc).at(0).at(0) == -1 &&
            (ptn1_acol == Tempjunctions.at(ijunc).at(1).at(1) ||
             ptn1_acol == Tempjunctions.at(ijunc).at(2).at(1) ||
             ptn1_acol == Tempjunctions.at(ijunc).at(3).at(1))) {
          anti_colors_match = true;
        }
      }
    }
    if (!colors_match && SP_remnants[ptn1].id() > 0) {
      int loc = 0;
      if (SP_remnants[ptn1].id() == 21) {
        SP_remnants[ptn1].id(1);
        loc = findcloserepl(SP_remnants[ptn1], ptn1 + 1, false, true,
                            HH_showerptns, HH_thermal);
        SP_remnants[ptn1].id(21);
      } else {
        loc = findcloserepl(SP_remnants[ptn1], ptn1 + 1, false, true,
                            HH_showerptns, HH_thermal);
      }

      if (loc == 999999999 || loc > 0) {
        int fid = (ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3;
        HHparton fakep = SP_remnants[ptn1];
        fakep.id(fid);
        fakep.set_color(0);
        fakep.set_anti_color(ptn1_col);
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double dir = (ran() < 0.5) ? 1. : -1.;
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        HH_thermal[-loc - 1].acol(ptn1_col);
        HH_thermal[-loc - 1].col(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }
    if (!anti_colors_match &&
        (SP_remnants[ptn1].id() < 0 || SP_remnants[ptn1].id() == 21)) {
      int loc = 0;
      if (SP_remnants[ptn1].id() == 21) {
        SP_remnants[ptn1].id(-1);
        loc = findcloserepl(SP_remnants[ptn1], ptn1 + 1, false, true,
                            HH_showerptns, HH_thermal);
        SP_remnants[ptn1].id(21);
      } else {
        loc = findcloserepl(SP_remnants[ptn1], ptn1 + 1, false, true,
                            HH_showerptns, HH_thermal);
      }

      if (loc == 999999999 || loc > 0) {
        int fid = (ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3;
        HHparton fakep = SP_remnants[ptn1];
        fakep.id(fid);
        fakep.set_color(ptn1_acol);
        fakep.set_anti_color(0);
        if (std::abs(fakep.id()) < 3) {
          fakep.mass(xmq);
        } else {
          fakep.mass(xms);
        }
        double dir = (ran() < 0.5) ? 1. : -1.;
        double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
        double fake_phi = 2. * 3.14159265358979 * ran();
        fakep.px(fake_pT * cos(fake_phi));
        fakep.py(fake_pT * sin(fake_phi));
        if (number_p_fake == 0 || number_p_fake == 2) {
          fakep.pz(-p_fake);
        } else if (number_p_fake == 1 || number_p_fake == 3) {
          fakep.pz(p_fake);
        } else {
          fakep.pz(0.);
        }
        number_p_fake++;
        fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                          fakep.pz() * fakep.pz() +
                          fakep.mass() * fakep.mass()));
        fakep.orig(-1);
        fakep.is_remnant(true);
        fakep.is_fakeparton(true);
        SP_remnants.add(fakep);
      } else if (loc < 0) {
        HH_thermal[-loc - 1].col(ptn1_acol);
        HH_thermal[-loc - 1].acol(0);
        HH_thermal[-loc - 1].is_used(true);
        HH_thermal[-loc - 1].is_remnant(true);
        SP_remnants.add(HH_thermal[-loc - 1]);
        SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
      }
    }
  }

  set_initial_parton_masses(SP_remnants);

  for (int ijunc = 0; ijunc < Tempjunctions.size();
       ++ijunc) {  // Let's make Legs for Final string by verifying
                   // Tempjunctions! starting from the kind of junction(-1 or
                   // +1)
    // std::cout <<endl<<endl<<"  Let's see candidates for (anti)junction
    // "<<endl; std::cout <<" Candidate " << ijunc+1 <<endl; std::cout <<" [ (
    // " << Tempjunctions.at(ijunc).at(0).at(0) <<" , " <<
    // Tempjunctions.at(ijunc).at(0).at(1) <<" ), " ; std::cout <<"  (  " <<
    // Tempjunctions.at(ijunc).at(1).at(0) <<" , " <<
    // Tempjunctions.at(ijunc).at(1).at(1) <<" ), " ; std::cout <<"  (  " <<
    // Tempjunctions.at(ijunc).at(2).at(0) <<" , " <<
    // Tempjunctions.at(ijunc).at(2).at(1) <<" ), " ; std::cout <<"  (  " <<
    // Tempjunctions.at(ijunc).at(3).at(0) <<" , " <<
    // Tempjunctions.at(ijunc).at(3).at(1) <<" ), " ; std::cout <<" ] "<<endl;

    // std::cout <<endl<<" I'm Working !! from line 2027" <<endl;
    vector<int> correction;  // for the case of 1 or 2 initiating particles, we
                             // need to remember the partons and get the
                             // used_junction tag to the original
    if (Tempjunctions.at(ijunc).at(0).at(0) ==
        -1) {  // check anti-color tags in SP_remnants partons!! to form
               // anti-junction
      for (int irem = 0; irem < SP_remnants.num();
           ++irem) {  // searching through all remnant particles
        // std::cout <<endl<<" Let's see " << irem+1 <<"th remnant particle with
        // "<< "pid : " << SP_remnants[irem].pid() <<" color : "<<
        // SP_remnants[irem].col()<< " anti-color : " <<
        // SP_remnants[irem].acol() <<endl;
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(1).at(1)) {
          Leg1.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          // std::cout <<endl<<" Leg 1 has initiating particle " <<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].acol()<<endl<<endl;
          Legconsidering.push_back(irem);
        }
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(2).at(1)) {
          Leg2.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          // std::cout <<endl<<" Leg 2 has initiating particle " <<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].acol()<<endl<<endl;
          Legconsidering.push_back(irem);
        }
        if (SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(3).at(1)) {
          Leg3.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          // std::cout <<endl<<" Leg 3 has initiating particle " <<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].acol()<<endl<<endl;
          Legconsidering.push_back(irem);
        }
        if (Legconsidering.size() !=
            3) {  // if only one or two particles found in remnants, we need to
                  // set the tag to the original
          for (int icor = 0; icor < correction.size(); icor++) {
            SP_remnants[correction.at(icor)].used_junction(false);
            // dignostic measure
            // std::cout <<endl<<"the "<<icor<<" th particle with color tag of
            // "<<SP_remnants[icor].col()<<" , "<<SP_remnants[icor].acol()<<"
            // are turned into unused particle"<<endl;
          }
          correction.clear();
        }
      }
    } else if (Tempjunctions.at(ijunc).at(0).at(0) ==
               1) {  // check color tags in SP_remnants partons!! to form
                     // junction
      for (int irem = 0; irem < SP_remnants.num();
           ++irem) {  // searching through all remnant particles
        // std::cout <<endl<<" Let's see " << irem+1 <<"th remnant particle with
        // "<< "pid : " << SP_remnants[irem].pid() <<" color : "<<
        // SP_remnants[irem].col()<< " anti-color : " <<
        // SP_remnants[irem].acol() <<endl;
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(1).at(1)) {
          Leg1.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          Legconsidering.push_back(irem);
          // std::cout <<endl<<" Leg 1 has initiating particle "<<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].col()<<endl<<endl;
        }
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(2).at(1)) {
          Leg2.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          // std::cout <<endl<<" Leg 2 has initiating particle " <<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].col()<<endl<<endl;
          Legconsidering.push_back(irem);
        }
        if (SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(3).at(1)) {
          Leg3.push_back(SP_remnants[irem]);
          correction.push_back(irem);
          SP_remnants[irem].used_junction(true);
          // std::cout <<endl<<" Leg 3 has initiating particle " <<" and the pID
          // is "<< SP_remnants[irem].pid()<<" and the color tag :
          // "<<SP_remnants[irem].col()<<endl<<endl;
          Legconsidering.push_back(irem);
        }
        if (Legconsidering.size() != 3) {
          for (int icor = 0; icor < correction.size(); icor++) {
            SP_remnants[correction.at(icor)].used_junction(false);
            // dignostic measure
            // std::cout <<endl<<"the "<<icor<<" th particle with color tag of
            // "<<SP_remnants[icor].col()<<" , "<<SP_remnants[icor].acol()<<"
            // are turned into unused particle"<<endl;
          }
          correction.clear();
        }
      }
    }

    if (Legconsidering.size() == 3) {
      realjuncindice.push_back(ijunc);
      // std::cout <<endl<<"the indice added for string repair : "<<ijunc<<endl;
      if (Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg1.at(0).col() != 0 &&
          Leg1.at(0).acol() !=
              0) {  // starting to form the Leg1 for anti_junction by Color Tag
                    // Tracing. If first particle is (anti)quark, no need to
                    // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            if ((Leg1.back().col() != 0) &&
                (Leg1.back().col() == SP_remnants[icf].acol())) {
              // std::cout <<endl<<"particle added to Leg 1 with color tag  (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg1.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg1.back().col() == 0 || Leg1.back().acol() == 0) {
              break;
            }
          }
        }
      } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                 Leg1.at(0).col() != 0 &&
                 Leg1.at(0).acol() !=
                     0) {  // starting to form the Leg1 junction by Color Tag
                           // Tracing. If first particle is quark, no need to
                           // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            if ((Leg1.back().acol() != 0) &&
                (Leg1.back().acol() == SP_remnants[icf].col())) {
              // std::cout <<endl<<"particle added to Leg 1 with color tag (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg1.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg1.back().col() == 0 || Leg1.back().acol() == 0) {
              break;
            }
          }
        }
      }

      if (Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg2.at(0).col() != 0 &&
          Leg2.at(0).acol() !=
              0) {  // starting to form the Leg2 for anti_junction by Color Tag
                    // Tracing. If first particle is (anti)quark, no need to
                    // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            if ((Leg2.back().col() != 0) &&
                (Leg2.back().col() == SP_remnants[icf].acol())) {
              // std::cout <<endl<<"particle added to Leg 2 with color tag (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg2.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg2.back().col() == 0 || Leg2.back().acol() == 0) {
              break;
            }
          }
        }
      } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                 Leg2.at(0).col() != 0 &&
                 Leg2.at(0).acol() !=
                     0) {  // starting to form the Leg2 junction by Color Tag
                           // Tracing. If first particle is quark, no need to
                           // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            if ((Leg2.back().acol() != 0) &&
                (Leg2.back().acol() == SP_remnants[icf].col())) {
              // std::cout <<endl<<"particle added to Leg 2 with color tag (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg2.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg2.back().col() == 0 || Leg2.back().acol() == 0) {
              break;
            }
          }
        }
      }

      if (Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg3.at(0).col() != 0 &&
          Leg3.at(0).acol() !=
              0) {  // starting to form the Leg3 for anti_junction by Color Tag
                    // Tracing. If first particle is (anti)quark, no need to
                    // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            if ((Leg3.back().col() != 0) &&
                (Leg3.back().col() == SP_remnants[icf].acol())) {
              // std::cout <<endl<<"particle added to Leg 3 with color tag (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg3.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg3.back().col() == 0 || Leg3.back().acol() == 0) {
              break;
            }
          }
        }
      } else if (Tempjunctions.at(ijunc).at(0).at(0) == 1 &&
                 Leg3.at(0).col() != 0 &&
                 Leg3.at(0).acol() !=
                     0) {  // starting to form the Leg3 junction by Color Tag
                           // Tracing. If first particle is quark, no need to
                           // trace
        for (int iloop = 0; iloop < SP_remnants.num(); iloop++) {
          for (int icf = 0; icf < SP_remnants.num(); ++icf) {
            // std::cout <<"candidate particle is "<<SP_remnants[icf].pid()<<" ,
            // "<<SP_remnants[icf].col()<<" , "<<SP_remnants[icf].acol()<<" ,
            // "<<SP_remnants[icf].used_junction()<<endl; std::cout <<"and
            // standard particle is "<<Leg3.back().pid()<<" ,
            // "<<Leg3.back().col()<<" , "<<Leg3.back().acol()<<" ,
            // "<<Leg3.back().used_junction()<<endl;
            if ((Leg3.back().acol() != 0) &&
                (Leg3.back().acol() == SP_remnants[icf].col())) {
              // std::cout <<endl<<"particle added to Leg 3 with color tag (
              // "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<"
              // ) "<<endl;
              Leg3.push_back(SP_remnants[icf]);
              SP_remnants[icf].used_junction(true);
            }
            if (Leg3.back().col() == 0 || Leg3.back().acol() == 0) {
              break;
            }
          }
        }
      }
      JuncLegs.push_back(Leg1);
      JuncLegs.push_back(Leg2);
      JuncLegs.push_back(Leg3);
      JuncStructure.push_back(
          JuncLegs);  // now junction structure with three complete legs is
                      // saved, so clear previous infos for considering next
                      // junction candidate.
    }
    JuncLegs.clear();
    Leg1.clear();
    Leg2.clear();
    Leg3.clear();
    Legconsidering.clear();
  }

  // So far, all the information for (Anti)Junction and Legs are formed by color
  // tag tracing. Let's Check this out dignostic measure
  /*for(int ijuncstr=0; ijuncstr < JuncStructure.size(); ++ijuncstr) {
    std::cout <<endl<<endl<<"  Let's see Temporary junction structure   "<<endl;
    std::cout <<" Temp Junction :" << ijuncstr+1 <<endl;
    std::cout <<"Leg 1 : [ ";
    for(int ileg1=0; ileg1 < JuncStructure.at(ijuncstr).at(0).size(); ++ileg1){
      std::cout <<" (  " << JuncStructure.at(ijuncstr).at(0)[ileg1].col() <<" ,
  " << JuncStructure.at(ijuncstr).at(0)[ileg1].acol() <<" ), " ;
    }
    std::cout <<" ] " <<endl;
    std::cout <<"Leg 2 : [ ";
    for(int ileg2=0; ileg2 < JuncStructure.at(ijuncstr).at(1).size(); ++ileg2){
      std::cout <<" (  " << JuncStructure.at(ijuncstr).at(1)[ileg2].col() <<" ,
  " << JuncStructure.at(ijuncstr).at(1)[ileg2].acol() <<" ), " ;
    }
    std::cout <<" ] " <<endl;
    std::cout <<"Leg 3 : [ ";
    for(int ileg3=0; ileg3 < JuncStructure.at(ijuncstr).at(2).size(); ++ileg3){
      std::cout <<" (  " << JuncStructure.at(ijuncstr).at(2)[ileg3].col() <<" ,
  " << JuncStructure.at(ijuncstr).at(2)[ileg3].acol() <<" ), " ;
    }
    std::cout <<" ] " <<endl;
  }*/

  // Now, we will find the shared legs and sorting them into the corresponding
  // vectors, the information from JuncStructure vector and realjuncindice
  // vector work together here. if the tags of first particle in a string and
  // end particle in another string are same, they are shared Leg!!
  for (int irep1 = 0; irep1 < JuncStructure.size(); irep1++) {
    IMStructure3.push_back(realjuncindice.at(irep1));
    IMStructure3.push_back(
        Tempjunctions.at(realjuncindice.at(irep1)).at(0).at(0));
    IMStructure3.push_back(irep1);  // room for other usage
    IMStructure3.push_back(0);      // tag 1: dijinction, tag0: other cases
    IMStructure2.push_back(IMStructure3);
    IMStructure3.clear();
    for (int irep2 = 0; irep2 < 3; irep2++) {
      for (int irep3 = 0; irep3 < JuncStructure.size(); irep3++) {
        for (int irep4 = 0; irep4 < 3; irep4++) {
          if ((irep1 != irep3) &&
              (JuncStructure.at(irep1).at(irep2).at(0).col() ==
               JuncStructure.at(irep3).at(irep4).back().col()) &&
              (JuncStructure.at(irep1).at(irep2).at(0).acol() ==
               JuncStructure.at(irep3).at(irep4).back().acol())) {
            IMStructure3.push_back(irep1);
            IMStructure3.push_back(irep2);
            IMStructure3.push_back(irep3);
            IMStructure3.push_back(
                irep4);  // the first pair means irep2_th Leg in irep1_th
                         // junction is linked with irep4_th Leg in irep3_th
                         // Tempjunction
            IMStructure2.push_back(IMStructure3);
            IMStructure3.clear();
          }
        }
      }
    }
    IMStructure1.push_back(IMStructure2);
    IMStructure2.clear();
  }

  // dignostic measure!
  /*JSINFO << "First IMStructure printout";
  std::cout <<endl;
  for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
    for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
      std::cout <<" ( ";
      for(int icheck3 = 0; icheck3 <
  IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout
  <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
      }
      std::cout <<" ) ";
    }
    std::cout <<endl;
  }*/

  // create a vector containing 0 if the (anti-)junction contains only shower
  // partons, 1 if at least one thermal this is only for the ones in
  // Recombearly1
  std::vector<int> thermal_parton_junction;

  // since the relation between junctions is defined, Let's sorting Junctions
  // based on shared legs
  for (int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++) {
    // structure for 3-shared particles find shared legs and tossing them into
    // the baryon formation
    if (IMStructure1.at(iloop1).size() == 4) {
      for (int iloop2 = 1; iloop2 < IMStructure1.at(iloop1).size(); iloop2++) {
        for (int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++) {
          for (int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size();
               iloop4++) {
            if ((IMStructure1.at(iloop1).at(iloop2).at(2) ==
                 IMStructure1.at(iloop3).at(iloop4).at(0)) &&
                (IMStructure1.at(iloop1).at(iloop2).at(3) ==
                 IMStructure1.at(iloop3).at(iloop4).at(1))) {
              vector<int> testing;
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(
                  1));  // now temporary vector is created by this procedure.
              std::vector<vector<int>>::iterator it =
                  std::find(IMStructure1.at(iloop3).begin(),
                            IMStructure1.at(iloop3).end(), testing);
              IMStructure1.at(iloop3).erase(it);
              testing.clear();  // now the one of the leg is eliminated in the
                                // list, so there is no overlapped leg to be
                                // tossed into baryon formation list.
            }
            if ((JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                     .at(0)
                     .col() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                     .at(0)
                     .acol() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                     .back()
                     .col() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                     .back()
                     .acol() != 0)) {
              parton_collection tempqpair;
              gluon_decay(
                  JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                      .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                      .at(0),
                  tempqpair);  // the first particle has col tag and the second
                               // had acol tag
              // std::vector<HHparton>::iterator tempit1 =
              // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).begin();
              // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).erase(tempit1);
              if (IMStructure1.at(iloop1).at(0).at(1) == -1) {
                Recombearly2.push_back(
                    tempqpair[1]);  // add anti-particle to form anti baryon
                // std::cout <<endl<<" the particle with tags { "<<
                // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().col()<<"
                // , "
                //<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().acol()<<"
                //} is changed into ";
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .pop_back();
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .push_back(
                        tempqpair[0]);  // add remain parton to the original leg
                // std::cout <<" { " << tempqpair[0].col()<<" , "<<
                // tempqpair[0].acol()<<" }"<<endl; std::cout <<" and "<<" { "
                // << tempqpair[1].col()<<" , "<< tempqpair[1].acol()<<" } "<<"
                // is gone to be recombined"<<endl;
              }
              if (IMStructure1.at(iloop1).at(0).at(1) == 1) {
                Recombearly2.push_back(
                    tempqpair[0]);  // add particle to form baryon
                // std::cout <<endl<<" the particle with tags { "<<
                // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().col()<<"
                // , "
                //<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().acol()<<"
                //} is changed into ";
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .pop_back();
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .push_back(
                        tempqpair[1]);  // add remain parton to the original leg
                // std::cout <<" { " << tempqpair[0].col()<<" , "<<
                // tempqpair[0].acol()<<" }"<<endl; std::cout <<" and "<<" { "
                // << tempqpair[1].col()<<" , "<< tempqpair[1].acol()<<" } "<<"
                // is gone to be recombined"<<endl;
              }
              tempqpair.clear();
            }
          }
        }
      }
      thermal_parton_junction.push_back(
          Tempjunctions.at(IMStructure1.at(iloop1).at(0).at(2)).at(0).at(1));
      IMStructure1.erase(IMStructure1.begin() +
                         iloop1);  // erase the junction with three shared legs
                                   // from the IMStructure
      iloop1--;
      Recombearly1.push_back(Recombearly2);
      Recombearly2.clear();
    }
  }

  // One thing to remind : JuncStructure has information of parton class, and
  // IMStructure has informations about these partons or junctions with same
  // indice with JuncStructure
  for (int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++) {
    if (IMStructure1.at(iloop1).size() ==
        3) {  // 2-shared leg case {basic info, shared leg1's info, shared
              // leg2's info}
      // loop for pop out the unshared leg in the junction structure.
      int leftleg;  // indice of unshared leg which will remain after the Early
                    // recombinaton
      int SI = 0;   // Sum of the Indice of Leg{if 1 = 0+1, Leg3 is not shared,
                    // if 3 = 1+2, Leg 0 is not shared}
      for (int itail = 1; itail < IMStructure1.at(iloop1).size(); itail++) {
        SI = SI + IMStructure1.at(iloop1).at(itail).at(1);
      }
      // std::cout <<endl<<"SI = "<<SI<<endl;
      vector<int> indicevec;  // vector for indice, we could use if statement,
                              // but when if and for statement are repeated,
                              // possibly forced execution occurs regardless of
                              // if statement, so take most defensive measure.
      indicevec.push_back(2);
      indicevec.push_back(1);
      indicevec.push_back(0);
      leftleg = indicevec[SI - 1];
      // below was tested, but showed forced execution{that is all of the Three
      // statements executed so that left leg was always 0} if(SI = 1 ){ leftleg
      // = 2; SI = 0;} if(SI = 2 ){ leftleg = 1; SI = 0;} if(SI = 3 ){ leftleg =
      // 0; SI = 0;} std::cout <<endl<<"leftleg number is "<<leftleg<<endl;

      int jidentity = IMStructure1.at(iloop1).at(0).at(
          1);  // junction identity 1 or -1 {junction or anti juncttion}
      // std::cout <<endl<<"junction identity is "<<jidentity<<endl;
      int juncnum = IMStructure1.at(iloop1).at(0).at(
          0);  // junction number to be dealt with
      vector<HHparton> tempstring =
          JuncStructure.at(juncnum).at(leftleg);  // declare the unshared leg.
      if (jidentity == 1 && tempstring.at(0).col() != 0 &&
          tempstring.at(0).acol() != 0) {
        parton_collection tempqpair;
        gluon_decay(tempstring.at(0), tempqpair);
        tempstring.erase(tempstring.begin());
        tempstring.insert(tempstring.begin(), tempqpair[1]);
        Tailoredstring1.push_back(tempstring);
        Recombearly2.push_back(tempqpair[0]);
      } else if (jidentity == -1 && tempstring.at(0).col() != 0 &&
                 tempstring.at(0).acol() != 0) {
        parton_collection tempqpair;
        gluon_decay(tempstring.at(0), tempqpair);
        tempstring.erase(tempstring.begin());
        tempstring.insert(tempstring.begin(), tempqpair[0]);
        Tailoredstring1.push_back(tempstring);
        Recombearly2.push_back(tempqpair[1]);
      } else {
        Recombearly2.push_back(tempstring.at(0));
      }
      indicevec.clear();

      for (int iloop2 = 1; iloop2 < IMStructure1.at(iloop1).size(); iloop2++) {
        for (int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++) {
          for (int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size();
               iloop4++) {
            if ((IMStructure1.at(iloop1).at(iloop2).at(2) ==
                 IMStructure1.at(iloop3).at(iloop4).at(0)) &&
                (IMStructure1.at(iloop1).at(iloop2).at(3) ==
                 IMStructure1.at(iloop3).at(iloop4).at(1))) {
              vector<int> testing;
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(
                  1));  // now temporary vector is created by this procedure.
              std::vector<vector<int>>::iterator it =
                  std::find(IMStructure1.at(iloop3).begin(),
                            IMStructure1.at(iloop3).end(), testing);
              IMStructure1.at(iloop3).erase(it);
              testing.clear();  // now the one of the leg is eliminated in the
                                // list, so there is no overlapped leg to be
                                // tossed into baryon formation list.
            }
            if ((JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                     .at(0)
                     .col() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                     .at(0)
                     .acol() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                     .back()
                     .col() != 0) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                     .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                     .back()
                     .acol() != 0)) {
              parton_collection tempqpair;
              gluon_decay(
                  JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                      .at(IMStructure1.at(iloop1).at(iloop2).at(1))
                      .at(0),
                  tempqpair);  // the first particle has col tag and the second
                               // had acol tag
              // std::vector<HHparton>::iterator tempit1 =
              // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).begin();
              // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).erase(tempit1);
              if (IMStructure1.at(iloop1).at(0).at(1) == -1) {
                Recombearly2.push_back(
                    tempqpair[1]);  // add anti-particle to form anti baryon
                // std::cout <<endl<<" the particle with tags { "<<
                // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().col()<<"
                // , "
                //<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().acol()<<"
                //} is changed into ";
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .pop_back();
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .push_back(
                        tempqpair[0]);  // add remain parton to the original leg
                // std::cout <<" { " << tempqpair[0].col()<<" , "<<
                // tempqpair[0].acol()<<" }"<<endl; std::cout <<" and "<<" { "
                // << tempqpair[1].col()<<" , "<< tempqpair[1].acol()<<" } "<<"
                // is gone to be recombined"<<endl;
              }
              if (IMStructure1.at(iloop1).at(0).at(1) == 1) {
                Recombearly2.push_back(
                    tempqpair[0]);  // add particle to form baryon
                // std::cout <<endl<<" the particle with tags { "<<
                // JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().col()<<"
                // , "
                //<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().acol()<<"
                //} is changed into ";
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .pop_back();
                JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2))
                    .at(IMStructure1.at(iloop1).at(iloop2).at(3))
                    .push_back(
                        tempqpair[1]);  // add remain parton to the original leg
                // std::cout <<" { " << tempqpair[0].col()<<" , "<<
                // tempqpair[0].acol()<<" }"<<endl; std::cout <<" and "<<" { "
                // << tempqpair[1].col()<<" , "<< tempqpair[1].acol()<<" } "<<"
                // is gone to be recombined"<<endl;
              }
            }
          }
        }
      }
      thermal_parton_junction.push_back(
          Tempjunctions.at(IMStructure1.at(iloop1).at(0).at(2)).at(0).at(1));
      IMStructure1.erase(IMStructure1.begin() +
                         iloop1);  // erase the junction with three shared legs
                                   // from the IMStructure
      iloop1--;
      Recombearly1.push_back(Recombearly2);
      Recombearly2.clear();
    }
  }

  /*std::cout <<endl;
  std::cout <<" number of early recombined baryons 1:
  "<<Recombearly1.size()<<endl; for(int ibary1 = 0; ibary1 <
  Recombearly1.size(); ibary1++){ std::cout <<"the "<< ibary1 <<" st baryon is
  made of "; for(int ibary2 = 0; ibary2 < Recombearly1.at(ibary1).size();
  ibary2++){ std::cout <<" { "<<Recombearly1.at(ibary1).at(ibary2).col()<<" ,
  "<< Recombearly1.at(ibary1).at(ibary2).acol()  <<" } , ";
    }
    std::cout <<endl;
  }*/

  // dignostic measure{IMStructure1}
  /*JSINFO << "1.5 IMStructure printout";
  std::cout <<endl;
  for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
    for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
      std::cout <<" ( ";
      for(int icheck3 = 0; icheck3 <
  IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout
  <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
      }
      std::cout <<" ) ";
    }
    std::cout <<endl;
  }*/

  for (int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++) {
    if (IMStructure1.at(iloop1).size() == 2 &&
        IMStructure1.at(iloop1).at(0).at(3) == 0) {  // Dijunction Structure
      for (int iloop2 = 1; iloop2 < IMStructure1.at(iloop1).size(); iloop2++) {
        for (int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++) {
          for (int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size();
               iloop4++) {
            if ((IMStructure1.at(iloop1).at(iloop2).at(2) ==
                 IMStructure1.at(iloop3).at(iloop4).at(0)) &&
                (IMStructure1.at(iloop1).at(iloop2).at(3) ==
                 IMStructure1.at(iloop3).at(iloop4).at(
                     1))) {  // the case where we found the shared legs pair.
              if (IMStructure1.at(iloop1).at(0).at(3) == 0) {
                for (int idijunc1 = 0; idijunc1 < 3;
                     idijunc1++) {  // first, put three legs of first junction
                  Dijunction2.push_back(
                      JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0))
                          .at(idijunc1));
                }
                for (int idijunc2 = 0; idijunc2 < 3;
                     idijunc2++) {  // second, exclude shared leg in second
                                    // junction to be merged
                  if (idijunc2 != IMStructure1.at(iloop1).at(iloop2).at(3)) {
                    Dijunction2.push_back(
                        JuncStructure
                            .at(IMStructure1.at(iloop1).at(iloop2).at(2))
                            .at(idijunc2));
                  }
                }
              }
              Dijunction1.push_back(
                  Dijunction2);  // finally, five legs are added in the vector
                                 // to form dijunction structure.
              Dijunction2.clear();  // clear temp vector for next procedure.
              IMStructure1.at(iloop1).at(0).pop_back();
              IMStructure1.at(iloop1).at(0).push_back(
                  1);  // tag representing dijunction structure.
              IMStructure1.at(iloop3).at(0).pop_back();
              IMStructure1.at(iloop3).at(0).push_back(
                  1);  // tag representing dijunction structure.

              vector<int> testing;
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
              testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(
                  1));  // now temporary vector is created by this procedure.
              std::vector<vector<int>>::iterator it =
                  std::find(IMStructure1.at(iloop3).begin(),
                            IMStructure1.at(iloop3).end(), testing);
              IMStructure1.at(iloop3).erase(it);
              testing.clear();
              // now the one of the leg is eliminated in the list, so there is
              // no overlapped leg to be tossed into baryon formation list.
            }
          }
        }
      }
    }
  }

  // dignostic measure{IMStructure1}
  /*JSINFO << "1.75 IMStructure printout";
  std::cout <<endl;
  for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
    for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
      std::cout <<" ( ";
      for(int icheck3 = 0; icheck3 <
  IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout
  <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
      }
      std::cout <<" ) ";
    }
    std::cout <<endl;
  }*/

  for (int irepair = 0; irepair < IMStructure1.size(); irepair++) {
    if (IMStructure1.at(irepair).size() == 1 &&
        IMStructure1.at(irepair).at(0).at(3) ==
            0) {  // size one means there are no shared legs to the junction, so
                  // put them directly to the single junction vector
      Singlejunction2.push_back(
          JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(0));
      Singlejunction2.push_back(
          JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(1));
      Singlejunction2.push_back(
          JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(2));
      Singlejunction1.push_back(Singlejunction2);
      Singlejunction2.clear();
    }  // After all, all these vectors of junction legs will be tossed into
       // PYTHIA to be hadronized. before that we need to set Mother Daughter
       // tag for PYTHIA
  }    // collecting all single junctions

  // dignostic measure{IMStructure1}
  /*JSINFO << "Second IMStructure printout";
  std::cout <<endl;
  for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
    for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
      std::cout <<" ( ";
      for(int icheck3 = 0; icheck3 <
  IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout
  <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
      }
      std::cout <<" ) ";
    }
    std::cout <<endl;
  }*/

  // now we need to form fake partons and baryons for appending particles to
  // PYTHIA{for Dijunction Structure} declaration of required Variables
  bool J = false;
  bool antiJ = false;
  bool bothJ = false;
  for (int idj1 = 0; idj1 < Dijunction1.size();
       idj1++) {  // searching dijunction
    for (int idj2 = 0; idj2 < Dijunction1.at(idj1).size();
         idj2++) {  // searching dijunction's leg
      // std::cout <<endl<<"let's check these particles "<<endl;
      // std::cout <<" { "<<Dijunction1.at(idj1).at(idj2).at(0).col()<<" ,
      // "<<Dijunction1.at(idj1).at(idj2).at(0).acol()<<" } "<<endl; std::cout
      // <<" { "<<Dijunction1.at(idj1).at(idj2).back().col()<<" ,
      // "<<Dijunction1.at(idj1).at(idj2).back().acol()<<" } "<<endl;
      for (int itpj1 = 0; itpj1 < IMStructure1.size();
           itpj1++) {  // searching Tempjunction
        for (int itpj2 = 1; itpj2 < Tempjunctions.at(itpj1).size();
             itpj2++) {  // searching color tags in tempjunction{ {_+1, 0 },
                         // {_+1, tag1}, {_+1, tag2}, {_+1, tag3} }
          if (Dijunction1.at(idj1).at(idj2).at(0).col() ==
                  Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0))
                      .at(itpj2)
                      .at(1) ||  // check whether first or last particle has
                                 // correponding color tags for J
              Dijunction1.at(idj1).at(idj2).back().col() ==
                  Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0))
                      .at(itpj2)
                      .at(1)) {
            J = true;
          }
          if (Dijunction1.at(idj1).at(idj2).at(0).acol() ==
                  Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0))
                      .at(itpj2)
                      .at(1) ||  // check whether first or last particle has
                                 // correponding color tags for J
              Dijunction1.at(idj1).at(idj2).back().acol() ==
                  Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0))
                      .at(itpj2)
                      .at(1)) {
            antiJ = true;
          }
        }
      }
      if (J == true && antiJ == false) {  // this leg is a part of junction
        DijunctionInfo2.push_back(1);
      }
      if (J == false && antiJ == true) {  // this leg is a part of junction
        DijunctionInfo2.push_back(-1);
      }
      if (J == true && antiJ == true) {
        bothJ = true;
      }
      if (bothJ) {  // this leg is shared leg, so that fake parton shoud be
                    // added after all
        DijunctionInfo2.push_back(0);
      }
      J = false;  // reseting variables to check other legs in Dijunction
                  // structure
      antiJ = false;
      bothJ = false;
    }
    DijunctionInfo1.push_back(DijunctionInfo2);
    DijunctionInfo2.clear();
  }  // so far, vector about the information of dijunction structure is formed,{
     // so, we could put identity 1-leg first, and put identity-0 leg and -1
     // leg}
  // based on the information above, we need to add fake partons, which are used
  // for telling pythia about the color flow{fake particle id is ignored, since
  // we'll set
  //  negative status flag}, with this color flow information, PYTHIA can detect
  //  J and anti-J

  // additional proceudre to change the configuration of shared leg by searching
  // through dijunction1 and dijunction1info vectors
  //  make sure that the shared leg is ordered from junction to antijunction,
  //  used later
  for (int ileg1 = 0; ileg1 < Dijunction1.size(); ileg1++) {
    for (int ileg2 = 0; ileg2 < Dijunction1.at(ileg1).size(); ileg2++) {
      bool needflip = false;
      if ((DijunctionInfo1.at(ileg1).at(ileg2) == 0 &&
           Dijunction1.at(ileg1).at(ileg2).size() != 1) &&
          (Dijunction1.at(ileg1).at(ileg2).at(0).col() ==
           Dijunction1.at(ileg1)
               .at(ileg2)
               .at(1)
               .acol())) {  // In this case, we need to flip the order of the
                            // partons in this shared leg
        needflip = true;
      }
      if (needflip) {
        std::reverse(Dijunction1.at(ileg1).at(ileg2).begin(),
                     Dijunction1.at(ileg1).at(ileg2).end());
      }
      needflip = false;
    }
  }

  // dignostic measure{DijunctionInfo1}
  /*std::cout <<endl;
  std::cout <<" DijunctionInfo Checking"<<endl;
  for(int icheck1 = 0; icheck1 < DijunctionInfo1.size(); icheck1++){
    std::cout <<" Let's Check all the elements! : ";
    for(int icheck2 = 0; icheck2 < DijunctionInfo1.at(icheck1).size();
  icheck2++){ std::cout <<DijunctionInfo1.at(icheck1).at(icheck2)<<" , ";
    }
    std::cout <<endl;
  }*/

  // dignostic measure{Dijunction1}
  /*std::cout <<endl;
  std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Dijunction1.size(); icheck1++){
    std::cout <<" Dijunction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Dijunction1.at(icheck1).size(); icheck2++){
      std::cout <<" Leg "<<icheck2<<" : "<<"{ identity : "<<
  DijunctionInfo1.at(icheck1).at(icheck2)<<" }  "; for(int icheck3 = 0; icheck3
  < Dijunction1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout <<" (
  "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" ,
  "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      std::cout <<endl;
    }
  }*/

  // dignostic measure{Singlejunction1}
  /*std::cout <<endl;
  std::cout <<" List of Single Junction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Singlejunction1.size(); icheck1++){
    std::cout <<" Single junction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Singlejunction1.at(icheck1).size();
  icheck2++){ std::cout <<" Leg "<<icheck2<<" : "; for(int icheck3 = 0; icheck3
  < Singlejunction1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout <<" (
  "<<Singlejunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" ,
  "<<Singlejunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      std::cout <<endl;
    }
  }*/

  // For the validity with PYTHIA running, each Leg in Single or Dijunction
  // Structure should have q or q-bar as endpoint particle, so checking whether
  // the gluon is located at the endpoint of the Leg Starting from Leg's in
  // single junction.
  for (int is1 = 0; is1 < Singlejunction1.size(); is1++) {
    HHparton p1 =
        Singlejunction1.at(is1).at(0).back();  // endpoint particle of first leg
                                               // in is1_st Singlejunction1
    HHparton p2 =
        Singlejunction1.at(is1).at(1).back();  // endpoint particle of second
                                               // leg in is1_st Singlejunction1
    HHparton p3 =
        Singlejunction1.at(is1).at(2).back();  // endpoint particle of third leg
                                               // in is1_st Singlejunction1

    // basically, assume the junction. but check the id's of them and change
    // identity of the condition for antijunction is meeted
    bool ajunction = false;

    if (p1.id() < 0 || p2.id() < 0 || p3.id() < 0) {
      ajunction = true;
    }

    for (int is2 = 0; is2 < Singlejunction1.at(is1).size(); is2++) {
      HHparton endpoint = Singlejunction1.at(is1).at(is2).back();

      if (endpoint.col() != 0 && endpoint.acol() != 0 && ajunction) {
        int loc = 0;
        HHparton tempparton1 = endpoint;
        if (tempparton1.id() == 21) {
          tempparton1.id(1);
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
          tempparton1.id(21);
        } else {
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
        }
        if (loc == 999999999 || loc > 0) {
          int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3);
          HHparton fakep = tempparton1;
          fakep.id(fid);
          fakep.set_color(0);
          fakep.set_anti_color(endpoint.col());
          if (std::abs(fakep.id()) < 3) {
            fakep.mass(xmq);
          } else {
            fakep.mass(xms);
          }
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          SP_remnants.add(fakep);
          Singlejunction1.at(is1).at(is2).push_back(fakep);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].acol(endpoint.col());
          HH_thermal[-loc - 1].col(0);
          HH_thermal[-loc - 1].is_used(true);
          HH_thermal[-loc - 1].is_remnant(true);
          SP_remnants.add(HH_thermal[-loc - 1]);
          SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
          Singlejunction1.at(is1).at(is2).push_back(HH_thermal[-loc - 1]);
        }
      }
    }

    for (int is3 = 0; is3 < Singlejunction1.at(is1).size(); is3++) {
      HHparton endpoint = Singlejunction1.at(is1).at(is3).back();

      if (endpoint.col() != 0 && endpoint.acol() != 0 && !ajunction) {
        int loc = 0;
        HHparton tempparton1 = endpoint;
        if (tempparton1.id() == 21) {
          tempparton1.id(-1);
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
          tempparton1.id(21);
        } else {
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
        }
        if (loc == 999999999 || loc > 0) {
          int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3);
          HHparton fakep = tempparton1;
          fakep.id(fid);
          fakep.set_color(endpoint.acol());
          fakep.set_anti_color(0);
          if (std::abs(fakep.id()) < 3) {
            fakep.mass(xmq);
          } else {
            fakep.mass(xms);
          }
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          SP_remnants.add(fakep);
          Singlejunction1.at(is1).at(is3).push_back(fakep);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].col(endpoint.acol());
          HH_thermal[-loc - 1].acol(0);
          HH_thermal[-loc - 1].is_used(true);
          HH_thermal[-loc - 1].is_remnant(true);
          SP_remnants.add(HH_thermal[-loc - 1]);
          SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
          Singlejunction1.at(is1).at(is3).push_back(HH_thermal[-loc - 1]);
        }
      }
    }
    ajunction = false;  // resetting the variable.
  }

  // Working for Leg's in dijunction structure.
  for (int idj1 = 0; idj1 < Dijunction1.size(); idj1++) {
    for (int idj2 = 0; idj2 < Dijunction1.at(idj1).size(); idj2++) {
      HHparton endpoint = Dijunction1.at(idj1).at(idj2).back();

      if (DijunctionInfo1.at(idj1).at(idj2) == -1 && endpoint.col() != 0 &&
          endpoint.acol() != 0) {
        int loc = 0;
        HHparton tempparton1 = endpoint;
        if (tempparton1.id() == 21) {
          tempparton1.id(1);
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
          tempparton1.id(21);
        } else {
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
        }
        if (loc == 999999999 || loc > 0) {
          int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? -1 : -2) : -3);
          HHparton fakep = tempparton1;
          fakep.id(fid);
          fakep.set_color(0);
          fakep.set_anti_color(endpoint.col());
          if (std::abs(fakep.id()) < 3) {
            fakep.mass(xmq);
          } else {
            fakep.mass(xms);
          }
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          SP_remnants.add(fakep);
          Dijunction1.at(idj1).at(idj2).push_back(fakep);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].acol(endpoint.col());
          HH_thermal[-loc - 1].col(0);
          HH_thermal[-loc - 1].is_used(true);
          HH_thermal[-loc - 1].is_remnant(true);
          SP_remnants.add(HH_thermal[-loc - 1]);
          SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
          Dijunction1.at(idj1).at(idj2).push_back(HH_thermal[-loc - 1]);
        }
      }
    }

    for (int idj3 = 0; idj3 < Dijunction1.at(idj1).size(); idj3++) {
      HHparton endpoint = Dijunction1.at(idj1).at(idj3).back();

      if (DijunctionInfo1.at(idj1).at(idj3) == 1 && endpoint.col() != 0 &&
          endpoint.acol() != 0) {
        int loc = 0;
        HHparton tempparton1 = endpoint;
        if (tempparton1.id() == 21) {
          tempparton1.id(-1);
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
          tempparton1.id(21);
        } else {
          loc = findcloserepl(tempparton1, 1, false, true, HH_showerptns,
                              HH_thermal);
        }
        if (loc == 999999999 || loc > 0) {
          int fid = ((ran() > 0.33333333) ? ((ran() > 0.5) ? 1 : 2) : 3);
          HHparton fakep = tempparton1;
          fakep.id(fid);
          fakep.set_color(endpoint.acol());
          fakep.set_anti_color(0);
          if (std::abs(fakep.id()) < 3) {
            fakep.mass(xmq);
          } else {
            fakep.mass(xms);
          }
          double fake_pT = (p_fake >= 1.) ? 0.282842712474619 : 0.;
          double fake_phi = 2. * 3.14159265358979 * ran();
          fakep.px(fake_pT * cos(fake_phi));
          fakep.py(fake_pT * sin(fake_phi));
          if (number_p_fake == 0 || number_p_fake == 2) {
            fakep.pz(-p_fake);
          } else if (number_p_fake == 1 || number_p_fake == 3) {
            fakep.pz(p_fake);
          } else {
            fakep.pz(0.);
          }
          number_p_fake++;
          fakep.e(std::sqrt(fakep.px() * fakep.px() + fakep.py() * fakep.py() +
                            fakep.pz() * fakep.pz() +
                            fakep.mass() * fakep.mass()));
          fakep.orig(-1);
          fakep.is_remnant(true);
          fakep.is_fakeparton(true);
          SP_remnants.add(fakep);
          Dijunction1.at(idj1).at(idj3).push_back(fakep);
        } else if (loc < 0) {
          HH_thermal[-loc - 1].col(endpoint.acol());
          HH_thermal[-loc - 1].acol(0);
          HH_thermal[-loc - 1].is_used(true);
          HH_thermal[-loc - 1].is_remnant(true);
          SP_remnants.add(HH_thermal[-loc - 1]);
          SP_remnants[SP_remnants.num() - 1].par(-loc - 1);
          Dijunction1.at(idj1).at(idj3).push_back(HH_thermal[-loc - 1]);
        }
      }
    }
  }

  // since the repeating is needed for completing the junction structure, It is
  // inevitable that there are repeated particle in Recombeary1, It's
  // complicated to explain, and useless to focus on.
  //  so that just delete repeated particles by check color Tags
  for (int irb1 = 0; irb1 < Recombearly1.size(); irb1++) {  // pick one baryon
    for (int irb2 = 0; irb2 < Recombearly1.at(irb1).size();
         irb2++) {  // pick one particle in baryon
      for (int irb3 = irb2 + 1; irb3 < Recombearly1.at(irb1).size();
           irb3++) {  // checking through all particles in baryon
        if ((Recombearly1.at(irb1).at(irb2).col() ==
             Recombearly1.at(irb1).at(irb3).col()) &&
            Recombearly1.at(irb1).at(irb2).acol() ==
                Recombearly1.at(irb1).at(irb3).acol()) {
          Recombearly1.at(irb1).erase(Recombearly1.at(irb1).begin() + irb3);
          irb3--;
        }
      }
    }
  }

  // dignositic measure{early recombined baryon}
  /*std::cout <<endl;
  std::cout <<" number of early recombined baryons :
  "<<Recombearly1.size()<<endl; for(int ibary1 = 0; ibary1 <
  Recombearly1.size(); ibary1++){ std::cout <<"the "<< ibary1 <<" st baryon is
  made of "; for(int ibary2 = 0; ibary2 < Recombearly1.at(ibary1).size();
  ibary2++){ std::cout <<" { "<<Recombearly1.at(ibary1).at(ibary2).col()<<" ,
  "<< Recombearly1.at(ibary1).at(ibary2).acol()<<" , "<<
  Recombearly1.at(ibary1).at(ibary2).id()<<" , "<<
  Recombearly1.at(ibary1).at(ibary2).e()  <<" } , ";
    }
    std::cout <<endl;
  }*/

  // force recombine baryons in Recombearly1
  for (int irb = 0; irb < Recombearly1.size(); irb++) {
    double hbac2 = hbarc * hbarc;

    FourVector Pbaryon;
    Pbaryon.Set(Recombearly1.at(irb)[0].px() + Recombearly1.at(irb)[1].px() +
                    Recombearly1.at(irb)[2].px(),
                Recombearly1.at(irb)[0].py() + Recombearly1.at(irb)[1].py() +
                    Recombearly1.at(irb)[2].py(),
                Recombearly1.at(irb)[0].pz() + Recombearly1.at(irb)[1].pz() +
                    Recombearly1.at(irb)[2].pz(),
                0.);

    // baryon(CM) velocity
    FourVector betaB;  // really p[i]/e below
    betaB.Set(Pbaryon.x() /
                  (Recombearly1.at(irb)[0].e() + Recombearly1.at(irb)[1].e() +
                   Recombearly1.at(irb)[2].e()),
              Pbaryon.y() /
                  (Recombearly1.at(irb)[0].e() + Recombearly1.at(irb)[1].e() +
                   Recombearly1.at(irb)[2].e()),
              Pbaryon.z() /
                  (Recombearly1.at(irb)[0].e() + Recombearly1.at(irb)[1].e() +
                   Recombearly1.at(irb)[2].e()),
              0.);
    betaB.Set(betaB.x(), betaB.y(), betaB.z(),
              1. / (sqrt(1. - (betaB.x() * betaB.x() + betaB.y() * betaB.y() +
                               betaB.z() * betaB.z()))));

    // boosting into CM frame
    FourVector pos_BCM[3], p_BCM[3];
    pos_BCM[0] = Recombearly1.at(irb)[0].boost_pos(betaB);
    pos_BCM[1] = Recombearly1.at(irb)[1].boost_pos(betaB);
    pos_BCM[2] = Recombearly1.at(irb)[2].boost_pos(betaB);
    p_BCM[0] = Recombearly1.at(irb)[0].boost_P(betaB);
    p_BCM[1] = Recombearly1.at(irb)[1].boost_P(betaB);
    p_BCM[2] = Recombearly1.at(irb)[2].boost_P(betaB);

    // velocities in CM frame
    FourVector v_BCM[3];
    v_BCM[0].Set(p_BCM[0].x() / p_BCM[0].t(), p_BCM[0].y() / p_BCM[0].t(),
                 p_BCM[0].z() / p_BCM[0].t(), 0.);  // these are really p[i]/e
    v_BCM[1].Set(p_BCM[1].x() / p_BCM[1].t(), p_BCM[1].y() / p_BCM[1].t(),
                 p_BCM[1].z() / p_BCM[1].t(), 0.);
    v_BCM[2].Set(p_BCM[2].x() / p_BCM[2].t(), p_BCM[2].y() / p_BCM[2].t(),
                 p_BCM[2].z() / p_BCM[2].t(), 0.);

    // propagating quarks until time of youngest quark
    double curtime =
        std::max(std::max(pos_BCM[0].t(), pos_BCM[1].t()), pos_BCM[2].t());
    FourVector cur_pos[3];
    cur_pos[0].Set(pos_BCM[0].x() + v_BCM[0].x() * (curtime - pos_BCM[0].t()),
                   pos_BCM[0].y() + v_BCM[0].y() * (curtime - pos_BCM[0].t()),
                   pos_BCM[0].z() + v_BCM[0].z() * (curtime - pos_BCM[0].t()),
                   curtime);
    cur_pos[1].Set(pos_BCM[1].x() + v_BCM[1].x() * (curtime - pos_BCM[1].t()),
                   pos_BCM[1].y() + v_BCM[1].y() * (curtime - pos_BCM[1].t()),
                   pos_BCM[1].z() + v_BCM[1].z() * (curtime - pos_BCM[1].t()),
                   curtime);
    cur_pos[2].Set(pos_BCM[2].x() + v_BCM[2].x() * (curtime - pos_BCM[2].t()),
                   pos_BCM[2].y() + v_BCM[2].y() * (curtime - pos_BCM[2].t()),
                   pos_BCM[2].z() + v_BCM[2].z() * (curtime - pos_BCM[2].t()),
                   curtime);

    // finding position of CM at curtime
    FourVector pos_CM;
    pos_CM.Set(
        (cur_pos[0].x() * Recombearly1.at(irb)[0].mass() +
         cur_pos[1].x() * Recombearly1.at(irb)[1].mass() +
         cur_pos[2].x() * Recombearly1.at(irb)[2].mass()) /
            (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass() +
             Recombearly1.at(irb)[2].mass()),
        (cur_pos[0].y() * Recombearly1.at(irb)[0].mass() +
         cur_pos[1].y() * Recombearly1.at(irb)[1].mass() +
         cur_pos[2].y() * Recombearly1.at(irb)[2].mass()) /
            (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass() +
             Recombearly1.at(irb)[2].mass()),
        (cur_pos[0].z() * Recombearly1.at(irb)[0].mass() +
         cur_pos[1].z() * Recombearly1.at(irb)[1].mass() +
         cur_pos[2].z() * Recombearly1.at(irb)[2].mass()) /
            (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass() +
             Recombearly1.at(irb)[2].mass()),
        curtime);

    // finding position of baryon in lab frame
    betaB.Set(-betaB.x(), -betaB.y(), -betaB.z(), betaB.t());
    FourVector pos_lab = HHboost(betaB, pos_CM);

    // finding relative positions of partons in CM frame
    FourVector pos_rel_square[2];
    pos_rel_square[0].Set((cur_pos[0].x() - cur_pos[1].x()) / sqrt(2.),
                          (cur_pos[0].y() - cur_pos[1].y()) / sqrt(2.),
                          (cur_pos[0].z() - cur_pos[1].z()) / sqrt(2.), 0.);
    pos_rel_square[1].Set(
        ((cur_pos[0].x() * Recombearly1.at(irb)[0].mass() +
          cur_pos[1].x() * Recombearly1.at(irb)[1].mass()) /
             (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass()) -
         cur_pos[2].x()) *
            sqrt(2. / 3.),
        ((cur_pos[0].y() * Recombearly1.at(irb)[0].mass() +
          cur_pos[1].y() * Recombearly1.at(irb)[1].mass()) /
             (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass()) -
         cur_pos[2].y()) *
            sqrt(2. / 3.),
        ((cur_pos[0].z() * Recombearly1.at(irb)[0].mass() +
          cur_pos[1].z() * Recombearly1.at(irb)[1].mass()) /
             (Recombearly1.at(irb)[0].mass() + Recombearly1.at(irb)[1].mass()) -
         cur_pos[2].z()) *
            sqrt(2. / 3.),
        0.);

    // so far, values to form a baryon are set, now make hadron objects based on
    // this information now we're forming the hadron
    HHhadron formedhadron;
    // setting the hadron values: is a recombined hadron, mass, and parents
    formedhadron.is_recohad(true);
    formedhadron.mass(p_BCM[0].t() + p_BCM[1].t() + p_BCM[2].t());

    // setting hadron position and momentum vectors
    Pbaryon.Set(Pbaryon.x(), Pbaryon.y(), Pbaryon.z(),
                sqrt(Pbaryon.x() * Pbaryon.x() + Pbaryon.y() * Pbaryon.y() +
                     Pbaryon.z() * Pbaryon.z() +
                     formedhadron.mass() * formedhadron.mass()));
    formedhadron.pos(pos_lab);
    formedhadron.P(Pbaryon);

    // need to choose *what* hadron we've formed... base this on the parton
    // id's, mass, & if excited might want to do this differently? void
    // f'n(partoncollection, formedhadron)? since set_baryon_id function works
    // with object of parton_colection class, make temporary parton_collection
    parton_collection tempobject;  // since the generator resets all components,
                                   // don't need to reset the vactors in it
    tempobject.add(Recombearly1.at(irb)[0]);
    tempobject.add(Recombearly1.at(irb)[1]);
    tempobject.add(Recombearly1.at(irb)[2]);

    set_baryon_id(tempobject, formedhadron);

    formedhadron.is_shth(thermal_parton_junction[irb]);

    // need to add the hadron to the collection
    HH_hadrons.add(formedhadron);
  }

  // from now we need re indexing a of the particles in SP_remnants by comparing
  // color tag in Juncstructure component. Maybe there would not be a problem in
  // the used_junction value, but take most defensive measure for the possible
  // error.
  std::vector<int> checking;
  bool indicecheck = false;
  for (int irem = 0; irem < SP_remnants.num(); irem++) {
    for (int ijs1 = 0; ijs1 < JuncStructure.size(); ijs1++) {
      for (int ijs2 = 0; ijs2 < JuncStructure.at(ijs1).size(); ijs2++) {
        for (int ijs3 = 0; ijs3 < JuncStructure.at(ijs1).at(ijs2).size();
             ijs3++) {
          if ((SP_remnants[irem].col() ==
               JuncStructure.at(ijs1).at(ijs2).at(ijs3).col()) &&
              (SP_remnants[irem].acol() ==
               JuncStructure.at(ijs1).at(ijs2).at(ijs3).acol())) {
            SP_remnants[irem].used_junction(true);
            checking.push_back(irem);
          }
        }
      }
    }
  }
  for (int irem = 0; irem < SP_remnants.num(); irem++) {
    for (int icheck = 0; icheck < checking.size(); icheck++) {
      if (irem == checking.at(icheck)) {
        indicecheck = true;
      }
    }
    if (indicecheck = false) {
      SP_remnants[irem].used_junction(false);
    }
  }
  // Now set up for junction is done, let's check the left partons to establish
  // string structure

  // std::cout <<endl<<"let's check left particles!"<<endl;
  /*for(int irp = 0 ; irp < SP_remnants.num(); irp++) {
    if(SP_remnants[irp].used_junction() == false){
      std::cout << "( "<<SP_remnants[irp].id()<<" , "<<SP_remnants[irp].col()<<"
  , "<<SP_remnants[irp].acol()<<" ) "<<endl;
    }
  }*/

  for (int ileft = 0; ileft < SP_remnants.num(); ileft++) {
    if (SP_remnants[ileft].used_junction() == false) {
      finalstring.push_back(SP_remnants[ileft]);
    }
  }

  // contains the chopped off legs of junctions
  for (int itail1 = 0; itail1 < Tailoredstring1.size(); itail1++) {
    for (int itail2 = 0; itail2 < Tailoredstring1.at(itail1).size(); itail2++) {
      finalstring.push_back(Tailoredstring1.at(itail1).at(itail2));
    }
  }

  // std::cout <<endl<<"check final particles!"<<endl;
  /*for(int ifin = 0; ifin < finalstring.size(); ifin++){
    std::cout << finalstring.at(ifin).id() << "," << finalstring.at(ifin).col()
  << "," << finalstring.at(ifin).acol() <<endl;
  }*/

  /*std::cout <<endl<<"check final particles!"<<endl;
  for(int ifin = 0; ifin < finalstring.size(); ifin++){
    std::cout <<" { "<< finalstring.at(ifin).col() << " , " <<
  finalstring.at(ifin).acol() << " } "<<endl;
  }*/

  // so far, particles in finalstring are ready for being transfered into
  // PYTHIA, since we don't need to set mother, daughter tag. but, still need to
  // care about single junction and dijunction system, especially for the mother
  // daughter tag!
  //  mother, daughter tags are working with the execution order{= order of
  //  being declared in pythia} so before tossing the datum to PYTHIA, rearrange
  //  order in the vector declare the transit space for these particles, while
  //  moving onto vector to vector, partons are assigned mother, daughter tag
  //  tor PYTHIA running{definition of that tags are given well in main21.cc
  //  file in   /installed pythia folder/bin/example/ }
  vector<vector<HHparton>> Transitdijunction1;
  vector<HHparton>
      Transitdijunction2;  // two fake mother B, B-bar and fake q,qbar added and
                           // give mother daughter tags!
  vector<vector<vector<HHparton>>>
      Tempsorting1;  // vectors for rearranging legs in dijunction{Desired
                     // arrangement : identity 1, 1, 0 , -1, -1}
  vector<vector<HHparton>> Tempsorting2;
  vector<vector<HHparton>> Transitsinglejunction1;
  vector<HHparton> Transitsinglejunction2;
  vector<HHparton> WaitingLineforPY;  // final list of partons with complete M,D
                                      // tag and col tags

  // set the value for incrementation in mother daughter tag.
  // and start from dijunction, first, check the identity{1: junction leg, -1:
  // anti_junction leg, 0:shared leg}, since we'll read col tag first, consider
  // identity=1 leg first and then 0, -1
  for (int dijuncfin1 = 0; dijuncfin1 < Dijunction1.size();
       dijuncfin1++) {  // Going to find shared leg and attach fake partons to
                        // the begin and the end of the leg
    for (int dijuncfin2 = 0; dijuncfin2 < 5;
         dijuncfin2++) {  // inspect through five legs
      if (DijunctionInfo1.at(dijuncfin1).at(dijuncfin2) ==
          0) {  // check whether it's shared leg{identity =0} and add
                // fakepartons to the first and second posittion in the vector
        HHparton fakeq;
        HHparton fakeqbar;
        fakeq.id(1);
        fakeq.set_color(Dijunction1.at(dijuncfin1).at(dijuncfin2).at(0).col());
        fakeq.PY_stat(-21);
        fakeq.mass(xmq);
        fakeq.e(xmq);
        fakeq.orig(-1);
        fakeq.px(0.);
        fakeq.py(0.);
        fakeq.pz(0.);
        int endpoint = Dijunction1.at(dijuncfin1).at(dijuncfin2).size() - 1;
        fakeqbar.id(-1);
        fakeqbar.set_anti_color(
            Dijunction1.at(dijuncfin1).at(dijuncfin2).at(endpoint).acol());
        fakeqbar.PY_stat(-21);
        fakeqbar.mass(xmq);
        fakeqbar.e(xmq);
        fakeqbar.orig(-1);
        fakeqbar.px(0.);
        fakeqbar.py(0.);
        fakeqbar.pz(0.);
        std::vector<HHparton>::iterator it2 =
            Dijunction1.at(dijuncfin1).at(dijuncfin2).begin();
        Dijunction1.at(dijuncfin1).at(dijuncfin2).insert(it2, fakeqbar);
        std::vector<HHparton>::iterator it1 =
            Dijunction1.at(dijuncfin1).at(dijuncfin2).begin();
        Dijunction1.at(dijuncfin1).at(dijuncfin2).insert(it1, fakeq);
      }
    }
    for (int isort1 = 0; isort1 < 5; isort1++) {
      if (DijunctionInfo1.at(dijuncfin1).at(isort1) == 1) {
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    for (int isort1 = 0; isort1 < 5; isort1++) {
      if (DijunctionInfo1.at(dijuncfin1).at(isort1) == 0) {
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    for (int isort1 = 0; isort1 < 5; isort1++) {
      if (DijunctionInfo1.at(dijuncfin1).at(isort1) == -1) {
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    Tempsorting1.push_back(Tempsorting2);
    Tempsorting2.clear();
  }  // first looping to add fake partons is Finished, now sort the order of the
     // legs to be identity 1,1,0,-1,-1

  // dignostic measure{Dijunction1}
  /*std::cout <<endl;
  std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Dijunction1.size(); icheck1++){
    std::cout <<" Dijunction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Dijunction1.at(icheck1).size(); icheck2++){
      std::cout <<" Leg "<<icheck2<<" : "<<"{ identity : "<<
  DijunctionInfo1.at(icheck1).at(icheck2)<<" }  "; for(int icheck3 = 0; icheck3
  < Dijunction1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout <<" (
  "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" ,
  "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      std::cout <<endl;
    }
  }*/

  // dignostic measure{Tempsorting}
  /*std::cout <<endl;
  std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Tempsorting1.size(); icheck1++){
    std::cout <<" Tempsorted Dijunction1 : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Tempsorting1.at(icheck1).size(); icheck2++){
      for(int icheck3 = 0; icheck3 <
  Tempsorting1.at(icheck1).at(icheck2).size(); icheck3++){ std::cout <<" (
  "<<Tempsorting1.at(icheck1).at(icheck2).at(icheck3).col()<<" ,
  "<<Tempsorting1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  "<<" ,
  "<<Tempsorting1.at(icheck1).at(icheck2).at(icheck3).x_t();
      }
      std::cout <<endl;
    }
  }*/

  // adding fake partons and sorting them with right order{Junctionleg1,2 shared
  // leg, antijunctionleg1,2} is finished so far, Now form a fake mothers to
  // prepare for PYTHIA
  for (int itag1 = 0; itag1 < Tempsorting1.size();
       itag1++) {  /// we need to link Tempsorting1 with Transitdijunction1,
    // First, put fake B,Bbar into 1st and 2nd position of Transitdijunction2
    parton_collection FakeBaryonElements;  // array of quarks in fake B
    parton_collection
        FakeAntibaryonElements;  // array of anti quarks in fake Bbar

    FakeBaryonElements.add(
        Tempsorting1.at(itag1).at(0).back());  // quarks to form fake FakeBaryon
    FakeBaryonElements.add(Tempsorting1.at(itag1).at(1).back());
    FakeBaryonElements.add(Tempsorting1.at(itag1).at(2).at(0));
    // for the corrspondence with pdg particle id, sorting the these ids based
    // on absolute value.

    for (int iswap = 0; iswap < 3; iswap++) {
      if (abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[1].id())) {
        std::swap(FakeBaryonElements[2], FakeBaryonElements[1]);
      }
      if (abs(FakeBaryonElements[1].id()) > abs(FakeBaryonElements[0].id())) {
        std::swap(FakeBaryonElements[1], FakeBaryonElements[0]);
      }
      if (abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[0].id())) {
        std::swap(FakeBaryonElements[2], FakeBaryonElements[0]);
      }
    }

    FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(2).at(
        1));  // anti quarks to form FakeAntiBaryon
    FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(3).back());
    FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(4).back());

    for (int iswap = 0; iswap < 3; iswap++) {
      if (abs(FakeAntibaryonElements[2].id()) >
          abs(FakeAntibaryonElements[1].id())) {
        std::swap(FakeAntibaryonElements[2], FakeAntibaryonElements[1]);
      }
      if (abs(FakeAntibaryonElements[1].id()) >
          abs(FakeAntibaryonElements[0].id())) {
        std::swap(FakeAntibaryonElements[1], FakeAntibaryonElements[0]);
      }
      if (abs(FakeAntibaryonElements[2].id()) >
          abs(FakeAntibaryonElements[0].id())) {
        std::swap(FakeAntibaryonElements[2], FakeAntibaryonElements[0]);
      }
    }

    HHhadron store_id_hadron1;
    set_baryon_id(FakeBaryonElements, store_id_hadron1);
    HHhadron store_id_hadron2;
    set_baryon_id(FakeAntibaryonElements, store_id_hadron2);

    FakeBaryonElements.clear();
    FakeAntibaryonElements.clear();

    HHparton fakeB;
    HHparton fakeBbar;
    fakeB.id(store_id_hadron1.id());
    fakeB.PY_stat(-11);
    fakeB.mass(3 * xmq);
    fakeB.e(3 * xmq);
    fakeB.px(0.);
    fakeB.py(0.);
    fakeB.pz(0.);

    fakeBbar.id(store_id_hadron2.id());
    fakeBbar.PY_stat(-11);
    fakeBbar.mass(3 * xmq);
    fakeBbar.e(3 * xmq);
    fakeBbar.px(0.);
    fakeBbar.py(0.);
    fakeBbar.pz(0.);

    fakeB.PY_tag1(Tempsorting1.at(itag1).at(0).back().col());
    fakeB.PY_tag2(Tempsorting1.at(itag1).at(1).back().col());
    fakeB.PY_tag3(Tempsorting1.at(itag1).at(2).at(0).col());

    fakeBbar.PY_tag1(Tempsorting1.at(itag1).at(2).at(1).acol());
    fakeBbar.PY_tag2(Tempsorting1.at(itag1).at(3).back().acol());
    fakeBbar.PY_tag3(Tempsorting1.at(itag1).at(4).back().acol());

    // so far, we formed two fake partons, which would be at the center of
    // Junction and Antijunction, Now put them 0th and 1st position of
    // Transirdijunction1 vector
    Transitdijunction2.push_back(fakeB);
    Transitdijunction2.push_back(fakeBbar);
    Transitdijunction1.push_back(Transitdijunction2);
    Transitdijunction2.clear();
  }

  // the next step is to set M,D Tags
  int Intag = 0;  // tag for internal M,D tagging process
  for (int iMD1 = 0; iMD1 < Tempsorting1.size();
       iMD1++) {  // after working in these Legs, we will return the partons to
                  // Transitdijunction1,
    // so far, the order in the vectors{Dijunction1, DijunctionInfo1,
    // Tempsorting, Transitdijunction1} are unified, but the difference is the
    // value in each vectors
    vector<HHparton> Leg1 = Tempsorting1.at(iMD1).at(
        0);  // for convenience, reassigning legs in sorted dijunction
             // structure.{identity 1,1,0,-1,-1}
    vector<HHparton> Leg2 = Tempsorting1.at(iMD1).at(1);
    vector<HHparton> Leg3 = Tempsorting1.at(iMD1).at(2);
    vector<HHparton> Leg4 = Tempsorting1.at(iMD1).at(3);
    vector<HHparton> Leg5 = Tempsorting1.at(iMD1).at(4);
    // starting from Leg1{identity1}
    for (int ileg1 = 0; ileg1 < Leg1.size(); ileg1++) {
      Leg1.at(ileg1).PY_par1(
          0);  // the mother of first leg1, 2 should be first fakebaryon{located
               // 1st in the Transitdijunction2 vector!}
      Leg1.at(ileg1).PY_par2(0);
      Transitdijunction1.at(iMD1).push_back(Leg1.at(
          ileg1));  // toss Leg1 to the Transirdijunction1.at{iMD}, which would
                    // be the last step before WaitingLineforPY
    }
    for (int ileg2 = 0; ileg2 < Leg2.size(); ileg2++) {
      Leg2.at(ileg2).PY_par1(
          0);  // the mother of first leg1, 2 should be first fakebaryon{located
               // 1st in the Transitdijunction2 vector!}
      Leg2.at(ileg2).PY_par2(0);
      Transitdijunction1.at(iMD1).push_back(Leg2.at(ileg2));
    }
    // for leg3, because of first and end fake particles, it goes differently
    // from other legs
    Intag = 2 + Leg1.size() +
            Leg2.size();  // Previous tag + two fakemothers + Leg1 + Leg2
    // First of all, fake parton pairs are located for the convenience of MD
    // tagging
    Leg3.at(0).PY_par1(0);  // daughter of fakeB
    Leg3.at(0).PY_par2(0);
    Leg3.at(0).PY_dau1(Intag + Leg4.size() + Leg5.size() +
                       2);  // daughter tags for gluons between two fake qqbar
    Leg3.at(0).PY_dau2(Intag + Leg4.size() + Leg5.size() + Leg3.size() -
                       1);  // since
    Leg3.at(1).PY_par1(1);  // daughter of fakeBbar
    Leg3.at(1).PY_par2(1);
    Leg3.at(1).PY_dau1(Intag + Leg4.size() + Leg5.size() + 2);
    Leg3.at(1).PY_dau2(Intag + Leg4.size() + Leg5.size() + Leg3.size() - 1);

    Transitdijunction1.at(iMD1).push_back(Leg3.at(0));
    Transitdijunction1.at(iMD1).push_back(
        Leg3.at(1));  // push first two fake partons and append remnants at the
                      // last, this is for convenience

    Transitdijunction1.at(iMD1).at(0).PY_dau1(
        2);  // correcting daughter tags of fake B located at 0
    Transitdijunction1.at(iMD1).at(0).PY_dau2(Intag);
    Transitdijunction1.at(iMD1).at(1).PY_dau1(
        Intag + 1);  // setting starting daughter tag of fakeBbar located at 1
    Transitdijunction1.at(iMD1).at(1).PY_dau2(Intag + 1 + Leg4.size() +
                                              Leg5.size());

    for (int ileg4 = 0; ileg4 < Leg4.size(); ileg4++) {
      Leg4.at(ileg4).PY_par1(1);
      Leg4.at(ileg4).PY_par2(1);
      Transitdijunction1.at(iMD1).push_back(Leg4.at(ileg4));
    }
    for (int ileg5 = 0; ileg5 < Leg5.size(); ileg5++) {
      Leg5.at(ileg5).PY_par1(1);
      Leg5.at(ileg5).PY_par2(1);
      Transitdijunction1.at(iMD1).push_back(Leg5.at(ileg5));
    }
    for (int ileg3 = 2; ileg3 < Leg3.size();
         ileg3++) {  // mother tags for the gluons between two fake qqbar pair
      Leg3.at(ileg3).PY_par1(Intag);
      Leg3.at(ileg3).PY_par2(Intag + 1);
      Transitdijunction1.at(iMD1).push_back(Leg3.at(ileg3));
    }
  }  // Inner tagging loop is finished for dijunction system,{That is,
     // Information in Tempsorting1 is transfered to Transitdijunction1}

  // dignostic measure{Transitdijunction1}
  /*std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2
  >> dau1 >> dau2 "<< endl; for(int icheck1 = 0; icheck1 <
  Transitdijunction1.size(); icheck1++ ){ for(int i = 0; i <
  Transitdijunction1.at(icheck1).size(); i++){ vector<HHparton> temp =
  Transitdijunction1.at(icheck1); std::cout  << i <<"   "<< temp.at(i).id() <<"
  "<< temp.at(i).PY_stat() <<"   "<< temp.at(i).PY_par1() <<"   "<<
  temp.at(i).PY_par2() <<"   "
      << temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<
  temp.at(i).col()<< " , " << temp.at(i).acol() << " ) " <<endl;
    }
  }*/

  // Now,start working on Single Junction System.
  for (int iSJ = 0; iSJ < Singlejunction1.size();
       iSJ++) {  // this process is for trasferring information from
                 // Singlejunction1 into Transitsinglejunction1, and finally
                 // they will be located at WaitingLineforPY
    // based on the partons at the end of each Leg, add fake mother baryon
    // first.
    parton_collection FakeBaryonElements;  // array of quarks in fake mother

    HHparton q1 = Singlejunction1.at(iSJ).at(0).back();
    HHparton q2 = Singlejunction1.at(iSJ).at(1).back();
    HHparton q3 = Singlejunction1.at(iSJ).at(2).back();

    HHparton q1fin = q1;
    HHparton q2fin = q2;
    HHparton q3fin = q3;

    FakeBaryonElements.add(
        q1fin);  // add first particle in 1st leg in iSJ_st singlejunction
    FakeBaryonElements.add(
        q2fin);  // add first particle in 2nd leg in iSJ_st singlejunction
    FakeBaryonElements.add(
        q3fin);  // add first particle in 3rd leg in iSJ_st singlejunction

    // for the corrspondence with pdg particle id, sorting the these ids based
    // on absolute value
    for (int iswap = 0; iswap < 3; iswap++) {
      if (abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[1].id())) {
        std::swap(FakeBaryonElements[2], FakeBaryonElements[1]);
      }
      if (abs(FakeBaryonElements[1].id()) > abs(FakeBaryonElements[0].id())) {
        std::swap(FakeBaryonElements[1], FakeBaryonElements[0]);
      }
      if (abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[0].id())) {
        std::swap(FakeBaryonElements[2], FakeBaryonElements[0]);
      }
    }

    HHhadron store_id_hadron1;
    set_baryon_id(FakeBaryonElements, store_id_hadron1);

    int m_id = store_id_hadron1.id();
    if ((q1.id() < 0 && m_id > 0) || (q1.id() > 0 && m_id < 0)) {
      m_id *= -1;
    }

    HHparton fakeB;
    fakeB.id(m_id);
    fakeB.PY_stat(-11);
    fakeB.mass(3 * xmq);
    fakeB.e(3 * xmq);
    fakeB.px(0.);
    fakeB.py(0.);
    fakeB.pz(0.);
    if (q1.id() > 0) {
      fakeB.PY_tag1(q1.col());
      fakeB.PY_tag2(q2.col());
      fakeB.PY_tag3(q3.col());
    } else {
      fakeB.PY_tag1(q1.acol());
      fakeB.PY_tag2(q2.acol());
      fakeB.PY_tag3(q3.acol());
    }
    Transitsinglejunction2.push_back(fakeB);
    Transitsinglejunction1.push_back(Transitsinglejunction2);
    Transitsinglejunction2.clear();
  }  // so farfake mother is added to all single junction vectors{vectors of
     // partons}, so we need to remember all execution number is incremented by
     // one because of this fake baryons

  for (int iSJ = 0; iSJ < Singlejunction1.size(); iSJ++) {
    vector<HHparton> Leg1 = Singlejunction1.at(iSJ).at(0);
    vector<HHparton> Leg2 = Singlejunction1.at(iSJ).at(1);
    vector<HHparton> Leg3 = Singlejunction1.at(iSJ).at(2);

    for (int ileg1 = 0; ileg1 < Leg1.size(); ileg1++) {
      Leg1.at(ileg1).PY_par1(0);  // setting the mother tag as zero, which
                                  // indicates the first particle in the Transit
      Leg1.at(ileg1).PY_par2(0);
      Transitsinglejunction1.at(iSJ).push_back(Leg1.at(ileg1));
    }
    for (int ileg2 = 0; ileg2 < Leg2.size(); ileg2++) {
      Leg2.at(ileg2).PY_par1(0);  // setting the mother tag as zero, which
                                  // indicates the first particle in the Transit
      Leg2.at(ileg2).PY_par2(0);
      Transitsinglejunction1.at(iSJ).push_back(Leg2.at(ileg2));
    }
    for (int ileg3 = 0; ileg3 < Leg3.size(); ileg3++) {
      Leg3.at(ileg3).PY_par1(0);  // setting the mother tag as zero, which
                                  // indicates the first particle in the Transit
      Leg3.at(ileg3).PY_par2(0);
      Transitsinglejunction1.at(iSJ).push_back(Leg3.at(ileg3));
    }
    Transitsinglejunction1.at(iSJ).at(0).PY_dau1(
        1);  // dau tag starting from 1 {after zero, which means fake mother
             // baryon}
    Transitsinglejunction1.at(iSJ).at(0).PY_dau2(
        Leg1.size() + Leg2.size() +
        Leg3.size());  // dau tag starting from 1 {after zero, which means fake
                       // mother baryon}
  }  // At last, iSJ_st Singlejunction1 partons are moved to  iSJ_st vector
     // component in Transitsinglejunction1

  // dignostic measure{Transitsinglejunction1}
  /*std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2
  >> dau1 >> dau2 "<< endl; for(int icheck1 = 0; icheck1 <
  Transitsinglejunction1.size(); icheck1++ ){ std::cout <<endl<<"SingleJunction
  : "<<icheck1<<endl; for(int i = 0; i <
  Transitsinglejunction1.at(icheck1).size(); i++){ vector<HHparton> temp =
  Transitsinglejunction1.at(icheck1); std::cout  << i <<"   "<< temp.at(i).id()
  <<"   "<< temp.at(i).PY_stat() <<"   "<< temp.at(i).PY_par1() <<"   "<<
  temp.at(i).PY_par2() <<"   "
      << temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<
  temp.at(i).col()<< " , " << temp.at(i).acol() << " ) " <<endl;
    }
  }*/

  // dignostic measure{Dijunction1}
  /*std::cout <<endl;
  std::cout <<" List of Tempsorting for Dijunction"<<endl;
  for(int icheck1 = 0; icheck1 < Tempsorting1.size(); icheck1++){
    std::cout <<" Sorted Leg : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Tempsorting1.at(icheck1).size(); icheck2++){
        std::cout <<" ( "<<Tempsorting1.at(icheck1).at(icheck2).col()<<" ,
  "<<Tempsorting1.at(icheck1).at(icheck2).acol()<< " )  ";
    }
    std::cout <<endl;
  }*/

  // So far, M,D tags for each junction and dijunction structure are designated.
  // and EP_conservation checking will be done through Transitdijunction1 and
  // Transitsinglejunction1 vectors
  for (int itd1 = 0; itd1 < Transitdijunction1.size(); itd1++) {
    /*JSINFO << "This is dijunction " << itd1+1 << " before:";
    for(int i = 0; i <  Transitdijunction1.at(itd1).size(); i++){
      vector<HHparton> temp = Transitdijunction1.at(itd1);
      std::cout  << i <<"   "<< temp.at(i).id() <<"   "<< temp.at(i).PY_stat()
    <<"   "<< temp.at(i).PY_par1() <<"   "<< temp.at(i).PY_par2() <<"   "
      << temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<
    temp.at(i).col()<< " , " << temp.at(i).acol() << " ) "
      << "  ( "<<  temp.at(i).e()<< " , " << temp.at(i).px()<< " , " <<
    temp.at(i).py()<< " , " << temp.at(i).pz() << " ) "
      << "  ( "<<  temp.at(i).t()<< " , " << temp.at(i).x()<< " , " <<
    temp.at(i).y()<< " , " << temp.at(i).z() << " ) " << std::endl;
    }*/
    bool EP_conserved = false;
    while (!EP_conserved) {
      EP_conserved = true;
      for (int i = 0; i < Transitdijunction1.at(itd1).size(); ++i) {
        if (Transitdijunction1[itd1][i].PY_stat() >= 0 &&
            std::abs(Transitdijunction1[itd1][i].id()) < 1000) {
          continue;
        }
        FourVector P_new(0., 0., 0., 0.);
        FourVector pos_new(0., 0., 0., 0.);
        int jmax = (Transitdijunction1[itd1][i].PY_dau2() >
                    Transitdijunction1[itd1][i].PY_dau1())
                       ? Transitdijunction1[itd1][i].PY_dau2()
                       : Transitdijunction1[itd1][i].PY_dau1();
        // JSINFO << "jmax = " << jmax;
        for (int j = Transitdijunction1[itd1][i].PY_dau1(); j < jmax + 1; ++j) {
          double n = double(j - Transitdijunction1[itd1][i].PY_dau1()) + 1.;
          // JSINFO << "n = " << n;
          // JSINFO << "transdijunct.x_t = "<<Transitdijunction1[itd1][j].x_t();
          // JSINFO <<
          // pos_new.t()+(Transitdijunction1[itd1][j].x_t()-pos_new.t())/n;
          P_new.Set(P_new.x() + Transitdijunction1[itd1][j].px(),
                    P_new.y() + Transitdijunction1[itd1][j].py(),
                    P_new.z() + Transitdijunction1[itd1][j].pz(),
                    P_new.t() + Transitdijunction1[itd1][j].e());
          pos_new.Set(
              pos_new.x() + (Transitdijunction1[itd1][j].x() - pos_new.x()) / n,
              pos_new.y() + (Transitdijunction1[itd1][j].y() - pos_new.y()) / n,
              pos_new.z() + (Transitdijunction1[itd1][j].z() - pos_new.z()) / n,
              pos_new.t() +
                  (Transitdijunction1[itd1][j].x_t() - pos_new.t()) / n);
        }
        // JSINFO << "p = (" << P_new.t() << ", " << P_new.y() << ", " <<
        // P_new.y() << ", " << P_new.z() << ")"; JSINFO << "x = (" <<
        // pos_new.t() << ", " << pos_new.y() << ", " << pos_new.y() << ", " <<
        // pos_new.z() << ")";
        if ((dif2(P_new, Transitdijunction1[itd1][i].P()) +
                 (P_new.t() - Transitdijunction1[itd1][i].e()) *
                     (P_new.t() - Transitdijunction1[itd1][i].e()) >
             0.00000001 /*0.0001^2*/) ||
            (dif2(pos_new, Transitdijunction1[itd1][i].pos()) +
                 (pos_new.t() - Transitdijunction1[itd1][i].x_t()) *
                     (pos_new.t() - Transitdijunction1[itd1][i].x_t()) >
             0.00000001)) {
          Transitdijunction1[itd1][i].P(P_new);
          Transitdijunction1[itd1][i].pos(pos_new);
          Transitdijunction1[itd1][i].mass(
              Transitdijunction1[itd1][i].e() *
                  Transitdijunction1[itd1][i].e() -
              Transitdijunction1[itd1][i].px() *
                  Transitdijunction1[itd1][i].px() -
              Transitdijunction1[itd1][i].py() *
                  Transitdijunction1[itd1][i].py() -
              Transitdijunction1[itd1][i].pz() *
                  Transitdijunction1[itd1][i].pz());
          Transitdijunction1[itd1][i].mass(
              (Transitdijunction1[itd1][i].mass() >= 0.)
                  ? sqrt(Transitdijunction1[itd1][i].mass())
                  : -sqrt(-Transitdijunction1[itd1][i]
                               .mass()));  // don't need mass >0 for reco now
          EP_conserved = false;
        }
      }
    }
    for (int i = 0; i < Transitdijunction1.at(itd1).size(); ++i) {
      if ((Transitdijunction1[itd1][i].PY_stat() ==
           -21)) {  // && (Transitdijunction1[itd1][i].PY_dau2() >
                    // Transitdijunction1[itd1][i].PY_dau1())){
        Transitdijunction1[itd1][i].px(Transitdijunction1[itd1][i].px() / 2.);
        Transitdijunction1[itd1][i].py(Transitdijunction1[itd1][i].py() / 2.);
        Transitdijunction1[itd1][i].pz(Transitdijunction1[itd1][i].pz() / 2.);
        Transitdijunction1[itd1][i].e(Transitdijunction1[itd1][i].e() / 2.);
        Transitdijunction1[itd1][i].mass(Transitdijunction1[itd1][i].mass() /
                                         2.);
      }
    }

    /*JSINFO << "This is dijunction " << itd1+1 << " after:";
    for(int i = 0; i <  Transitdijunction1.at(itd1).size(); i++){
      vector<HHparton> temp = Transitdijunction1.at(itd1);
      std::cout  << i <<"   "<< temp.at(i).id() <<"   "<< temp.at(i).PY_stat()
    <<"   "<< temp.at(i).PY_par1() <<"   "<< temp.at(i).PY_par2() <<"   "
      << temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<
    temp.at(i).col()<< " , " << temp.at(i).acol() << " ) "
      << "  ( "<<  temp.at(i).e()<< " , " << temp.at(i).px()<< " , " <<
    temp.at(i).py()<< " , " << temp.at(i).pz() << " ) "
      << "  ( "<<  temp.at(i).t()<< " , " << temp.at(i).x()<< " , " <<
    temp.at(i).y()<< " , " << temp.at(i).z() << " ) " << std::endl;
    }*/
  }

  // for transitsinglejunction1
  for (int its1 = 0; its1 < Transitsinglejunction1.size(); its1++) {
    bool EP_conserved = false;
    while (!EP_conserved) {
      EP_conserved = true;
      for (int i = 0; i < Transitsinglejunction1.at(its1).size(); ++i) {
        if (Transitsinglejunction1[its1][i].PY_stat() >= 0 &&
            std::abs(Transitsinglejunction1[its1][i].id()) < 1000) {
          continue;
        }
        FourVector P_new(0., 0., 0., 0.);
        FourVector pos_new(0., 0., 0., 0.);
        int jmax = (Transitsinglejunction1[its1][i].PY_dau2() >
                    Transitsinglejunction1[its1][i].PY_dau1())
                       ? Transitsinglejunction1[its1][i].PY_dau2()
                       : Transitsinglejunction1[its1][i].PY_dau1();
        for (int j = Transitsinglejunction1[its1][i].PY_dau1(); j < jmax + 1;
             ++j) {
          double n = double(j - Transitsinglejunction1[its1][i].PY_dau1()) + 1.;
          P_new.Set(P_new.x() + Transitsinglejunction1[its1][j].px(),
                    P_new.y() + Transitsinglejunction1[its1][j].py(),
                    P_new.z() + Transitsinglejunction1[its1][j].pz(),
                    P_new.t() + Transitsinglejunction1[its1][j].e());
          pos_new.Set(
              pos_new.x() +
                  (Transitsinglejunction1[its1][j].x() - pos_new.x()) / n,
              pos_new.y() +
                  (Transitsinglejunction1[its1][j].y() - pos_new.y()) / n,
              pos_new.z() +
                  (Transitsinglejunction1[its1][j].z() - pos_new.z()) / n,
              pos_new.t() +
                  (Transitsinglejunction1[its1][j].x_t() - pos_new.t()) / n);
        }
        if ((dif2(P_new, Transitsinglejunction1[its1][i].P()) +
                 (P_new.t() - Transitsinglejunction1[its1][i].e()) *
                     (P_new.t() - Transitsinglejunction1[its1][i].e()) >
             0.00000001 /*0.0001^2*/) ||
            (dif2(pos_new, Transitsinglejunction1[its1][i].pos()) +
                 (pos_new.t() - Transitsinglejunction1[its1][i].x_t()) *
                     (pos_new.t() - Transitsinglejunction1[its1][i].x_t()) >
             0.00000001)) {
          Transitsinglejunction1[its1][i].P(P_new);
          Transitsinglejunction1[its1][i].pos(pos_new);
          Transitsinglejunction1[its1][i].mass(
              Transitsinglejunction1[its1][i].e() *
                  Transitsinglejunction1[its1][i].e() -
              Transitsinglejunction1[its1][i].px() *
                  Transitsinglejunction1[its1][i].px() -
              Transitsinglejunction1[its1][i].py() *
                  Transitsinglejunction1[its1][i].py() -
              Transitsinglejunction1[its1][i].pz() *
                  Transitsinglejunction1[its1][i].pz());
          Transitsinglejunction1[its1][i].mass(
              (Transitsinglejunction1[its1][i].mass() >= 0.)
                  ? sqrt(Transitsinglejunction1[its1][i].mass())
                  : -sqrt(-Transitsinglejunction1[its1][i]
                               .mass()));  // don't need mass >0 for reco now
          EP_conserved = false;
        }
      }
    }
    for (int i = 0; i < Transitsinglejunction1.at(its1).size(); ++i) {
      if ((Transitsinglejunction1[its1][i].PY_stat() ==
           -21)) {  // && (Transitsinglejunction1[its1][i].PY_dau2() >
                    // Transitsinglejunction1[its1][i].PY_dau1())){
        Transitsinglejunction1[its1][i].px(
            Transitsinglejunction1[its1][i].px() / 2.);
        Transitsinglejunction1[its1][i].py(
            Transitsinglejunction1[its1][i].py() / 2.);
        Transitsinglejunction1[its1][i].pz(
            Transitsinglejunction1[its1][i].pz() / 2.);
        Transitsinglejunction1[its1][i].e(Transitsinglejunction1[its1][i].e() /
                                          2.);
        Transitsinglejunction1[its1][i].mass(
            Transitsinglejunction1[its1][i].mass() / 2.);
      }
    }
  }

  // make sure that there is enough energy in the system, such that the
  // dijunctions can be hadronized this is a rare case
  for (int itagdij1 = 0; itagdij1 < Transitdijunction1.size(); itagdij1++) {
    double baryon_mass_scaling = 1.25;
    double mass_baryons = Transitdijunction1.at(itagdij1).at(0).mass() +
                          Transitdijunction1.at(itagdij1).at(1).mass();
    double m1 = pythia.particleData.m0(
        std::abs(Transitdijunction1.at(itagdij1).at(0).id()));
    double m2 = pythia.particleData.m0(
        std::abs(Transitdijunction1.at(itagdij1).at(1).id()));
    if (mass_baryons < baryon_mass_scaling * (m1 + m2)) {
      // find fake gluon
      double new_gluon_mass = 2. * xmq;
      for (int itagdij2 = 0; itagdij2 < Transitdijunction1.at(itagdij1).size();
           itagdij2++) {
        HHparton p = Transitdijunction1.at(itagdij1).at(itagdij2);
        if (p.id() == 21) {
          new_gluon_mass =
              p.mass() + baryon_mass_scaling * (m1 + m2) - mass_baryons;
          p.mass(new_gluon_mass);
          double energy_new = std::sqrt(p.px() * p.px() + p.py() * p.py() +
                                        p.pz() * p.pz() + p.mass() * p.mass());
          p.e(energy_new);
          Transitdijunction1.at(itagdij1).at(itagdij2).mass(p.mass());
          Transitdijunction1.at(itagdij1).at(itagdij2).e(p.e());
        }
      }
      // find fake (anti-)quarks
      for (int itagdij2 = 0; itagdij2 < Transitdijunction1.at(itagdij1).size();
           itagdij2++) {
        HHparton p = Transitdijunction1.at(itagdij1).at(itagdij2);
        if (p.PY_stat() == -21 && std::abs(p.id()) < 6) {
          p.mass(new_gluon_mass / 2.);
          double energy_new = std::sqrt(p.px() * p.px() + p.py() * p.py() +
                                        p.pz() * p.pz() + p.mass() * p.mass());
          p.e(energy_new);
        }
        if (p.PY_stat() == -11 && std::abs(p.id()) > 1000) {
          p.mass(p.mass() +
                 (baryon_mass_scaling * (m1 + m2) - mass_baryons) / 2.);
          double energy_new = std::sqrt(p.px() * p.px() + p.py() * p.py() +
                                        p.pz() * p.pz() + p.mass() * p.mass());
          p.e(energy_new);
        }
      }
    }
  }

  for (int itagdij1 = 0; itagdij1 < Transitdijunction1.size(); itagdij1++) {
    for (int itagdij2 = 0; itagdij2 < Transitdijunction1.at(itagdij1).size();
         itagdij2++) {
      WaitingLineforPY.push_back(Transitdijunction1.at(itagdij1).at(itagdij2));
    }
  }

  // dignostic measure{Transitdijunction1}
  /*int icheck2 = 0;
  std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2
  >> dau1 >> dau2 "<< endl; for(int icheck1 = 0; icheck1 <
  Transitdijunction1.size(); icheck1++ ){ for(int i = 0; i <
  Transitdijunction1.at(icheck1).size(); i++){ vector<HHparton> temp =
  Transitdijunction1.at(icheck1); std::cout  << icheck2 <<"   "<<
  temp.at(i).id() <<"   "<< temp.at(i).PY_stat() <<"   "<< temp.at(i).PY_par1()
  <<"   "<< temp.at(i).PY_par2() <<"   "
      << temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<
  temp.at(i).col()<< " , " << temp.at(i).acol() << " ) " <<endl; icheck2++;
    }
  }*/

  // make sure that there is enough energy in the system, such that the
  // junctions can be hadronized this is a rare case
  for (int itagsj1 = 0; itagsj1 < Transitsinglejunction1.size(); itagsj1++) {
    double baryon_mass_scaling = 1.15;
    double mass_baryon = Transitsinglejunction1.at(itagsj1).at(0).mass();
    double m1 = pythia.particleData.m0(
        std::abs(Transitsinglejunction1.at(itagsj1).at(0).id()));

    if (mass_baryon < baryon_mass_scaling * m1) {
      double baryon_mass_correction = baryon_mass_scaling * m1 - mass_baryon;
      for (int itagsj2 = 0; itagsj2 < Transitsinglejunction1.at(itagsj1).size();
           itagsj2++) {
        if (std::abs(Transitsinglejunction1.at(itagsj1).at(itagsj2).id()) < 6) {
          HHparton p = Transitsinglejunction1.at(itagsj1).at(itagsj2);
          p.mass(p.mass() + baryon_mass_correction / 3.);
          double energy_new = std::sqrt(p.px() * p.px() + p.py() * p.py() +
                                        p.pz() * p.pz() + p.mass() * p.mass());
          p.e(energy_new);
          Transitsinglejunction1.at(itagsj1).at(itagsj2).mass(p.mass());
          Transitsinglejunction1.at(itagsj1).at(itagsj2).e(p.e());
        }
      }

      HHparton p = Transitsinglejunction1.at(itagsj1).at(0);
      p.mass(baryon_mass_scaling * m1);
      double energy_new = std::sqrt(p.px() * p.px() + p.py() * p.py() +
                                    p.pz() * p.pz() + p.mass() * p.mass());
      p.e(energy_new);
      Transitsinglejunction1.at(itagsj1).at(0).mass(p.mass());
      Transitsinglejunction1.at(itagsj1).at(0).e(p.e());
    }
  }

  for (int itagsj1 = 0; itagsj1 < Transitsinglejunction1.size(); itagsj1++) {
    for (int itagsj2 = 0; itagsj2 < Transitsinglejunction1.at(itagsj1).size();
         itagsj2++) {
      WaitingLineforPY.push_back(
          Transitsinglejunction1.at(itagsj1).at(itagsj2));
    }
  }

  /*for(int i=0; i<WaitingLineforPY.size(); i++) {
    std::cout << WaitingLineforPY[i].id() << "," << WaitingLineforPY[i].col() <<
  "," << WaitingLineforPY[i].acol() << "," << WaitingLineforPY[i].PY_par1() <<
  "," << WaitingLineforPY[i].PY_par2() << std::endl;
  }*/

  for (int ifins = 0; ifins < finalstring.size(); ifins++) {
    WaitingLineforPY.push_back(finalstring.at(ifins));
  }

  /*std::cout <<endl<<" Let's check final entry before PY!!"<<endl;
  std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2
  >> dau1 >> dau2 >> (col, acol) >> x_t"<< endl; for(int fincheck = 0; fincheck
  < WaitingLineforPY.size(); fincheck++){ vector<HHparton> temp =
  WaitingLineforPY; std::cout  << fincheck <<"   "<< temp.at(fincheck).id() <<"
  "<< temp.at(fincheck).PY_stat() <<"   "<< temp.at(fincheck).PY_par1() <<" "<<
  temp.at(fincheck).PY_par2() <<"   "
    << temp.at(fincheck).PY_dau1() <<"   "<< temp.at(fincheck).PY_dau2() << "  (
  "<<  temp.at(fincheck).col()<< " , " << temp.at(fincheck).acol() << " )   " <<
  temp.at(fincheck).x_t()<<endl;
  }*/

  for (int ifin = 0; ifin < WaitingLineforPY.size(); ifin++) {
    SP_prepremn.add(WaitingLineforPY[ifin]);
  }

  Tempjunctions.clear();  // clear these informations for next running
  JuncStructure.clear();
  realjuncindice.clear();
  IMStructure1.clear();
  Dijunction1.clear();
  DijunctionInfo1.clear();
  Recombearly1.clear();
  Tailoredstring1.clear();
  finalstring.clear();
  Transitdijunction1.clear();
  Transitsinglejunction1.clear();
  WaitingLineforPY.clear();
  Tempsorting1.clear();
}

// function to hand partons/strings and hadron resonances (and various other
// color neutral and colored objects) to Pythia8
bool HybridHadronization::invoke_py() {
  Event& event = pythia.event;

  // should have been checked before call, but if there are no partons/hadrons
  // to deal with, just exit without invoking pythia
  if (HH_pyremn.num() + HH_hadrons.num() == 0) {
    return true;
  }

  // first things first, need to reindex py_remn; just need to increment ALL by
  // one also restoring original id for color octet particles (or other 'odd'
  // colored particles)
  for (int i = 0; i < HH_pyremn.num(); ++i) {
    HH_pyremn[i].PY_par1(HH_pyremn[i].PY_par1() + 1);
    HH_pyremn[i].PY_par2(HH_pyremn[i].PY_par2() + 1);
    HH_pyremn[i].PY_dau1(HH_pyremn[i].PY_dau1() + 1);
    HH_pyremn[i].PY_dau2(HH_pyremn[i].PY_dau2() + 1);
    if (HH_pyremn[i].PY_origid() != 0) {
      HH_pyremn[i].id(HH_pyremn[i].PY_origid());
    }
  }

  bool need_hadronization = true;
  bool success = true;
  int attempt = 0;
  while (need_hadronization) {
    // incrementing attempt number
    ++attempt;

    // resetting PYTHIA event record, so that this event can be filled
    event.reset();
    HH_pythia_hadrons.clear();

    // number of partons/hadrons/particles handed to pythia
    int size_input = 0;

    // keeping track of partons from py_remn and hadrons into the event
    std::vector<int> eve_to_had;
    eve_to_had.push_back(0);

    // filling PYTHIA event record with the partons from this event
    int dijuncflag = 0;
    bool case1 = true;
    bool case2 = true;
    bool case3 = true;
    bool case4 = true;
    bool case5 = true;
    bool case6 = true;
    for (int i = 0; i < HH_pyremn.num(); ++i) {
      // code part to hadronize junctions and dijunctions separately
      // if fake pythia baryon (=junction or anti-junction) increment by one
      if (abs(HH_pyremn[i].id()) > 1112 && HH_pyremn[i].PY_tag1() != 0 &&
          HH_pyremn[i].PY_tag2() != 0 && HH_pyremn[i].PY_tag3() != 0) {
        dijuncflag++;
      }
      // this first part is for the hadronization of junctions/anti-junctions
      if ((dijuncflag == 1 && abs(HH_pyremn[i].id()) < 100 &&
           HH_pyremn[i].PY_par1() == 0 && HH_pyremn[i].PY_par2() == 0) ||
          (dijuncflag == 2 && abs(HH_pyremn[i - 1].id()) < 100)) {
        size_input = event.size() - 1;
        // event.listJunctions();
        // event.list();
        case1 &= pythia.next();
        // event.list();
        set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                         case1, false, false);
        event.reset();
        eve_to_had.clear();
        if (dijuncflag == 1) {
          dijuncflag = 0;
        }
        if (dijuncflag == 2) {
          dijuncflag = 1;
        }
      }
      // if we are past two junctions, check if they are back-to-back;
      if (dijuncflag == 2) {
        if (abs(HH_pyremn[i - 1].id()) > 1112) {
          dijuncflag = 3;
        }  // if yes, increment to 3; we are now inside a dijunction
        else {
          dijuncflag = 1;
        }  // no dijunction, must be a single one, step back to counter 1
      }
      // code red: we only reach this if we have reached the end of a
      // dijunction; it's either another fake baryon or something without mother
      // or daughter tags this second part is for the hadronization of
      // di-junctions
      if (dijuncflag == 4 ||
          (HH_pyremn[i].PY_par1() == 0 && HH_pyremn[i].PY_par2() == 0 &&
           HH_pyremn[i].PY_dau1() == 0 && HH_pyremn[i].PY_dau2() == 0 &&
           dijuncflag == 3)) {
        // hadronize event here (which is only one dijunction); need to copy all
        // the code that calls pythia AND refers to event.xxx up here
        size_input = event.size() - 1;
        // event.listJunctions();
        // event.list();
        case2 &= pythia.next();
        // event.list();
        set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                         case2, false, false);
        event.reset();
        eve_to_had.clear();
        if (dijuncflag == 4) {
          dijuncflag = 1;
        } else {
          dijuncflag = 0;
        }
      }

      // append( id, status, mother1, mother2, daughter1, daughter2, col, acol,
      // px, py, pz, e, m)
      event.append(HH_pyremn[i].id(), HH_pyremn[i].PY_stat(),
                   HH_pyremn[i].PY_par1(), HH_pyremn[i].PY_par2(),
                   HH_pyremn[i].PY_dau1(), HH_pyremn[i].PY_dau2(),
                   HH_pyremn[i].col(), HH_pyremn[i].acol(), HH_pyremn[i].px(),
                   HH_pyremn[i].py(), HH_pyremn[i].pz(), HH_pyremn[i].e(),
                   HH_pyremn[i].mass());
      // JSINFO <<
      // HH_pyremn[i].id()<<","<<HH_pyremn[i].PY_stat()<<","<<HH_pyremn[i].PY_par1()<<","<<HH_pyremn[i].PY_par2()<<","<<HH_pyremn[i].PY_dau1()<<","<<HH_pyremn[i].PY_dau2()<<","<<
      //   HH_pyremn[i].col()<<","<<HH_pyremn[i].acol()<<","<<HH_pyremn[i].px()<<","<<HH_pyremn[i].py()<<","<<HH_pyremn[i].pz()<<","<<HH_pyremn[i].e()<<","<<HH_pyremn[i].mass();
      event[event.size() - 1].vProd(HH_pyremn[i].x(), HH_pyremn[i].y(),
                                    HH_pyremn[i].z(), HH_pyremn[i].x_t());
      // for junction mother, adding junction to list manually for color tracing
      // appendJunction(int kind, int col0, int col1, int col2);
      if (std::abs(event[event.size() - 1].id()) > 1112) {
        event.appendJunction(((event[event.size() - 1].id() > 0) ? 1 : 2),
                             HH_pyremn[i].PY_tag1(), HH_pyremn[i].PY_tag2(),
                             HH_pyremn[i].PY_tag3());
      }
      eve_to_had.push_back(-i - 1);
    }

    size_input = event.size() - 1;
    // make PYTHIA hadronize this event, if it fails then retry N=10 times...
    // (PYTHIA can and will rarely fail, without concern) if this fails more
    // than 10 times, we may retry this event starting back before recombination
    // (some number of times)
    // event.listJunctions();
    // event.list();
    case3 &= pythia.next();
    // event.list();
    set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                     case3, false, false);
    event.reset();
    eve_to_had.clear();

    // if we want to decay particles in pythia we have to do it separately for
    // hadrons from recombination, otherwise the information about the origin is
    // lost adding in hadrons/leptons/other colorless particles too...
    if (reco_hadrons_pythia) {
      // reco hadrons
      for (int i = 0; i < HH_hadrons.num(); ++i) {
        if (HH_hadrons[i].is_final() && HH_hadrons[i].is_recohad() &&
            HH_hadrons[i].is_shsh()) {
          // to make sure mass is set appropriately
          double massnow = HH_hadrons[i].e() * HH_hadrons[i].e() -
                           (HH_hadrons[i].px() * HH_hadrons[i].px() +
                            HH_hadrons[i].py() * HH_hadrons[i].py() +
                            HH_hadrons[i].pz() * HH_hadrons[i].pz());
          massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
          event.append(HH_hadrons[i].id(), 81, 0, 0, HH_hadrons[i].px(),
                       HH_hadrons[i].py(), HH_hadrons[i].pz(),
                       HH_hadrons[i].e(), massnow);
          event[event.size() - 1].vProd(HH_hadrons[i].x(), HH_hadrons[i].y(),
                                        HH_hadrons[i].z(), HH_hadrons[i].x_t());
          eve_to_had.push_back(i + 1);
        }
      }
      case4 &= pythia.next();
      set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                       case4, true, true);
      event.reset();
      eve_to_had.clear();

      for (int i = 0; i < HH_hadrons.num(); ++i) {
        if (HH_hadrons[i].is_final() && HH_hadrons[i].is_recohad() &&
            HH_hadrons[i].is_shth()) {
          // to make sure mass is set appropriately
          double massnow = HH_hadrons[i].e() * HH_hadrons[i].e() -
                           (HH_hadrons[i].px() * HH_hadrons[i].px() +
                            HH_hadrons[i].py() * HH_hadrons[i].py() +
                            HH_hadrons[i].pz() * HH_hadrons[i].pz());
          massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
          event.append(HH_hadrons[i].id(), 81, 0, 0, HH_hadrons[i].px(),
                       HH_hadrons[i].py(), HH_hadrons[i].pz(),
                       HH_hadrons[i].e(), massnow);
          event[event.size() - 1].vProd(HH_hadrons[i].x(), HH_hadrons[i].y(),
                                        HH_hadrons[i].z(), HH_hadrons[i].x_t());
          eve_to_had.push_back(i + 1);
        }
      }
      case5 &= pythia.next();
      set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                       case5, true, false);
      event.reset();
      eve_to_had.clear();

      for (int i = 0; i < HH_hadrons.num(); ++i) {
        if (HH_hadrons[i].is_final() && !HH_hadrons[i].is_recohad()) {
          // to make sure mass is set appropriately
          double massnow = HH_hadrons[i].e() * HH_hadrons[i].e() -
                           (HH_hadrons[i].px() * HH_hadrons[i].px() +
                            HH_hadrons[i].py() * HH_hadrons[i].py() +
                            HH_hadrons[i].pz() * HH_hadrons[i].pz());
          massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
          event.append(HH_hadrons[i].id(), 81, 0, 0, HH_hadrons[i].px(),
                       HH_hadrons[i].py(), HH_hadrons[i].pz(),
                       HH_hadrons[i].e(), massnow);
          event[event.size() - 1].vProd(HH_hadrons[i].x(), HH_hadrons[i].y(),
                                        HH_hadrons[i].z(), HH_hadrons[i].x_t());
          eve_to_had.push_back(i + 1);
        }
      }
      case6 &= pythia.next();
      set_spacetime_for_pythia_hadrons(event, size_input, eve_to_had, attempt,
                                       case6, false, false);
      event.reset();
      eve_to_had.clear();
    }

    if (!case1 || !case2 || !case3 || !case4 || !case5 || !case6) {
      if (attempt > 4) {
        need_hadronization = false;
        success = false;
        break;
      }
      continue;
    }
    // this event has been successfully hadronized (hopefully), put all final
    // state particles into HH_pythia_hadrons; set is_final=true for these
    // this is done in set_spacetime_for_pythia_hadrons() and in DoHadronization
    // all these Hadrons are written to HH_hadrons
    need_hadronization = false;
  }
  return success;
}

void HybridHadronization::set_spacetime_for_pythia_hadrons(
    Pythia8::Event& event, int& size_input, std::vector<int>& eve_to_had,
    int pythia_attempt, bool find_positions, bool is_recohadron,
    bool recohadron_shsh) {
  // directly return if hadronization in pythia failed, such that nothing is
  // added to HH_pythia_hadrons
  if (!find_positions) {
    return;
  }

  vector<HHhadron>
      final_hadrons_from_pythia;  // hadrons with preliminary set positions
                                  // (average positions along string segment)

  vector<int> Case1_hadron_idx;   // index of the hadron from case 1 below
  vector<int> Case1_parton1_idx;  // index of the parton1 from case 1 below
  vector<int> Case1_parton2_idx;  // index of the parton2 from case 1 below

  vector<int> Case2_hadron_idx;   // index of the hadron from case 2 below
  vector<int> Case2_parton1_idx;  // index of the parton1 from case 2 below
  vector<int> Case2_parton2_idx;  // index of the parton2 from case 2 below
  vector<int> Case2_parton3_idx;  // index of the parton3 from case 2 below
  vector<double> Case2_junction_center;  // store the center position of the
                                         // junction in element 0=t,1=x,2=y,3=z
  // this works only for the current implementation, where the junctions are
  // handed over to pythia one after the other

  vector<int> Case3_hadron_idx;   // index of the hadron from case 3 below
  vector<int> Case3_parton1_idx;  // index of the parton1 from case 3 below
  vector<int> Case3_parton2_idx;  // index of the parton2 from case 3 below
  vector<double>
      Case3_partons_center;  // store the center position of the partons

  bool compute_more_precise_positions = true;
  bool warn_could_not_find_positions = false;

  for (int hadron_idx = 1; hadron_idx < event.size(); hadron_idx++) {
    if (event[hadron_idx].isFinal()) {  // take only final hadrons from pythia
      // convert pythia hadron to HHhadron
      HHhadron hadron_out;  // a number of other tags need to be set from the
                            // parent (either hadron or string)
      hadron_out.is_final(true);
      hadron_out.id(event[hadron_idx].id());
      hadron_out.mass(event[hadron_idx].m());
      hadron_out.px(event[hadron_idx].px());
      hadron_out.py(event[hadron_idx].py());
      hadron_out.pz(event[hadron_idx].pz());
      hadron_out.e(event[hadron_idx].e());

      // since using inbuilt pythia mother/daughter functions will segfault
      // 'occasionally', going to code it in by hand. this could probably be
      // done more efficiently, but it's good enough for now...
      std::vector<int> mothers;
      // using a stack system to fill mothers
      std::vector<int> stack;
      // filling stack with initial mothers
      if ((event[hadron_idx].mother1() < event[hadron_idx].mother2()) &&
          (event[hadron_idx].mother1() > 0) &&
          (std::abs(event[hadron_idx].status()) >= 81) &&
          (std::abs(event[hadron_idx].status()) <= 86)) {
        for (int hadron_parent = event[hadron_idx].mother1();
             hadron_parent <= event[hadron_idx].mother2(); ++hadron_parent) {
          stack.push_back(hadron_parent);
        }
      } else if ((event[hadron_idx].mother2() > 0) &&
                 (event[hadron_idx].mother1() != event[hadron_idx].mother2())) {
        stack.push_back(event[hadron_idx].mother1());
        stack.push_back(event[hadron_idx].mother2());
      } else if (event[hadron_idx].mother1() > 0) {
        stack.push_back(event[hadron_idx].mother1());
      } else {
        mothers.push_back(hadron_idx);
      }  // setting it as its own mother if there are no mothers (pythia didn't
         // decay a directly input hadron...)

      // filling the stack with any valid mothers of the 'current' stack element
      // then we check the 'current' stack element to see if it's a valid mother
      // (0<element<=n_input), if so, we write it to mothers
      while (stack.size() > 0) {
        int current = stack.back();
        stack.pop_back();
        if ((event[current].mother1() < event[current].mother2()) &&
            (event[current].mother1() > 0) &&
            (std::abs(event[current].status()) >= 81) &&
            (std::abs(event[current].status()) <= 86)) {
          for (int hadron_parent = event[current].mother1();
               hadron_parent <= event[current].mother2(); ++hadron_parent) {
            stack.push_back(hadron_parent);
          }
        } else if ((event[current].mother2() > 0) &&
                   (event[current].mother1() != event[current].mother2())) {
          stack.push_back(event[current].mother1());
          stack.push_back(event[current].mother2());
        } else if (event[current].mother1() > 0) {
          stack.push_back(event[current].mother1());
        }
        if ((current > 0) && (current <= size_input)) {
          mothers.push_back(current);
        }
      }

      // just in case...
      if (mothers.size() == 0) {
        mothers.push_back(hadron_idx);
      }

      // sorting and removing duplicate entries
      std::sort(mothers.begin(), mothers.end());
      mothers.erase(std::unique(mothers.begin(), mothers.end()), mothers.end());

      // first, using mothers to determine if this hadron was formed via
      // recombination, or by string fragmentation
      if (mothers[0] <= HH_pyremn.num()) {
        hadron_out.is_strhad(true);
      } else {
        hadron_out.is_recohad(true);
      }

      // lastly, using mothers (except fake) in original input to determine if
      // this is a shower-shower or shower-thermal hadron
      bool is_therm(false);
      for (int hadron_parent = 0; hadron_parent < mothers.size();
           ++hadron_parent) {
        if (mothers[hadron_parent] <= HH_pyremn.num()) {
          if (HH_pyremn[mothers[hadron_parent] - 1].orig() != -1) {
            hadron_out.add_par(mothers[hadron_parent] - 1);
            if (HH_pyremn[mothers[hadron_parent] - 1].is_thermal()) {
              is_therm = true;
            }
          }
        } else if (mothers[hadron_parent] <
                   size_input) {  // shouldn't actually need to check, but doing
                                  // so just in case.
          hadron_out.parh(eve_to_had[mothers[hadron_parent]] - 1);
          if (HH_hadrons[hadron_out.parh()].is_shth()) {
            is_therm = true;
          }
        }
      }
      if (is_therm) {
        hadron_out.is_shth(true);
      } else {
        hadron_out.is_shsh(true);
      }

      if (is_recohadron && recohadron_shsh) {
        hadron_out.is_recohad(true);
        hadron_out.is_shsh(true);
        hadron_out.is_shth(false);
      } else if (is_recohadron && !recohadron_shsh) {
        hadron_out.is_recohad(true);
        hadron_out.is_shsh(false);
        hadron_out.is_shth(true);
      }

      int hadron_col = event[hadron_idx].col();
      int hadron_acol = event[hadron_idx].acol();
      bool info_found = false;
      bool col_known = false;
      bool acol_known = false;

      // 3 main cases: col!=acol, col==acol!=0, col==acol==0
      if (hadron_col !=
          hadron_acol) {  // col!=acol -> this hadron should trace back to a
                          // single gluon; find it and grab its spacetime info
        for (int irem = 0; irem < HH_pyremn.num(); ++irem) {
          if (((HH_pyremn[irem].col() == hadron_col) &&
               (HH_pyremn[irem].acol() == hadron_acol)) ||
              ((HH_pyremn[irem].col() == hadron_acol) &&
               (HH_pyremn[irem].acol() ==
                hadron_col))) {  // alter for glu loops
            if (HH_pyremn[irem].PY_stat() <= 0) {
              continue;
            }
            hadron_out.x(HH_pyremn[irem].x());
            hadron_out.y(HH_pyremn[irem].y());
            hadron_out.z(HH_pyremn[irem].z());
            hadron_out.x_t(HH_pyremn[irem].x_t());
            info_found = true;
          }
        }
        if (!info_found) {
          // this is bad; could either be a hadron from *multiple* segments, or
          // still from a single gluon under color reconnections - OR BOTH! this
          // is a first-order handling to find all the intermediate partons
          // using the input partons - it really should be done using pythia's
          // history, but should suffice for most cases, for now... will apply
          // Dijkstra to the partons in the current string to find the shortest
          // path from hadron_col to hadron_acol

          // start by finding the partons with the col/acol tags
          int ptn1 = -1;
          int ptn2 = -1;
          for (int irem = 0; irem < HH_pyremn.num(); ++irem) {
            if ((HH_pyremn[irem].col() == hadron_col) &&
                (HH_pyremn[irem].PY_stat() > 0)) {
              ptn1 = irem;
            }
            if ((HH_pyremn[irem].acol() == hadron_acol) &&
                (HH_pyremn[irem].PY_stat() > 0)) {
              ptn2 = irem;
            }
            if ((ptn1 >= 0) && (ptn2 >= 0)) {
              break;
            }
          }

          // if ptn1 or ptn2 < 0, then we're giving up on this hadron (needs
          // pythia history to properly reconstruct the mother parton(s))
          // otherwise, we'll trace along the string from one parton end to the
          // other (using Dijkstra) to find all mother partons
          if ((ptn1 >= 0) && (ptn2 >= 0)) {
            // for 0th order approx, we'll just average the positions of the
            // known parton ends
            double pos_x, pos_y, pos_z, pos_t;
            pos_x = 0.;
            pos_y = 0.;
            pos_z = 0.;
            pos_t = 0.;
            pos_x += HH_pyremn[ptn1].x();
            pos_y += HH_pyremn[ptn1].y();
            pos_z += HH_pyremn[ptn1].z();
            pos_t += HH_pyremn[ptn1].x_t();
            pos_x += HH_pyremn[ptn2].x();
            pos_y += HH_pyremn[ptn2].y();
            pos_z += HH_pyremn[ptn2].z();
            pos_t += HH_pyremn[ptn2].x_t();
            pos_x /= 2.;
            pos_y /= 2.;
            pos_z /= 2.;
            pos_t /= 2.;
            hadron_out.x(pos_x);
            hadron_out.y(pos_y);
            hadron_out.z(pos_z);
            hadron_out.x_t(pos_t);
            info_found = true;
            Case1_hadron_idx.push_back(hadron_idx);
            Case1_parton1_idx.push_back(ptn1);
            Case1_parton2_idx.push_back(ptn2);
          }
        }
      } else if ((hadron_col == hadron_acol) &&
                 (hadron_col != 0)) {  // col==acol -> this hadron traces back
                                       // to a string segment; find the two
                                       // partons and interpolate position
        // there are 2 cases here - an easy case where the color tags match the
        // original partons, and a hard case where we need to trace history
        int ptn1 = 0;
        int ptn2 = 0;
        int ptn3 = 0;
        bool col_found = false;
        bool acol_found = false;
        bool ptn3_found = false;
        while (ptn1 < HH_pyremn.num()) {
          if ((HH_pyremn[ptn1].col() == hadron_col) &&
              (HH_pyremn[ptn1].PY_stat() > 0)) {
            col_found = true;
            break;
          }
          ++ptn1;
        }
        while (ptn2 < HH_pyremn.num()) {
          if ((HH_pyremn[ptn2].acol() == hadron_acol) &&
              (HH_pyremn[ptn2].PY_stat() > 0)) {
            acol_found = true;
            break;
          }
          ++ptn2;
        }

        col_known = col_found;
        acol_known = acol_found;
        // if we don't find color/anticolor tags in input partons, we'll need to
        // trace along pythia event/history to find either/both we can grab what
        // pythia claims are the mothers, then search that for the color tag,
        // which should return ONE candidate repeat this until event_i <=
        // HH_pyremn.num() - this will be our input would it be better to start
        // at the beginning, grab everything that doesn't have a terminal color
        // tag in a hadron, and trace those to respective hadron(s)?

        // since these *only* should happen in junction systems at the junction,
        // use the junction list to find the other 2 color tags, then find those
        // in orig. ptns!
        if (!col_found) {
          // find the junction with matching color, grab the other 2
          int coll[2] = {0, 0};
          for (int iJ = 0; iJ < pythia.event.sizeJunction(); ++iJ) {
            for (int iC = 0; iC < 3; ++iC) {
              if (pythia.event.colJunction(iJ, iC) == hadron_col) {
                if (iC == 0) {
                  coll[0] = pythia.event.colJunction(iJ, 1);
                  coll[1] = pythia.event.colJunction(iJ, 2);
                } else if (iC == 1) {
                  coll[0] = pythia.event.colJunction(iJ, 0);
                  coll[1] = pythia.event.colJunction(iJ, 2);
                } else {
                  coll[0] = pythia.event.colJunction(iJ, 0);
                  coll[1] = pythia.event.colJunction(iJ, 1);
                }
                col_found = true;
                break;
              }
            }
            if (col_found) {
              break;
            }
          }
          // since pythia refuses to keep track of initial junctions given to
          // it, we'll do it manually
          if (!col_found) {
            for (int irem = 0; irem < HH_pyremn.num(); ++irem) {
              if (std::abs(HH_pyremn[irem].id()) > 1112) {
                if (HH_pyremn[irem].PY_tag1() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag2();
                  coll[1] = HH_pyremn[irem].PY_tag3();
                  col_found = true;
                  break;
                } else if (HH_pyremn[irem].PY_tag2() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag1();
                  coll[1] = HH_pyremn[irem].PY_tag3();
                  col_found = true;
                  break;
                } else if (HH_pyremn[irem].PY_tag3() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag1();
                  coll[1] = HH_pyremn[irem].PY_tag2();
                  col_found = true;
                  break;
                }
              }
            }
          }
          // now we know the color tags needed, we look for them.
          if (col_found) {
            bool found_tags[2] = {0, 0};
            ptn1 = 0;
            while (ptn1 < HH_pyremn.num()) {
              if ((HH_pyremn[ptn1].acol() == coll[0]) &&
                  (HH_pyremn[ptn1].PY_stat() > 0)) {
                found_tags[0] = true;
                break;
              }
              ++ptn1;
            }
            while (ptn3 < HH_pyremn.num()) {
              if ((HH_pyremn[ptn3].acol() == coll[1]) &&
                  (HH_pyremn[ptn2].PY_stat() > 0)) {
                found_tags[1] = true;
                break;
              }
              ++ptn3;
            }
            if (found_tags[0] && found_tags[1]) {
              ptn3_found = true;
            }
          }
        }
        if (!acol_found) {
          // check to see if there's already 3 partons - if so, then we need a
          // warning/error!
          if (ptn3_found) {
            JSWARN << "A hadron was found with more than 3 partonic parents "
                      "for space-time info!";
          }
          // otherwise, works just like above...
          int coll[2] = {0, 0};
          for (int iJ = 0; iJ < pythia.event.sizeJunction(); ++iJ) {
            for (int iC = 0; iC < 3; ++iC) {
              if (pythia.event.colJunction(iJ, iC) == hadron_acol) {
                if (iC == 0) {
                  coll[0] = pythia.event.colJunction(iJ, 1);
                  coll[1] = pythia.event.colJunction(iJ, 2);
                } else if (iC == 1) {
                  coll[0] = pythia.event.colJunction(iJ, 0);
                  coll[1] = pythia.event.colJunction(iJ, 2);
                } else {
                  coll[0] = pythia.event.colJunction(iJ, 0);
                  coll[1] = pythia.event.colJunction(iJ, 1);
                }
                acol_found = true;
                break;
              }
            }
            if (acol_found) {
              break;
            }
          }
          // since pythia refuses to keep track of initial junctions given to
          // it, we'll do it manually
          if (!acol_found) {
            for (int irem = 0; irem < HH_pyremn.num(); ++irem) {
              if (std::abs(HH_pyremn[irem].id()) > 1112) {
                if (HH_pyremn[irem].PY_tag1() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag2();
                  coll[1] = HH_pyremn[irem].PY_tag3();
                  acol_found = true;
                  break;
                } else if (HH_pyremn[irem].PY_tag2() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag1();
                  coll[1] = HH_pyremn[irem].PY_tag3();
                  acol_found = true;
                  break;
                } else if (HH_pyremn[irem].PY_tag3() == hadron_col) {
                  coll[0] = HH_pyremn[irem].PY_tag1();
                  coll[1] = HH_pyremn[irem].PY_tag2();
                  acol_found = true;
                  break;
                }
              }
            }
          }
          // now we know the anti-color tags needed, we look for them.
          if (acol_found) {
            bool found_tags[2] = {0, 0};
            ptn2 = 0;
            while (ptn2 < HH_pyremn.num()) {
              if ((HH_pyremn[ptn2].col() == coll[0]) &&
                  (HH_pyremn[ptn1].PY_stat() > 0)) {
                found_tags[0] = true;
                break;
              }
              ++ptn2;
            }
            while (ptn3 < HH_pyremn.num()) {
              if ((HH_pyremn[ptn3].col() == coll[1]) &&
                  (HH_pyremn[ptn2].PY_stat() > 0)) {
                found_tags[1] = true;
                break;
              }
              ++ptn3;
            }
            if (found_tags[0] && found_tags[1]) {
              ptn3_found = true;
            }
          }
        }

        col_known = col_found;
        acol_known = acol_found;
        // turning off warning iff things work
        if (col_found && acol_found) {
          info_found = true;
        }

        // setting position from found partons
        double pos_x, pos_y, pos_z, pos_t;
        pos_x = 0.;
        pos_y = 0.;
        pos_z = 0.;
        pos_t = 0.;
        pos_x += HH_pyremn[ptn1].x();
        pos_y += HH_pyremn[ptn1].y();
        pos_z += HH_pyremn[ptn1].z();
        pos_t += HH_pyremn[ptn1].x_t();
        pos_x += HH_pyremn[ptn2].x();
        pos_y += HH_pyremn[ptn2].y();
        pos_z += HH_pyremn[ptn2].z();
        pos_t += HH_pyremn[ptn2].x_t();
        Case2_hadron_idx.push_back(hadron_idx);
        Case2_parton1_idx.push_back(ptn1);
        Case2_parton2_idx.push_back(ptn2);
        if (ptn3_found) {
          pos_x += HH_pyremn[ptn3].x();
          pos_y += HH_pyremn[ptn3].y();
          pos_z += HH_pyremn[ptn3].z();
          pos_t += HH_pyremn[ptn3].x_t();
          pos_x /= 3.;
          pos_y /= 3.;
          pos_z /= 3.;
          pos_t /= 3.;
          Case2_parton3_idx.push_back(ptn3);
          Case2_junction_center = {pos_t, pos_x, pos_y, pos_z};
        } else {
          pos_x /= 2.;
          pos_y /= 2.;
          pos_z /= 2.;
          pos_t /= 2.;
          Case2_parton3_idx.push_back(-1);
        }
        hadron_out.x(pos_x);
        hadron_out.y(pos_y);
        hadron_out.z(pos_z);
        hadron_out.x_t(pos_t);
      } else if ((hadron_col == hadron_acol) &&
                 (hadron_col == 0)) {  // col==acol==0 -> same as previous, but
                                       // pythia didn't give color tag since
                                       // this was from a short string q-qbar
        double avg_x, avg_y, avg_z, avg_t;
        avg_x = 0.;
        avg_y = 0.;
        avg_z = 0.;
        avg_t = 0.;
        int n_ptns = 0;  // should come out to 2 in the end, but keeping track
                         // just in case.
        bool found_first_q_or_qbar = false;
        bool found_second_q_or_qbar = false;
        for (int imot = 0; imot < mothers.size(); ++imot) {
          if (mothers[imot] <= HH_pyremn.num()) {
            info_found =
                true;  // if no mothers, then no position, and throw error!
            avg_x += HH_pyremn[mothers[imot] - 1].x();
            avg_y += HH_pyremn[mothers[imot] - 1].y();
            avg_z += HH_pyremn[mothers[imot] - 1].z();
            avg_t += HH_pyremn[mothers[imot] - 1].x_t();
            ++n_ptns;
            if (std::abs(HH_pyremn[mothers[imot] - 1].id()) < 7 &&
                !found_first_q_or_qbar) {
              Case3_parton1_idx.push_back(mothers[imot] - 1);
              found_first_q_or_qbar = true;
            } else if (std::abs(HH_pyremn[mothers[imot] - 1].id()) < 7 &&
                       !found_second_q_or_qbar && found_first_q_or_qbar) {
              Case3_parton2_idx.push_back(mothers[imot] - 1);
              found_second_q_or_qbar = true;
            }
          }
        }
        if (found_first_q_or_qbar && found_second_q_or_qbar) {
          Case3_hadron_idx.push_back(hadron_idx);
        }
        if (found_first_q_or_qbar && !found_second_q_or_qbar) {
          Case3_parton1_idx.pop_back();
        }
        if (n_ptns > 0) {
          avg_x /= double(n_ptns);
          avg_y /= double(n_ptns);
          avg_z /= double(n_ptns);
          avg_t /= double(n_ptns);
          Case3_partons_center = {avg_t, avg_x, avg_y, avg_z};
        }
        hadron_out.x(avg_x);
        hadron_out.y(avg_y);
        hadron_out.z(avg_z);
        hadron_out.x_t(avg_t);
      }
      if (!info_found) {
        compute_more_precise_positions = false;
        warn_could_not_find_positions = true;
        hadron_out.x(0.);
        hadron_out.y(0.);
        hadron_out.z(0.);
        hadron_out.x_t(0.);
      }

      // the mother procedure might skip some partons if there are junctions
      // involved this can be 'repaired' by taking a 'mother' parton, then
      // checking over all the partons in its string! (both adding to parents /
      // checking if thermal) this is done in hadronization calling function,
      // after invoke_py function is finished
      final_hadrons_from_pythia.push_back(hadron_out);
    }
  }

  // warn about not found position for hadrons only once per function call
  if (warn_could_not_find_positions &&
      GetXMLElementInt({"Afterburner", "include_fragmentation_hadrons"}) == 1) {
    VERBOSE(2) << "Could not find the spacetime information for hadron in "
                  "pythia event (string fragmentation attempt ="
               << pythia_attempt << "), set it to (0,0,0,0)";
  }

  // now let's do it a bit more precise and distribute the positions along the
  // string segments instead of using the average
  if (compute_more_precise_positions) {
    // find the more precise positions for case 1
    // have to find the hadrons from the same string segment
    vector<int> Case1_unique_parton1_idx = Case1_parton1_idx;
    std::sort(Case1_unique_parton1_idx.begin(), Case1_unique_parton1_idx.end());
    Case1_unique_parton1_idx.erase(std::unique(Case1_unique_parton1_idx.begin(),
                                               Case1_unique_parton1_idx.end()),
                                   Case1_unique_parton1_idx.end());
    // go through string segments
    for (int string_seg = 0; string_seg < Case1_unique_parton1_idx.size();
         string_seg++) {
      int segment_parton1_idx = Case1_unique_parton1_idx.at(string_seg);

      // select the hadrons from the string which need position modification
      vector<int> segment_hadron_idx;
      for (int iHad = 0; iHad < Case1_hadron_idx.size(); iHad++) {
        if (Case1_parton1_idx.at(iHad) == segment_parton1_idx) {
          segment_hadron_idx.push_back(iHad);
        }
      }

      // modify the hadron positions from the current string
      int number_hadrons_segment = segment_hadron_idx.size();
      if (number_hadrons_segment > 1) {
        for (int iHad = 0; iHad < number_hadrons_segment; iHad++) {
          int parton1_index = Case1_parton1_idx.at(segment_hadron_idx.at(iHad));
          int parton2_index = Case1_parton2_idx.at(segment_hadron_idx.at(iHad));
          int current_hadron_idx = segment_hadron_idx.at(iHad);

          // get the positions of the two partons at the ends of the string
          // segment
          double pos_x_ptn1 = HH_pyremn[parton1_index].x();
          double pos_y_ptn1 = HH_pyremn[parton1_index].y();
          double pos_z_ptn1 = HH_pyremn[parton1_index].z();
          double pos_t_ptn1 = HH_pyremn[parton1_index].x_t();
          double pos_x_ptn2 = HH_pyremn[parton2_index].x();
          double pos_y_ptn2 = HH_pyremn[parton2_index].y();
          double pos_z_ptn2 = HH_pyremn[parton2_index].z();
          double pos_t_ptn2 = HH_pyremn[parton2_index].x_t();

          double delta_x =
              (pos_x_ptn2 - pos_x_ptn1) / (number_hadrons_segment + 1);
          double delta_y =
              (pos_y_ptn2 - pos_y_ptn1) / (number_hadrons_segment + 1);
          double delta_z =
              (pos_z_ptn2 - pos_z_ptn1) / (number_hadrons_segment + 1);
          double delta_t =
              (pos_t_ptn2 - pos_t_ptn1) / (number_hadrons_segment + 1);

          double had_x = pos_x_ptn1 + (iHad + 1) * delta_x;
          double had_y = pos_y_ptn1 + (iHad + 1) * delta_y;
          double had_z = pos_z_ptn1 + (iHad + 1) * delta_z;
          double had_t = pos_t_ptn1 + (iHad + 1) * delta_t;

          // set the hadron to the new position along the string segment
          final_hadrons_from_pythia.at(current_hadron_idx).x(had_x);
          final_hadrons_from_pythia.at(current_hadron_idx).y(had_y);
          final_hadrons_from_pythia.at(current_hadron_idx).z(had_z);
          final_hadrons_from_pythia.at(current_hadron_idx).x_t(had_t);
        }
      }
    }

    // find the more precise positions for case 2
    // have to find the hadrons from the same string segment
    vector<int> Case2_unique_parton1_idx = Case2_parton1_idx;
    std::sort(Case2_unique_parton1_idx.begin(), Case2_unique_parton1_idx.end());
    Case2_unique_parton1_idx.erase(std::unique(Case2_unique_parton1_idx.begin(),
                                               Case2_unique_parton1_idx.end()),
                                   Case2_unique_parton1_idx.end());

    // go through string segments or junction legs
    for (int string_seg = 0; string_seg < Case2_unique_parton1_idx.size();
         string_seg++) {
      int segment_parton1_idx = Case2_unique_parton1_idx.at(string_seg);

      // select the hadrons from the string which need position modification
      vector<int> segment_hadron_idx;
      for (int iHad = 0; iHad < Case2_hadron_idx.size(); iHad++) {
        if (Case2_parton1_idx.at(iHad) == segment_parton1_idx) {
          segment_hadron_idx.push_back(iHad);
        }
      }

      // modify the hadron positions from the current string
      int number_hadrons_segment = segment_hadron_idx.size();
      if (number_hadrons_segment >
          1) {  // in case of a string segment use average position if only one
                // hadron, in case of junction use center (nothing to do here)
        for (int iHad = 0; iHad < number_hadrons_segment; iHad++) {
          int parton1_index = Case2_parton1_idx.at(segment_hadron_idx.at(iHad));
          int parton2_index = Case2_parton2_idx.at(segment_hadron_idx.at(iHad));
          int parton3_index = Case2_parton3_idx.at(segment_hadron_idx.at(iHad));
          int current_hadron_idx = segment_hadron_idx.at(iHad);

          if ((parton3_index == -1) &&
              (Case2_junction_center.size() ==
               0)) {  // this is the case where we have two partons at the
                      // string segment ends
            // get the positions of the two partons at the ends of the string
            // segment
            double pos_x_ptn1 = HH_pyremn[parton1_index].x();
            double pos_y_ptn1 = HH_pyremn[parton1_index].y();
            double pos_z_ptn1 = HH_pyremn[parton1_index].z();
            double pos_t_ptn1 = HH_pyremn[parton1_index].x_t();
            double pos_x_ptn2 = HH_pyremn[parton2_index].x();
            double pos_y_ptn2 = HH_pyremn[parton2_index].y();
            double pos_z_ptn2 = HH_pyremn[parton2_index].z();
            double pos_t_ptn2 = HH_pyremn[parton2_index].x_t();

            double delta_x =
                (pos_x_ptn2 - pos_x_ptn1) / (number_hadrons_segment + 1);
            double delta_y =
                (pos_y_ptn2 - pos_y_ptn1) / (number_hadrons_segment + 1);
            double delta_z =
                (pos_z_ptn2 - pos_z_ptn1) / (number_hadrons_segment + 1);
            double delta_t =
                (pos_t_ptn2 - pos_t_ptn1) / (number_hadrons_segment + 1);

            double had_x = pos_x_ptn1 + (iHad + 1) * delta_x;
            double had_y = pos_y_ptn1 + (iHad + 1) * delta_y;
            double had_z = pos_z_ptn1 + (iHad + 1) * delta_z;
            double had_t = pos_t_ptn1 + (iHad + 1) * delta_t;

            // set the hadron to the new position along the string segment
            final_hadrons_from_pythia.at(current_hadron_idx).x(had_x);
            final_hadrons_from_pythia.at(current_hadron_idx).y(had_y);
            final_hadrons_from_pythia.at(current_hadron_idx).z(had_z);
            final_hadrons_from_pythia.at(current_hadron_idx).x_t(had_t);
          } else if (parton3_index > -1 && Case2_junction_center.size() !=
                                               0) {  // this is the case with a
                                                     // junction / antijunction
            // compute the leg vectors (center is origin)
            vector<double> leg1_vec = {
                HH_pyremn[parton1_index].x_t() - Case2_junction_center[0],
                HH_pyremn[parton1_index].x() - Case2_junction_center[1],
                HH_pyremn[parton1_index].y() - Case2_junction_center[2],
                HH_pyremn[parton1_index].z() - Case2_junction_center[3]};
            vector<double> leg2_vec = {
                HH_pyremn[parton2_index].x_t() - Case2_junction_center[0],
                HH_pyremn[parton2_index].x() - Case2_junction_center[1],
                HH_pyremn[parton2_index].y() - Case2_junction_center[2],
                HH_pyremn[parton2_index].z() - Case2_junction_center[3]};
            vector<double> leg3_vec = {
                HH_pyremn[parton3_index].x_t() - Case2_junction_center[0],
                HH_pyremn[parton3_index].x() - Case2_junction_center[1],
                HH_pyremn[parton3_index].y() - Case2_junction_center[2],
                HH_pyremn[parton3_index].z() - Case2_junction_center[3]};
            double abs_leg1 = std::sqrt(
                leg1_vec[0] * leg1_vec[0] + leg1_vec[1] * leg1_vec[1] +
                leg1_vec[2] * leg1_vec[2] + leg1_vec[3] * leg1_vec[3]);
            double abs_leg2 = std::sqrt(
                leg2_vec[0] * leg2_vec[0] + leg2_vec[1] * leg2_vec[1] +
                leg2_vec[2] * leg2_vec[2] + leg2_vec[3] * leg2_vec[3]);
            double abs_leg3 = std::sqrt(
                leg3_vec[0] * leg3_vec[0] + leg3_vec[1] * leg3_vec[1] +
                leg3_vec[2] * leg3_vec[2] + leg3_vec[3] * leg3_vec[3]);
            if (abs_leg1 <= 1e-6 || abs_leg2 <= 1e-6 || abs_leg3 <= 1e-6) {
              continue;
            }
            double L = abs_leg1 + abs_leg2 + abs_leg3;
            double delta_l =
                L / (number_hadrons_segment + 1);  // inter hadron distance
            double distance = delta_l + iHad;  // center -> end leg1, center ->
                                               // end leg2, center -> end leg3
            double had_t, had_x, had_y, had_z;
            if (distance <= abs_leg1) {
              had_t = leg1_vec[0] * (distance / abs_leg1);
              had_x = leg1_vec[1] * (distance / abs_leg1);
              had_y = leg1_vec[2] * (distance / abs_leg1);
              had_z = leg1_vec[3] * (distance / abs_leg1);
            } else if ((abs_leg1 < distance) &&
                       (distance <= abs_leg1 + abs_leg2)) {
              had_t = leg2_vec[0] * ((distance - abs_leg1) / abs_leg2);
              had_x = leg2_vec[1] * ((distance - abs_leg1) / abs_leg2);
              had_y = leg2_vec[2] * ((distance - abs_leg1) / abs_leg2);
              had_z = leg2_vec[3] * ((distance - abs_leg1) / abs_leg2);
            } else {  // (abs_leg1+abs_leg2 < distance) && (distance <= L)
              had_t =
                  leg3_vec[0] * ((distance - abs_leg1 - abs_leg2) / abs_leg3);
              had_x =
                  leg3_vec[1] * ((distance - abs_leg1 - abs_leg2) / abs_leg3);
              had_y =
                  leg3_vec[2] * ((distance - abs_leg1 - abs_leg2) / abs_leg3);
              had_z =
                  leg3_vec[3] * ((distance - abs_leg1 - abs_leg2) / abs_leg3);
            }
            // set the hadron to the new position along the string segment
            final_hadrons_from_pythia.at(current_hadron_idx).x(had_x);
            final_hadrons_from_pythia.at(current_hadron_idx).y(had_y);
            final_hadrons_from_pythia.at(current_hadron_idx).z(had_z);
            final_hadrons_from_pythia.at(current_hadron_idx).x_t(had_t);
          }
        }
      }
    }

    // find the more precise positions for case 3
    // have to find the hadrons from the same string segment
    vector<int> Case3_unique_parton1_idx = Case3_parton1_idx;
    std::sort(Case3_unique_parton1_idx.begin(), Case3_unique_parton1_idx.end());
    Case3_unique_parton1_idx.erase(std::unique(Case3_unique_parton1_idx.begin(),
                                               Case3_unique_parton1_idx.end()),
                                   Case3_unique_parton1_idx.end());
    // go through string segments
    for (int string_seg = 0; string_seg < Case3_unique_parton1_idx.size();
         string_seg++) {
      int segment_parton1_idx = Case3_unique_parton1_idx.at(string_seg);

      // select the hadrons from the string which need position modification
      vector<int> segment_hadron_idx;
      for (int iHad = 0; iHad < Case3_hadron_idx.size(); iHad++) {
        if (Case3_parton1_idx.at(iHad) == segment_parton1_idx) {
          segment_hadron_idx.push_back(iHad);
        }
      }

      // modify the hadron positions from the current string
      int number_hadrons_segment = segment_hadron_idx.size();
      if (number_hadrons_segment > 1) {
        for (int iHad = 0; iHad < number_hadrons_segment; iHad++) {
          int parton1_index = Case3_parton1_idx.at(segment_hadron_idx.at(iHad));
          int parton2_index = Case3_parton2_idx.at(segment_hadron_idx.at(iHad));
          int current_hadron_idx = segment_hadron_idx.at(iHad);

          vector<double> seg1_vec = {
              HH_pyremn[parton1_index].x_t() - Case3_partons_center[0],
              HH_pyremn[parton1_index].x() - Case3_partons_center[1],
              HH_pyremn[parton1_index].y() - Case3_partons_center[2],
              HH_pyremn[parton1_index].z() - Case3_partons_center[3]};
          vector<double> seg2_vec = {
              HH_pyremn[parton2_index].x_t() - Case3_partons_center[0],
              HH_pyremn[parton2_index].x() - Case3_partons_center[1],
              HH_pyremn[parton2_index].y() - Case3_partons_center[2],
              HH_pyremn[parton2_index].z() - Case3_partons_center[3]};
          double abs_seg1 =
              std::sqrt(seg1_vec[0] * seg1_vec[0] + seg1_vec[1] * seg1_vec[1] +
                        seg1_vec[2] * seg1_vec[2] + seg1_vec[3] * seg1_vec[3]);
          double abs_seg2 =
              std::sqrt(seg2_vec[0] * seg2_vec[0] + seg2_vec[1] * seg2_vec[1] +
                        seg2_vec[2] * seg2_vec[2] + seg2_vec[3] * seg2_vec[3]);
          if (abs_seg1 <= 1e-6 || abs_seg2 <= 1e-6) {
            continue;
          }
          double L = abs_seg1 + abs_seg2;
          double delta_l =
              L / (number_hadrons_segment + 1);  // inter hadron distance
          double distance =
              delta_l + iHad;  // center -> end seg1, center -> end seg2
          double had_t, had_x, had_y, had_z;
          if (distance <= abs_seg1) {
            had_t = seg1_vec[0] * (distance / abs_seg1);
            had_x = seg1_vec[1] * (distance / abs_seg1);
            had_y = seg1_vec[2] * (distance / abs_seg1);
            had_z = seg1_vec[3] * (distance / abs_seg1);
          } else {  // (abs_leg1 < distance) && (distance <= abs_leg1+abs_leg2)
            had_t = seg2_vec[0] * ((distance - abs_seg1) / abs_seg2);
            had_x = seg2_vec[1] * ((distance - abs_seg1) / abs_seg2);
            had_y = seg2_vec[2] * ((distance - abs_seg1) / abs_seg2);
            had_z = seg2_vec[3] * ((distance - abs_seg1) / abs_seg2);
          }
          // set the hadron to the new position along the string segment
          final_hadrons_from_pythia.at(current_hadron_idx).x(had_x);
          final_hadrons_from_pythia.at(current_hadron_idx).y(had_y);
          final_hadrons_from_pythia.at(current_hadron_idx).z(had_z);
          final_hadrons_from_pythia.at(current_hadron_idx).x_t(had_t);
        }
      }
    }
  }

  // finally add the hadrons with position information to the HH_hadrons
  // the mother procedure might skip some partons if there are junctions
  // involved this can be 'repaired' by taking a 'mother' parton, then checking
  // over all the partons in its string! (both adding to parents / checking if
  // thermal)
  // this is done in hadronization calling function, after this invoke_py
  // function is finished
  for (int hadron_idx = 0; hadron_idx < final_hadrons_from_pythia.size();
       hadron_idx++) {
    HH_pythia_hadrons.add(final_hadrons_from_pythia[hadron_idx]);
  }
}

void HybridHadronization::bring_hadrons_to_mass_shell(
    hadron_collection& HH_hadrons) {
  // to calc. E conservation violation
  // double Ebefore = 0.; double Eafter = 0.; double pTbefore = 0.; double
  // pTafter = 0.; for(int iHad=0; iHad<HH_hadrons.num();
  // ++iHad){if(!HH_hadrons[iHad].is_final()){continue;} Ebefore +=
  // HH_hadrons[iHad].e(); pTbefore += HH_hadrons[iHad].pt();} std::cout <<
  // std::setprecision(5); std::ofstream eviol; std::ofstream echg;
  // std::ofstream ptchg; eviol.open("Evio_hadrons.dat", std::ios::app);
  // echg.open("Echg_hadrons.dat", std::ios::app);
  // ptchg.open("ptchg_hadrons.dat", std::ios::app);

  // hadron mass adjust
  double osf = 1e-6;
  for (int iHad = 0; iHad < HH_hadrons.num(); ++iHad) {
    // if this isn't a final state hadron, there's no point in fixing it.
    if (!HH_hadrons[iHad].is_final()) {
      continue;
    }

    // need hadron mass and pdg (pythia) mass to check
    double hadmass = HH_hadrons[iHad].mass();
    double m1 = pythia.particleData.m0(HH_hadrons[iHad].id());
    // if this hadron doesn't need to be fixed, then we skip it
    if (!(std::abs(hadmass - m1) / m1 > osf)) {
      continue;
    }

    // this hadron needs fixing, so we prepare to find a partner hadron to
    // adjust with it
    int partner = -1;

    // if this hadron has colors, then we can use that info to fix it
    if (HH_hadrons[iHad].cols.size() > 0) {
      // looking through the hadron list to find another hadron with the same
      // color tag.
      for (int jHad = 0; jHad < HH_hadrons.num(); ++jHad) {
        if (!HH_hadrons[jHad].is_final()) {
          continue;
        }  // do not want a non-final hadron
        if (iHad == jHad) {
          continue;
        }  // needs to be a different hadron
        double m2 = pythia.particleData.m0(HH_hadrons[jHad].id());
        double pair_mass =
            std::sqrt((HH_hadrons[iHad].e() + HH_hadrons[jHad].e()) *
                          (HH_hadrons[iHad].e() + HH_hadrons[jHad].e()) -
                      (HH_hadrons[iHad].px() + HH_hadrons[jHad].px()) *
                          (HH_hadrons[iHad].px() + HH_hadrons[jHad].px()) -
                      (HH_hadrons[iHad].py() + HH_hadrons[jHad].py()) *
                          (HH_hadrons[iHad].py() + HH_hadrons[jHad].py()) -
                      (HH_hadrons[iHad].pz() + HH_hadrons[jHad].pz()) *
                          (HH_hadrons[iHad].pz() + HH_hadrons[jHad].pz()));
        // variabele to compute the squared momentum difference of the partners
        // if this is too small the mass adjustment procedure would fail
        double momentum_diff =
            (HH_hadrons[iHad].px() - HH_hadrons[jHad].px()) *
                (HH_hadrons[iHad].px() - HH_hadrons[jHad].px()) +
            (HH_hadrons[iHad].py() - HH_hadrons[jHad].py()) *
                (HH_hadrons[iHad].py() - HH_hadrons[jHad].py()) +
            (HH_hadrons[iHad].pz() - HH_hadrons[jHad].pz()) *
                (HH_hadrons[iHad].pz() - HH_hadrons[jHad].pz());
        if (HH_hadrons[jHad].cols.size() >
            0) {  // if this hadron has color tags, check
          // JSINFO << "iHad = " << iHad << ", jHad = " << jHad << ",
          // momentum_diff = " << momentum_diff << ", pair_mass - m1 - m2 = " <<
          // pair_mass-m1-m2;
          for (int icol = 0; icol < HH_hadrons[iHad].cols.size(); ++icol) {
            for (int jcol = 0; jcol < HH_hadrons[jHad].cols.size(); ++jcol) {
              if (!(HH_hadrons[iHad].col(icol) == HH_hadrons[jHad].col(jcol)) ||
                  momentum_diff < 1e-6 || (pair_mass < m1 + m2)) {
                continue;
              }  // color not matching
              // if the previous condition isn't met, that means the colors
              // match!
              partner = jHad;
              break;
            }
            if (partner > -1) {
              break;
            }
          }
        }
        if (partner > -1) {
          break;
        }
      }
    }

    // ready to fix, if a partner has not been found, pick a close-by hadron
    // make sure that there is a relative momentum to reshuffle the masses
    if (partner == -1) {
      // grabbing closest final state hadron before and after, if it exists
      int iprev = iHad - 1;
      int inext = iHad + 1;
      while (iprev >= 0) {
        double m2 = pythia.particleData.m0(HH_hadrons[iprev].id());
        double pair_mass =
            std::sqrt((HH_hadrons[iHad].e() + HH_hadrons[iprev].e()) *
                          (HH_hadrons[iHad].e() + HH_hadrons[iprev].e()) -
                      (HH_hadrons[iHad].px() + HH_hadrons[iprev].px()) *
                          (HH_hadrons[iHad].px() + HH_hadrons[iprev].px()) -
                      (HH_hadrons[iHad].py() + HH_hadrons[iprev].py()) *
                          (HH_hadrons[iHad].py() + HH_hadrons[iprev].py()) -
                      (HH_hadrons[iHad].pz() + HH_hadrons[iprev].pz()) *
                          (HH_hadrons[iHad].pz() + HH_hadrons[iprev].pz()));
        double momentum_diff =
            (HH_hadrons[iHad].px() - HH_hadrons[iprev].px()) *
                (HH_hadrons[iHad].px() - HH_hadrons[iprev].px()) +
            (HH_hadrons[iHad].py() - HH_hadrons[iprev].py()) *
                (HH_hadrons[iHad].py() - HH_hadrons[iprev].py()) +
            (HH_hadrons[iHad].pz() - HH_hadrons[iprev].pz()) *
                (HH_hadrons[iHad].pz() - HH_hadrons[iprev].pz());
        if (HH_hadrons[iprev].is_final() && momentum_diff >= 1e-6 &&
            pair_mass >= m1 + m2) {
          break;
        }
        --iprev;
      }
      while (inext < HH_hadrons.num()) {
        double m2 = pythia.particleData.m0(HH_hadrons[inext].id());
        double pair_mass =
            std::sqrt((HH_hadrons[iHad].e() + HH_hadrons[inext].e()) *
                          (HH_hadrons[iHad].e() + HH_hadrons[inext].e()) -
                      (HH_hadrons[iHad].px() + HH_hadrons[inext].px()) *
                          (HH_hadrons[iHad].px() + HH_hadrons[inext].px()) -
                      (HH_hadrons[iHad].py() + HH_hadrons[inext].py()) *
                          (HH_hadrons[iHad].py() + HH_hadrons[inext].py()) -
                      (HH_hadrons[iHad].pz() + HH_hadrons[inext].pz()) *
                          (HH_hadrons[iHad].pz() + HH_hadrons[inext].pz()));
        double momentum_diff =
            (HH_hadrons[iHad].px() - HH_hadrons[inext].px()) *
                (HH_hadrons[iHad].px() - HH_hadrons[inext].px()) +
            (HH_hadrons[iHad].py() - HH_hadrons[inext].py()) *
                (HH_hadrons[iHad].py() - HH_hadrons[inext].py()) +
            (HH_hadrons[iHad].pz() - HH_hadrons[inext].pz()) *
                (HH_hadrons[iHad].pz() - HH_hadrons[inext].pz());
        if (HH_hadrons[inext].is_final() && momentum_diff >= 1e-6 &&
            pair_mass >= m1 + m2) {
          break;
        }
        ++inext;
      }

      // if there is no closest before, then grab closest after
      // or if no closest after, then grab closest before
      // if both exist, grab the one that's closest, unless both are the same in
      // which case grab the one after
      if (iprev < 0) {
        partner = inext;
      } else if (inext >= HH_hadrons.num()) {
        partner = iprev;
      } else {
        partner = (inext - iHad <= iHad - iprev) ? inext : iprev;
      }

      if ((partner >= HH_hadrons.num()) || (partner < 0)) {
        partner = -1;
      }
    }

    // choose the hadron with the smallest pair_mass - m1 - m2 value
    if (partner == -1) {
      double minimum = 1e6;
      int partner_temporary = -1;
      for (int jHad = 0; jHad < HH_hadrons.num(); ++jHad) {
        if (!HH_hadrons[jHad].is_final()) {
          continue;
        }  // do not want a non-final hadron
        if (iHad == jHad) {
          continue;
        }  // needs to be a different hadron
        double m2 = pythia.particleData.m0(HH_hadrons[jHad].id());
        double pair_mass =
            std::sqrt((HH_hadrons[iHad].e() + HH_hadrons[jHad].e()) *
                          (HH_hadrons[iHad].e() + HH_hadrons[jHad].e()) -
                      (HH_hadrons[iHad].px() + HH_hadrons[jHad].px()) *
                          (HH_hadrons[iHad].px() + HH_hadrons[jHad].px()) -
                      (HH_hadrons[iHad].py() + HH_hadrons[jHad].py()) *
                          (HH_hadrons[iHad].py() + HH_hadrons[jHad].py()) -
                      (HH_hadrons[iHad].pz() + HH_hadrons[jHad].pz()) *
                          (HH_hadrons[iHad].pz() + HH_hadrons[jHad].pz()));
        if (std::abs(pair_mass - m1 - m2) < minimum) {
          minimum = std::abs(pair_mass - m1 - m2);
          partner_temporary = jHad;
        }
      }
      if (partner_temporary != -1) {
        partner = partner_temporary;
      }
    }

    // by now, a partner *must* have been chosen - unless there's only 1 hadron
    // in the event, which is BAD. time to fix. if somehow everything failed,
    // just skip this hadron...
    if (partner == -1) {
      continue;
    }

    // std::cout << "H1_before: " << iHad << ", " << HH_hadrons[iHad].id() << ",
    // " << HH_hadrons[iHad].px() << ", " << HH_hadrons[iHad].py() << ", " <<
    // HH_hadrons[iHad].pz() << ", " << HH_hadrons[iHad].e() << "\n"; std::cout
    // << "H2_before: " << partner << ", " << HH_hadrons[partner].id() << ", "
    // << HH_hadrons[partner].px() << ", " << HH_hadrons[partner].py() << ", "
    // << HH_hadrons[partner].pz() << ", " << HH_hadrons[partner].e() << "\n";

    // Psys
    FourVector Psys;
    Psys.Set(HH_hadrons[iHad].px() + HH_hadrons[partner].px(),
             HH_hadrons[iHad].py() + HH_hadrons[partner].py(),
             HH_hadrons[iHad].pz() + HH_hadrons[partner].pz(),
             HH_hadrons[iHad].e() + HH_hadrons[partner].e());

    // CM velocity
    FourVector beta;
    beta.Set(Psys.x() / Psys.t(), Psys.y() / Psys.t(), Psys.z() / Psys.t(), 0.);
    beta.Set(beta.x(), beta.y(), beta.z(),
             1. / (sqrt(1. - (beta.x() * beta.x() + beta.y() * beta.y() +
                              beta.z() * beta.z()))));

    // boosting into CM frame
    FourVector p_CM[2];
    p_CM[0] = HH_hadrons[iHad].boost_P(beta);
    p_CM[1] = HH_hadrons[partner].boost_P(beta);

    // if E1 + E2 >= m1 + m2, shift momenta
    double m2 = pythia.particleData.m0(HH_hadrons[partner].id());
    double Etot = p_CM[0].t() + p_CM[1].t();
    if (Etot < m1 + m2 + 0.00001) {
      Etot = m1 + m2 +
             0.00001; /*eviol << Etot - (p_CM[0].t() + p_CM[1].t()) << "\n";*/
    }                 // can't shift, violating E/P cons.
    double E1 = Etot / 2. + ((m1 * m1) - (m2 * m2)) / (2. * Etot);
    double E2 = Etot / 2. - ((m1 * m1) - (m2 * m2)) / (2. * Etot);
    double pmag = sqrt(p_CM[0].x() * p_CM[0].x() + p_CM[0].y() * p_CM[0].y() +
                       p_CM[0].z() * p_CM[0].z());
    double fac = sqrt(Etot * Etot / 4. +
                      ((m1 * m1) - (m2 * m2)) * ((m1 * m1) - (m2 * m2)) /
                          (4. * Etot * Etot) -
                      ((m1 * m1) + (m2 * m2)) / 2.) /
                 pmag;

    // rescaling in CM frame
    p_CM[0].Set(fac * p_CM[0].x(), fac * p_CM[0].y(), fac * p_CM[0].z(), E1);
    p_CM[1].Set(fac * p_CM[1].x(), fac * p_CM[1].y(), fac * p_CM[1].z(), E2);

    // boosting back and setting hadron E,P to fixed values
    beta.Set(-beta.x(), -beta.y(), -beta.z(), 0.);
    beta.Set(beta.x(), beta.y(), beta.z(),
             1. / (sqrt(1. - (beta.x() * beta.x() + beta.y() * beta.y() +
                              beta.z() * beta.z()))));
    FourVector p_fin[2];
    p_fin[0] = HHboost(beta, p_CM[0]);
    p_fin[1] = HHboost(beta, p_CM[1]);
    HH_hadrons[iHad].P(p_fin[0]);
    HH_hadrons[partner].P(p_fin[1]);
    HH_hadrons[iHad].mass(m1);
    HH_hadrons[partner].mass(m2);

    // std::cout << "H1_after: " << iHad << ", " << HH_hadrons[iHad].id() << ",
    // " << HH_hadrons[iHad].px() << ", " << HH_hadrons[iHad].py() << ", " <<
    // HH_hadrons[iHad].pz() << ", " << HH_hadrons[iHad].e() << "\n"; std::cout
    // << "H2_after: " << partner << ", " << HH_hadrons[partner].id() << ", " <<
    // HH_hadrons[partner].px() << ", " << HH_hadrons[partner].py() << ", " <<
    // HH_hadrons[partner].pz() << ", " << HH_hadrons[partner].e() << "\n";
  }

  // calc. Etot after
  // for(int iHad=0; iHad<HH_hadrons.num();
  // ++iHad){if(!HH_hadrons[iHad].is_final()){continue;} Eafter +=
  // HH_hadrons[iHad].e(); pTafter += HH_hadrons[iHad].pt();} echg <<
  // Eafter-Ebefore << "\n"; ptchg << pTafter-pTbefore << "\n"; eviol.close();
  // echg.close(); ptchg.close();
}

void HybridHadronization::set_initial_parton_masses(
    parton_collection& HH_partons) {
  // to calc. E conservation violation
  /*double Ebefore = 0.; double Eafter = 0.; double pTbefore = 0.; double
  pTafter = 0.; for(int iPart=0; iPart<HH_partons.num(); ++iPart){Ebefore +=
  HH_partons[iPart].e(); pTbefore += HH_partons[iPart].pt();} std::cout <<
  std::setprecision(5); std::ofstream eviol; std::ofstream echg; std::ofstream
  ptchg; eviol.open("Evio_partons.dat", std::ios::app);
  echg.open("Echg_partons.dat", std::ios::app); ptchg.open("ptchg_partons.dat",
  std::ios::app);*/

  // case for only one parton in the event -> force to parton mass
  if (HH_partons.num() == 1) {
    HHparton parton = HH_partons[0];
    if (std::abs(parton.id()) < 6 &&
        parton.mass() <
            pythia.particleData.m0(std::abs(parton.id()))) {  // quark case
      parton.mass(pythia.particleData.m0(std::abs(parton.id())));
      parton.e(std::sqrt(parton.px() * parton.px() + parton.py() * parton.py() +
                         parton.pz() * parton.pz() +
                         parton.mass() * parton.mass()));
      HH_partons[0] = parton;
    } else if (parton.id() == 21 &&
               parton.mass() <
                   2. * xmq + 0.001) {  // gluon case (pythia has m_g = 0.)
      parton.mass(2. * xmq + 0.001);
      parton.e(std::sqrt(parton.px() * parton.px() + parton.py() * parton.py() +
                         parton.pz() * parton.pz() +
                         parton.mass() * parton.mass()));
      HH_partons[0] = parton;
    }
    return;
  }

  // parton mass adjust
  double osf = 1e-6;
  for (int iPart = 0; iPart < HH_partons.num(); ++iPart) {
    // if this isn't a final state parton, there's no point in fixing it.
    if ((HH_partons[iPart].id() != 21) &&
        (std::abs(HH_partons[iPart].id()) > 5)) {
      continue;
    }

    // need parton mass and pdg (pythia) mass to check
    double partmass = HH_partons[iPart].mass();
    double m1 = 0.;
    if (HH_partons[iPart].id() == 21) {
      m1 = 2. * xmq + 0.001;
    } else if (std::abs(HH_partons[iPart].id()) < 3) {
      m1 = xmq;
    } else if (std::abs(HH_partons[iPart].id()) == 3) {
      m1 = xms;
    } else if (std::abs(HH_partons[iPart].id()) == 4) {
      m1 = xmc;
    } else if (std::abs(HH_partons[iPart].id()) == 5) {
      m1 = xmb;
    }

    // if this parton doesn't need to be fixed, then we skip it
    if ((std::abs(partmass - m1) / m1 < osf) &&
        (std::abs(HH_partons[iPart].id()) < 6)) {
      continue;
    }

    if ((HH_partons[iPart].id() == 21) && (partmass > m1)) {
      continue;
    }

    // this parton needs fixing, so we prepare to find a partner to adjust with
    // it
    int partner = -1;

    // looking through the parton list to find another parton with the same
    // color tag.
    for (int jPart = 0; jPart < HH_partons.num(); ++jPart) {
      if (iPart == jPart || (std::abs(HH_partons[iPart].id()) > 5 &&
                             std::abs(HH_partons[iPart].id()) < 7)) {
        continue;
      }  // needs to be a different parton
      double m2 = 0.;
      if ((HH_partons[jPart].id() == 21) &&
          (HH_partons[jPart].mass() < 2. * xmq + 0.001)) {
        m2 = 2. * xmq + 0.001;
      } else if ((HH_partons[jPart].id() == 21) &&
                 (HH_partons[jPart].mass() >= 2. * xmq + 0.001)) {
        m2 = HH_partons[jPart].mass();
      } else if (std::abs(HH_partons[jPart].id()) < 3) {
        m2 = xmq;
      } else if (std::abs(HH_partons[jPart].id()) == 3) {
        m2 = xms;
      } else if (std::abs(HH_partons[jPart].id()) == 4) {
        m2 = xmc;
      } else if (std::abs(HH_partons[jPart].id()) == 5) {
        m2 = xmb;
      }
      double pair_mass =
          std::sqrt((HH_partons[iPart].e() + HH_partons[jPart].e()) *
                        (HH_partons[iPart].e() + HH_partons[jPart].e()) -
                    (HH_partons[iPart].px() + HH_partons[jPart].px()) *
                        (HH_partons[iPart].px() + HH_partons[jPart].px()) -
                    (HH_partons[iPart].py() + HH_partons[jPart].py()) *
                        (HH_partons[iPart].py() + HH_partons[jPart].py()) -
                    (HH_partons[iPart].pz() + HH_partons[jPart].pz()) *
                        (HH_partons[iPart].pz() + HH_partons[jPart].pz()));
      // variabele to compute the squared momentum difference of the partners
      // if this is too small the mass adjustment procedure would fail
      double momentum_diff =
          (HH_partons[iPart].px() - HH_partons[jPart].px()) *
              (HH_partons[iPart].px() - HH_partons[jPart].px()) +
          (HH_partons[iPart].py() - HH_partons[jPart].py()) *
              (HH_partons[iPart].py() - HH_partons[jPart].py()) +
          (HH_partons[iPart].pz() - HH_partons[jPart].pz()) *
              (HH_partons[iPart].pz() - HH_partons[jPart].pz());
      if (((HH_partons[iPart].col() == HH_partons[jPart].acol()) ||
           (HH_partons[iPart].acol() == HH_partons[jPart].col())) &&
          (momentum_diff > 1e-6) && (pair_mass >= m1 + m2)) {
        if ((HH_partons[jPart].id() == 21) &&
            (HH_partons[jPart].mass() < 2. * xmq + 0.001)) {
          partner = jPart;
          break;
        } else if ((std::abs(HH_partons[jPart].id()) < 3) &&
                   (std::abs(HH_partons[jPart].mass() - xmq) / xmq > osf)) {
          partner = jPart;
          break;
        } else if ((std::abs(HH_partons[jPart].id()) == 3) &&
                   (std::abs(HH_partons[jPart].mass() - xms) / xms > osf)) {
          partner = jPart;
          break;
        } else if ((std::abs(HH_partons[jPart].id()) == 4) &&
                   (std::abs(HH_partons[jPart].mass() - xmc) / xmc > osf)) {
          partner = jPart;
          break;
        } else if ((std::abs(HH_partons[jPart].id()) == 5) &&
                   (std::abs(HH_partons[jPart].mass() - xmb) / xmb > osf)) {
          partner = jPart;
          break;
        }
      }
    }

    // candidate for a parton with the smallest energy loss (if no partner is
    // found)
    int iPartnerCandidate = -1;
    double energy_loss_tracking = 1e6;  // just a large number

    // ready to fix, if a partner has not been found, pick a close-by parton
    // make sure that there is a relative momentum to reshuffle the masses
    if (partner == -1) {
      // grabbing closest parton before and after, if it exists
      int iprev = iPart - 1;
      int inext = iPart + 1;
      while (iprev >= 0) {
        double m2 = 0.;
        if ((HH_partons[iprev].id() == 21) &&
            (HH_partons[iprev].mass() < 2. * xmq + 0.001)) {
          m2 = 2. * xmq + 0.001;
        } else if ((HH_partons[iprev].id() == 21) &&
                   (HH_partons[iprev].mass() >= 2. * xmq + 0.001)) {
          m2 = HH_partons[iprev].mass();
        } else if (std::abs(HH_partons[iprev].id()) < 3) {
          m2 = xmq;
        } else if (std::abs(HH_partons[iprev].id()) == 3) {
          m2 = xms;
        } else if (std::abs(HH_partons[iprev].id()) == 4) {
          m2 = xmc;
        } else if (std::abs(HH_partons[iprev].id()) == 5) {
          m2 = xmb;
        }
        double pair_mass =
            std::sqrt((HH_partons[iPart].e() + HH_partons[iprev].e()) *
                          (HH_partons[iPart].e() + HH_partons[iprev].e()) -
                      (HH_partons[iPart].px() + HH_partons[iprev].px()) *
                          (HH_partons[iPart].px() + HH_partons[iprev].px()) -
                      (HH_partons[iPart].py() + HH_partons[iprev].py()) *
                          (HH_partons[iPart].py() + HH_partons[iprev].py()) -
                      (HH_partons[iPart].pz() + HH_partons[iprev].pz()) *
                          (HH_partons[iPart].pz() + HH_partons[iprev].pz()));
        double momentum_diff =
            (HH_partons[iPart].px() - HH_partons[iprev].px()) *
                (HH_partons[iPart].px() - HH_partons[iprev].px()) +
            (HH_partons[iPart].py() - HH_partons[iprev].py()) *
                (HH_partons[iPart].py() - HH_partons[iprev].py()) +
            (HH_partons[iPart].pz() - HH_partons[iprev].pz()) *
                (HH_partons[iPart].pz() - HH_partons[iprev].pz());
        if ((momentum_diff > 1e-6) && (pair_mass >= m1 + m2) &&
            (std::abs(HH_partons[iprev].id()) < 6 ||
             HH_partons[iprev].id() == 21)) {
          if (energy_loss_tracking > std::abs(pair_mass - m1 - m2)) {
            energy_loss_tracking = std::abs(pair_mass - m1 - m2);
            iPartnerCandidate = iprev;
          }
          break;
        }
        --iprev;
      }
      while (inext < HH_partons.num()) {
        double m2 = 0.;
        if ((HH_partons[inext].id() == 21) &&
            (HH_partons[inext].mass() < 2. * xmq + 0.001)) {
          m2 = 2. * xmq + 0.001;
        } else if ((HH_partons[inext].id() == 21) &&
                   (HH_partons[inext].mass() >= 2. * xmq + 0.001)) {
          m2 = HH_partons[inext].mass();
        } else if (std::abs(HH_partons[inext].id()) < 3) {
          m2 = xmq;
        } else if (std::abs(HH_partons[inext].id()) == 3) {
          m2 = xms;
        } else if (std::abs(HH_partons[inext].id()) == 4) {
          m2 = xmc;
        } else if (std::abs(HH_partons[inext].id()) == 5) {
          m2 = xmb;
        }
        double pair_mass =
            std::sqrt((HH_partons[iPart].e() + HH_partons[inext].e()) *
                          (HH_partons[iPart].e() + HH_partons[inext].e()) -
                      (HH_partons[iPart].px() + HH_partons[inext].px()) *
                          (HH_partons[iPart].px() + HH_partons[inext].px()) -
                      (HH_partons[iPart].py() + HH_partons[inext].py()) *
                          (HH_partons[iPart].py() + HH_partons[inext].py()) -
                      (HH_partons[iPart].pz() + HH_partons[inext].pz()) *
                          (HH_partons[iPart].pz() + HH_partons[inext].pz()));
        double momentum_diff =
            (HH_partons[iPart].px() - HH_partons[inext].px()) *
                (HH_partons[iPart].px() - HH_partons[inext].px()) +
            (HH_partons[iPart].py() - HH_partons[inext].py()) *
                (HH_partons[iPart].py() - HH_partons[inext].py()) +
            (HH_partons[iPart].pz() - HH_partons[inext].pz()) *
                (HH_partons[iPart].pz() - HH_partons[inext].pz());
        if ((momentum_diff > 1e-6) && (pair_mass >= m1 + m2) &&
            (std::abs(HH_partons[inext].id()) < 6 ||
             HH_partons[inext].id() == 21)) {
          if (energy_loss_tracking > std::abs(pair_mass - m1 - m2)) {
            energy_loss_tracking = std::abs(pair_mass - m1 - m2);
            iPartnerCandidate = inext;
          }
          break;
        }
        ++inext;
      }

      // if there is no closest before, then grab closest after
      // or if no closest after, then grab closest before
      // if both exist, grab the one that's closest, unless both are the same in
      // which case grab the one after
      if (iprev < 0) {
        partner = inext;
      } else if (inext >= HH_partons.num()) {
        partner = iprev;
      } else {
        partner = (inext - iPart <= iPart - iprev) ? inext : iprev;
      }

      if ((iPartnerCandidate >= HH_partons.num()) || (iPartnerCandidate < 0)) {
        iPartnerCandidate = -1;
      }
      if ((partner >= HH_partons.num()) || (partner < 0)) {
        partner = -1;
      }
    }

    // by now, a partner *must* have been chosen - unless there's only 1 parton
    // in the event, which is BAD. time to fix. if somehow everything failed,
    // just skip this parton and put it to its mass shell (here we have to
    // accept the energy conservation violation)
    if (partner == -1 && iPartnerCandidate >= 0) {
      partner = iPartnerCandidate;
    } else if (partner == -1 && iPartnerCandidate == -1) {
      HHparton parton = HH_partons[iPart];
      if (std::abs(parton.id()) < 6 &&
          parton.mass() <
              pythia.particleData.m0(std::abs(parton.id()))) {  // quark case
        parton.mass(pythia.particleData.m0(std::abs(parton.id())));
        parton.e(std::sqrt(
            parton.px() * parton.px() + parton.py() * parton.py() +
            parton.pz() * parton.pz() + parton.mass() * parton.mass()));
        HH_partons[iPart] = parton;
      } else if (parton.id() == 21 &&
                 parton.mass() <
                     2. * xmq + 0.001) {  // gluon case (pythia has m_g = 0.)
        parton.mass(2. * xmq + 0.001);
        parton.e(std::sqrt(
            parton.px() * parton.px() + parton.py() * parton.py() +
            parton.pz() * parton.pz() + parton.mass() * parton.mass()));
        HH_partons[iPart] = parton;
      }
      continue;
    }

    // std::cout << "P1_before: " << iPart << ", " << HH_partons[iPart].id() <<
    // ", " << HH_partons[iPart].px() << ", " << HH_partons[iPart].py() << ", "
    // << HH_partons[iPart].pz() << ", " << HH_partons[iPart].e() << ", " <<
    // HH_partons[iPart].mass() << "\n"; std::cout << "P2_before: " << partner
    // << ", " << HH_partons[partner].id() << ", " << HH_partons[partner].px()
    // << ", " << HH_partons[partner].py() << ", " << HH_partons[partner].pz()
    // << ", " << HH_partons[partner].e() << ", " << HH_partons[partner].mass()
    // << "\n";

    // Psys
    FourVector Psys;
    Psys.Set(HH_partons[iPart].px() + HH_partons[partner].px(),
             HH_partons[iPart].py() + HH_partons[partner].py(),
             HH_partons[iPart].pz() + HH_partons[partner].pz(),
             HH_partons[iPart].e() + HH_partons[partner].e());

    // CM velocity
    FourVector beta;
    beta.Set(Psys.x() / Psys.t(), Psys.y() / Psys.t(), Psys.z() / Psys.t(), 0.);
    beta.Set(beta.x(), beta.y(), beta.z(),
             1. / (sqrt(1. - (beta.x() * beta.x() + beta.y() * beta.y() +
                              beta.z() * beta.z()))));

    // boosting into CM frame
    FourVector p_CM[2];
    p_CM[0] = HH_partons[iPart].boost_P(beta);
    p_CM[1] = HH_partons[partner].boost_P(beta);

    // if E1 + E2 >= m1 + m2, shift momenta
    double m2 = 0.;
    if ((HH_partons[partner].id() == 21) &&
        (HH_partons[partner].mass() < 2. * xmq + 0.001)) {
      m2 = 2. * xmq + 0.001;
    } else if ((HH_partons[partner].id() == 21) &&
               (HH_partons[partner].mass() >= 2. * xmq + 0.001)) {
      m2 = HH_partons[partner].mass();
    } else if (std::abs(HH_partons[partner].id()) < 3) {
      m2 = xmq;
    } else if (std::abs(HH_partons[partner].id()) == 3) {
      m2 = xms;
    } else if (std::abs(HH_partons[partner].id()) == 4) {
      m2 = xmc;
    } else if (std::abs(HH_partons[partner].id()) == 5) {
      m2 = xmb;
    }
    double Etot = p_CM[0].t() + p_CM[1].t();
    if (Etot < m1 + m2 + 0.00001) {
      Etot = m1 + m2 +
             0.00001; /*eviol << Etot - (p_CM[0].t() + p_CM[1].t()) << "\n";*/
    }                 // can't shift, violating E/P cons.
    double E1 = Etot / 2. + ((m1 * m1) - (m2 * m2)) / (2. * Etot);
    double E2 = Etot / 2. - ((m1 * m1) - (m2 * m2)) / (2. * Etot);
    double pmag = sqrt(p_CM[0].x() * p_CM[0].x() + p_CM[0].y() * p_CM[0].y() +
                       p_CM[0].z() * p_CM[0].z());
    double fac = sqrt(Etot * Etot / 4. +
                      ((m1 * m1) - (m2 * m2)) * ((m1 * m1) - (m2 * m2)) /
                          (4. * Etot * Etot) -
                      ((m1 * m1) + (m2 * m2)) / 2.) /
                 pmag;
    // std::cout << "Etot_orig = " << p_CM[0].t() + p_CM[1].t() << " Etot = " <<
    // Etot << " pmag = " << pmag << " fac = " << fac << std::endl;

    // rescaling in CM frame
    p_CM[0].Set(fac * p_CM[0].x(), fac * p_CM[0].y(), fac * p_CM[0].z(), E1);
    p_CM[1].Set(fac * p_CM[1].x(), fac * p_CM[1].y(), fac * p_CM[1].z(), E2);

    // boosting back and setting parton E,P to fixed values
    beta.Set(-beta.x(), -beta.y(), -beta.z(), 0.);
    beta.Set(beta.x(), beta.y(), beta.z(),
             1. / (sqrt(1. - (beta.x() * beta.x() + beta.y() * beta.y() +
                              beta.z() * beta.z()))));
    FourVector p_fin[2];
    p_fin[0] = HHboost(beta, p_CM[0]);
    p_fin[1] = HHboost(beta, p_CM[1]);
    HH_partons[iPart].P(p_fin[0]);
    HH_partons[partner].P(p_fin[1]);
    HH_partons[iPart].mass(m1);
    HH_partons[partner].mass(m2);

    // std::cout << "P1_after: " << iPart << ", " << HH_partons[iPart].id() <<
    // ", " << HH_partons[iPart].px() << ", " << HH_partons[iPart].py() << ", "
    // << HH_partons[iPart].pz() << ", " << HH_partons[iPart].e() << ", " <<
    // HH_partons[iPart].mass() << "\n"; std::cout << "P2_after: " << partner <<
    // ", " << HH_partons[partner].id() << ", " << HH_partons[partner].px() <<
    // ", " << HH_partons[partner].py() << ", " << HH_partons[partner].pz() <<
    // ", " << HH_partons[partner].e() << ", " << HH_partons[partner].mass() <<
    // "\n";
  }

  // calc. Etot after
  /*for(int iPart=0; iPart<HH_partons.num(); ++iPart){Eafter +=
  HH_partons[iPart].e(); pTafter += HH_partons[iPart].pt();} echg <<
  Eafter-Ebefore << "\n"; ptchg << pTafter-pTbefore << "\n"; eviol.close();
  echg.close(); ptchg.close();*/
}
