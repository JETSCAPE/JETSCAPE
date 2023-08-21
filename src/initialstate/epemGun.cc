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

#include "epemGun.h"
#include <sstream>

#define MAGENTA "\033[35m"

using namespace std;

// Register the module with the base class
RegisterJetScapeModule<epemGun> epemGun::reg("epemGun");

epemGun::~epemGun() { VERBOSE(8); }

void epemGun::InitTask() {

  JSDEBUG << "Initialize epemGun";
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = off");
  readString("Init:showChangedParticleData = off");
  if (JetScapeLogger::Instance()->GetInfo()) {
    readString("Init:showProcesses = on");
    readString("Init:showChangedSettings = on");
    readString("Init:showMultipartonInteractions = on");
    readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  // Standard settings
  //readString(
      //"HardQCD:all = on"); // will repeat this line in the xml for demonstration
  //  readString("HardQCD:gg2ccbar = on"); // switch on heavy quark channel
  //readString("HardQCD:qqbar2ccbar = on");
  //readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  //readString("PartonLevel:ISR = on");
  //readString("PartonLevel:MPI = on");
  readString("PartonLevel:FSR = off");
  //readString("PromptPhoton:all=on");
  readString("PDF:lepton = off");                           //<-- Added by me.
  readString("WeakSingleBoson:ffbar2gmz=on"); //Scattering f fbar → gamma^*/Z^0, with full interference between the gamma^* and Z^0
  readString("WeakDoubleBoson:all=on"); //Common switch for the group of pair production of gamma^*/Z^0 and W^+-
  readString("WeakBosonExchange:all=on"); //Common switch for the group of gamma^*/Z^0 or W^+- exchange between two fermions

  //Stuff I added. Ask if I'm allowed to just do this.
  readString("23:onMode = off");
  readString("23:onIfAny = 1 2 3 4 5");
  //      //readString("ParticleDecays:limitTau0 = on");
  //        //readString("ParticleDecays:tau0Max = 299.792458");
  //          //readString("ColourReconnection:reconnect = off");

  // For parsing text
  stringstream numbf(stringstream::app | stringstream::in | stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);
  numbf.setf(ios::showpoint);
  numbf.precision(1);
  stringstream numbi(stringstream::app | stringstream::in | stringstream::out);

  std::string s = GetXMLElementText({"Hard", "epemGun", "name"});
  SetId(s);
  // cout << s << endl;

  // SC: read flag for FSR
  //FSR_on = GetXMLElementInt({"Hard", "epemGun", "FSR_on"});
  //if (FSR_on)
  //  readString("PartonLevel:FSR = on");
  //else
  //  readString("PartonLevel:FSR = off");

  //pTHatMin = GetXMLElementDouble({"Hard", "epemGun", "pTHatMin"});
  //pTHatMax = GetXMLElementDouble({"Hard", "epemGun", "pTHatMax"});

  //JSINFO << MAGENTA << "epem Gun with FSR_on: " << FSR_on;
  //JSINFO << MAGENTA << "epem Gun with " << pTHatMin << " < pTHat < "
  //       << pTHatMax;

  //numbf.str("PhaseSpace:pTHatMin = ");
  //numbf << pTHatMin;
  //readString(numbf.str());
  //numbf.str("PhaseSpace:pTHatMax = ");
  //numbf << pTHatMax;
  //readString(numbf.str());

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if (RandomXmlDescription) {
    tinyxml2::XMLElement *xmle =
        RandomXmlDescription->FirstChildElement("seed");
    if (!xmle)
      throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    JSWARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) << "Seeding pythia to " << seed;
  numbi << seed;
  readString(numbi.str());

  // Species
  readString("Beams:idA = 11");
  readString("Beams:idB = -11");

  // Energy
  eCM = GetXMLElementDouble({"Hard", "epemGun", "eCM"});
  numbf.str("Beams:eCM = ");
  numbf << eCM;
  readString(numbf.str());

  std::stringstream lines;
  lines << GetXMLElementText({"Hard", "epemGun", "LinesToRead"}, false);
  int i = 0;
  while (std::getline(lines, s, '\n')) {
    if (s.find_first_not_of(" \t\v\f\r") == s.npos)
      continue; // skip empty lines
    VERBOSE(7) << "Also reading in: " << s;
    readString(s);
  }

  // And initialize
  if (!init()) { // Pythia>8.1
    throw std::runtime_error("Pythia init() failed.");
  }

}

void epemGun::Exec() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  //Reading vir_factor from xml for MATTER
  //double vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

  bool flag62 = false;
  vector<Pythia8::Particle> p62;

  // sort by pt
  struct greater_than_pt {
    inline bool operator()(const Pythia8::Particle &p1,
                           const Pythia8::Particle &p2) {
      return (p1.pT() > p2.pT());
    }
  };

  do {
    next();
    //event.list();
    p62.clear();

    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for (int parid = 0; parid < event.size(); parid++) {
      if (parid < 3)
        continue; // 0, 1, 2: total event and beams
      Pythia8::Particle &particle = event[parid];

      //replacing diquarks with antiquarks (and anti-dq's with quarks)
      //the id is set to the heaviest quark in the diquark (except down quark)
      //this technically violates baryon number conservation over the entire event
      //also can violate electric charge conservation
      if( (std::abs(particle.id()) > 1100) && (std::abs(particle.id()) < 6000) && ((std::abs(particle.id())/10)%10 == 0) ){
        //if(particle.id() > 0){particle.id( -1*particle.id()/1000 );}
        //else{particle.id( particle.id()/1000 );}
        particle.id( -1*particle.id()/1000 );		//Changed from previous two lines.
      }

      //if (!FSR_on) {
        // only accept particles after MPI
        //if (particle.status() != 62) {continue;}
	if ( particle.status()<0 ) {continue;}
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 6 && (particle.id() != 21 && particle.id() != 22)) {
          continue;
	}

        // reject rare cases of very soft particles that don't have enough e to get
        // reasonable virtuality
        //if (particle.pT() < 1.0 / sqrt(vir_factor))
        //  continue;

        //if(particle.id()==22) cout<<"########this is a photon!######" <<endl;
        // accept
      /*} else { // FSR_on true: use Pythia vacuum shower instead of MATTER
        if (!particle.isFinal())
          continue;
        // only accept gluons and quarks
        // Also accept Gammas to put into the hadron's list
        if (fabs(particle.id()) > 5 &&
            (particle.id() != 21 && particle.id() != 22))
          continue;
      }*/
      p62.push_back(particle);
    }

    // if you want at least 2
    if (p62.size() < 2)
      continue;
    //if ( p62.size() < 1 ) continue;

    // Now have all candidates, sort them
    // sort by pt
    std::sort(p62.begin(), p62.end(), greater_than_pt());
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;

    flag62 = true;

  } while (!flag62);

  double p[4], xLoc[4];

  // This location should come from an initial state
  for (int i = 0; i <= 3; i++) {
    xLoc[i] = 0.0;
  };

  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

 /* if (!ini) {
    VERBOSE(1) << "No initial state module, setting the starting location to "
                  "0. Make sure to add e.g. trento before epemGun.";
  } else {
    auto num_bin_coll = ini->GetNumOfBinaryCollisions();
    if (num_bin_coll.size() == 0) {
      JSWARN << "num_of_binary_collisions is empty, setting the starting "
                "location to 0. Make sure to add e.g. trento before epemGun.";
    } else {
      std::discrete_distribution<> dist(
          begin(num_bin_coll), end(num_bin_coll)); // Create the distribution

      // Now generate values
      auto idx = dist(*GetMt19937Generator());
      auto coord = ini->CoordFromIdx(idx);
      xLoc[1] = get<0>(coord);
      xLoc[2] = get<1>(coord);
    }
  } */

  // Loop through particles

  // Only top two
  //for(int np = 0; np<2; ++np){

  // Accept them all

  int hCounter = 0;
  for (int np = 0; np < p62.size(); ++np) {
    Pythia8::Particle &particle = p62.at(np);

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << ", pT = " << particle.pT() << ", y = " << particle.y()
               << ", phi = " << particle.phi() << ", e = " << particle.e();

    VERBOSE(7) << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    auto ptn = make_shared<Parton>(0, particle.id(), 0, particle.pT(), particle.y(), particle.phi(), particle.e(), xLoc);
    ptn->set_color(particle.col());
    ptn->set_anti_color(particle.acol());
    ptn->set_max_color(1000 * (np + 1));
    AddParton(ptn);
  }

  //event.list();

  VERBOSE(8) << GetNHardPartons();
}
