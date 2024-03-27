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

// Create e+e- -> qqbar processes with Pythia and return the two inital partons with set virtualities

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
  readString("HadronLevel:all = off");
  //readString("PartonLevel:ISR = on");
  //readString("PartonLevel:MPI = on");
  readString("PartonLevel:FSR = off");
  //readString("PromptPhoton:all=on");
  readString("PDF:lepton = off");
  readString("WeakSingleBoson:ffbar2gmz=on"); //Scattering f fbar → gamma^*/Z^0, with full interference between the gamma^* and Z^0
  readString("WeakDoubleBoson:all=on"); //Common switch for the group of pair production of gamma^*/Z^0 and W^+-
  readString("WeakBosonExchange:all=on"); //Common switch for the group of gamma^*/Z^0 or W^+- exchange between two fermions

  //Stuff I added. Ask if I'm allowed to just do this.
  readString("23:onMode = off");
  readString("23:onIfAny = 1 2 3 4 5");

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

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double>{0.0, 1.0};

}

void epemGun::Exec() {
  VERBOSE(1) << "Run Hard Process : " << GetId() << " ...";
  VERBOSE(8) << "Current Event #" << GetCurrentEvent();
  //Reading vir_factor from xml for MATTER
  double vir_factor = GetXMLElementDouble({"Eloss", "Matter", "vir_factor"});

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

  // Now determine virtualities and rebalance
  // Need back-to-back kinematics with a q-qbar pair
  // This seems always to be the case for the epem gun but check nevertheless
  if(p62.size() == 2 && std::abs(p62[0].e() - 0.5*eCM) < 0.001 && std::abs(p62[1].e() - 0.5*eCM) < 0.001 && std::abs(p62[0].id()) < 6 && std::abs(p62[1].id()) < 6){

    // Virtualities of the two partons
    double q1 = 0.;
    double q2 = 0.;
    const double QS = 0.9;

    //Find initial virtuality one parton at a time
    for(int pass=0; pass<2; ++pass){

      double mass = p62[pass].m0();
      // this part is for the case that light quarks are considered massless
      /*if(std::abs(p62[pass].id()) < 4){
        mass = 0.;
      }*/
      double max_vir = (0.25*eCM*eCM - mass*mass) * std::abs(vir_factor);
      double min_vir = (0.5 * QS * QS ) * (1.0 + std::sqrt(1.0 + 4.0 * mass*mass / (QS*QS)));

      double tQ2 = 0.;

      if (max_vir <= QS * QS){
        tQ2 = 0.0;
      }else if(max_vir < min_vir){
        tQ2 = QS * QS;
      }else{

        double numer = 1.0;
        double random = ZeroOneDistribution(*GetMt19937Generator());
        double ratio = 1.0;
        double diff = ratio;
        if(random > rounding_error){
          diff = (ratio - random) / random;
        }

        if(max_vir >= (QS*QS / 2.) * (1.0 + std::sqrt(1.0 + 2.0 * mass * mass / (QS*QS / 2.)))){
          double g = (QS*QS / 2.) * (1.0 + std::sqrt(1.0 + 2.0 * mass * mass / (QS*QS / 2.)));
          numer = exp(-1.0 * (Cf / 2.0 / pi) * sud_val_QG_w_M(mass,(QS*QS / 2.), g, max_vir, 0.5*eCM));
        }

        if (numer > random){tQ2 = min_vir;}
        else{

          double t_hi = max_vir;
          double t_low = min_vir;
          double t_mid = t_low;

          double denom = 1.0;

          do{
              t_mid = 0.5*(t_low + t_hi);

              if (t_mid < (QS*QS / 2.) * (1.0 + std::sqrt(1.0 + 2.0 * mass * mass / (QS*QS / 2.)))){
                denom = 1.0;
              }else{
                double g = (QS*QS / 2.) * (1.0 + std::sqrt(1.0 + 2.0 * mass * mass / (QS*QS / 2.)));
                denom = exp(-1.0 * (Cf / 2.0 / pi) * sud_val_QG_w_M(mass,(QS*QS / 2.), g, t_mid, 0.5*eCM));
              }

              ratio = numer / denom;

              diff = (ratio - random) / random;

              if (diff < 0.0){t_low = t_mid;}
              else{t_hi = t_mid;}

          }while((abs(diff) > s_approx) && (abs(t_hi - t_low) / t_hi > s_error));

          tQ2 = t_mid;

        }
      }

      if(pass==0){q1=sqrt(tQ2);}
      else if(pass==1){q2=sqrt(tQ2);}
    }

    double modm_sq1 = q1*q1 + p62[0].m0()*p62[0].m0();
    double modm_sq2 = q2*q2 + p62[1].m0()*p62[1].m0();

    if(eCM > sqrt(modm_sq1)+sqrt(modm_sq2)){
      // Check viability condition; should always be satisfied in the massless case if virfactor < 1
      double pnew = 0.5*sqrt((eCM*eCM-modm_sq1-modm_sq2)*(eCM*eCM-modm_sq1-modm_sq2)-4.*modm_sq1*modm_sq2)/eCM;

      auto magnitude = [](const Pythia8::Particle &p) {
        return std::sqrt(p.px() * p.px() + p.py() * p.py() + p.pz() * p.pz());
      };

      double scale1 = pnew/magnitude(p62[0]);
      double scale2 = pnew/magnitude(p62[1]);
      double e1new = sqrt(pnew*pnew + modm_sq1);
      double e2new = sqrt(pnew*pnew + modm_sq2);
      p62[0].e(e1new);
      p62[0].px(p62[0].px()*scale1);
      p62[0].py(p62[0].py()*scale1);
      p62[0].pz(p62[0].pz()*scale1);
      p62[1].e(e2new);
      p62[1].px(p62[1].px()*scale2);
      p62[1].py(p62[1].py()*scale2);
      p62[1].pz(p62[1].pz()*scale2);
    }
  }

  //give partons to the framework
  for (int np = 0; np < p62.size(); ++np) {
    Pythia8::Particle &particle = p62.at(np);

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    VERBOSE(7) << "Adding particle with pid = " << particle.id()
               << ", pT = " << particle.pT() << ", eta = " << particle.eta()
               << ", phi = " << particle.phi() << ", e = " << particle.e();

    VERBOSE(7) << " at x=" << xLoc[1] << ", y=" << xLoc[2] << ", z=" << xLoc[3];

    auto ptn = make_shared<Parton>(0, particle.id(), 0, particle.pT(), particle.eta(), particle.phi(), particle.e(), xLoc);
    ptn->set_color(particle.col());
    ptn->set_anti_color(particle.acol());
    ptn->set_max_color(1000 * (np + 1));

    if(p62.size() == 2 && std::abs(particle.id()) < 6){
      double mean_form_time = (2.*ptn->e()) / (ptn->e()*ptn->e()
                            - ptn->px()*ptn->px() - ptn->py()*ptn->py()
                            - ptn->pz()*ptn->pz() - ptn->restmass()*ptn->restmass()
                            + rounding_error) / fmToGeVinv;
      ptn->set_form_time(mean_form_time);
      ptn->set_mean_form_time();

      double velocity[4];
      velocity[0] = 1.0;
      for (int j = 1; j <= 3; j++) {
        velocity[j] = ptn->p(j) / ptn->e();
      }
      ptn->set_jet_v(velocity);
    }
    AddParton(ptn);
  }

  //event.list();

  VERBOSE(8) << GetNHardPartons();
}

double epemGun::sud_val_QG_w_M(double M, double h0, double h1, double h2, double E1) {
  double val, h, intg, hL, hR, diff, intg_L, intg_R, span;

  val = 0.0;

  h = (h1 + h2) / 2.0;

  span = (h2 - h1) / h2;

  val = alpha_s(h) * sud_z_QG_w_M(M, h0, h, E1);

  intg = val * (h2 - h1);

  hL = (h1 + h) / 2.0;

  intg_L = alpha_s(hL) * sud_z_QG_w_M(M, h0, hL, E1) * (h - h1);

  hR = (h + h2) / 2.0;

  intg_R = alpha_s(hR) * sud_z_QG_w_M(M, h0, hR, E1) * (h2 - h);

  diff = std::abs((intg_L + intg_R - intg) / intg);

  if ((diff > approx) || (span > error)) {
    intg = sud_val_QG_w_M(M, h0, h1, h, E1) +
           sud_val_QG_w_M(M, h0, h, h2, E1);
  }

  return intg;
}

double epemGun::sud_z_QG_w_M(double M, double cg, double cg1, double E2){
  // this part is for the case that light quarks are considered massless
  /*if(M < 1.){
    if (cg1 < 2.0 * cg) {
      return 0.0;
    };

    double t2 = std::pow(cg1, 2);
    double t6 = std::log(cg);
    double t10 = std::abs(cg - cg1);
    double t11 = std::log(t10);
    double t17 = -1.0 / t2 * (3.0 * cg1 - 6.0 * cg + 4.0 * t6 * cg1 - 4.0 * t11 * cg1) /
        2.0;
    if (t17 < 0.0) {
      cerr << "ERROR: t17 negative in sud_z_QG_w_M = " << t17 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG_w_M");
    }
    return t17;
  }*/

  if (cg1 < 2.0 * cg + M * M / (1.0 + M * M / cg1)) {
    JSINFO << MAGENTA << " returning with cg, cg1 = " << cg << "   " << cg1
           << "    " << E_minimum << "  " << E2;
    return M * M;
  };

  double t1 = 1.0 / cg1;
  double t2 = t1 * cg;
  double t4 = std::pow(1.0 - t2, 2.0);
  double t7 = std::log(t2);
  double t9 = M * M;
  double t10 = t1 * t9;
  double t13 = 1.0 / (t10 + 1.0) * t10;
  double t15 = std::pow(t2 + t13, 2.0);
  double t18 = std::log(1.0 - t2 - t13);
  double t21 = t1 * (-t4 / 2.0 - 1.0 + 2.0 * t2 - 2.0 * t7 + t15 / 2.0 + t13 +
                     2.0 * t18);

  if (t21 < 0.0) {
    cerr << "ERROR: t21 negative in sud_z_QG_w_M = " << t21 << endl;
    throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG_w_M");
  }

  return t21;
}

double epemGun::alpha_s(double q2) {
  double a, L2, q24, c_nf;

  L2 = std::pow(Lambda_QCD, 2);

  q24 = q2 / 4.0;

  c_nf = nf;

  if (q24 > 4.0) {
    c_nf = 4;
  }

  if (q24 > 64.0) {
    c_nf = 5;
  }

  if (q24 > L2) {
    a = 12.0 * pi / (11.0 * Nc - 2.0 * c_nf) / std::log(q24 / L2);
  } else {
    JSWARN << " alpha too large ";
    a = 0.6;
  }

  return a;
  }