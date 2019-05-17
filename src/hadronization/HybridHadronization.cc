/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2019
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

#include "HybridHadronization.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"
#include "JetScapeConstants.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
//RegisterJetScapeModule<HybridHadronization> HybridHadronization::reg("HybridHadronization");

// Initialize static helper here
Pythia8::Pythia HybridHadronization::pythia ("IntentionallyEmpty",false);

//RNG - Mersenne Twist - 64 bit
std::mt19937_64 eng(std::random_device{}());
//std::mt19937_64 eng(1);
//returns a random number between 0 and 1, based on above engine
double ran() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(eng);}

HybridHadronization::HybridHadronization(){
  SetId("HybridHadronization");
  VERBOSE(8);
}

HybridHadronization::~HybridHadronization(){
  VERBOSE(8);
}

//meson&baryon width functions
double HybridHadronization::SigM2_calc(double R2chg, double qm1, double qm2, double qq1, double qq2){return R2chg*(2./3.)*(qm1+qm2)*(qm1+qm2)/(std::abs(qq1*qm2*qm2 + -1.*qq2*qm1*qm1));}
//CAREFUL, qmi's and qqi's need to be consistent for BOTH of these (mass1(2)(3) is the same for both!)...
double HybridHadronization::SigBR2_calc(double R2chg, double qm1, double qm2, double qm3, double qq1, double qq2, double qq3){
	return R2chg*2.*(2./3.)*(2./3.)*(qm1+qm2+qm3)*(qm1+qm2+qm3)/((qm1+qm2) * std::abs(qq1*(qm2+qm3)*(qm3/qm1)+qq2*(qm1+qm3)*(qm3/qm2)+qq3*(qm1+qm2)));
}
double HybridHadronization::SigBL2_calc(double SigBR2, double qm1, double qm2, double qm3){return SigBR2*((qm1*qm2)/(qm1+qm2))/(qm3*(qm1+qm2)/(qm1+qm2+qm3));}

void HybridHadronization::Init(){

  tinyxml2::XMLElement *hadronization= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("JetHadronization" );

  if ( !hadronization ) {
    JSWARN << "Couldn't find tag Jet Hadronization";
    throw std::runtime_error ("Couldn't find tag Jet Hadronization");
  }
  if (hadronization) {
    string s = hadronization->FirstChildElement( "name" )->GetText();
	//std::string s = GetXMLElementText({"JetHadronization", "name"});
	JSDEBUG << s << " to be initialized ...";
	JSINFO<<"Initialize Hybrid Hadronization ...";

    // Read sqrts to know remnants energies
    double p_read_xml = 10000;
    int flagInt=1;
    hadronization->FirstChildElement("eCMforHadronization")->QueryDoubleText(&p_read_xml);
    p_fake = p_read_xml;
    hadronization->FirstChildElement("take_recoil")->QueryIntText(&flagInt);
    take_recoil=flagInt;

    JSDEBUG<<"Initialize HybridHadronization";
    VERBOSE(8);
	
	maxE_level    = 3;			//maximum energy level considered for the recombination (was 3 in recent fortran code, prev. set to 8)
	gmax          = 1.25;		//maximum allowed mass of the gluon (for q-qbar split), in GeV
	xmq           = 0.33;		//light quark mass, in GeV
	xms           = 0.5;		//strange quark mass, in GeV
	hbarc         = 0.197327;	// GeV*fm - maybe just set this as a constant in common?
	dist2cut      = 25.;		//maximum distance [fm] squared for recombination (involving thermal partons) - in lab frame
	sh_recofactor = 1./3.;		//suppression/enhancement factor for shower-shower recombination
	th_recofactor = 1./3.;		//suppression/enhancement factor for shower-thermal recombination
	attempts_max  = 10;			//maximum number of failed attempts to hadronize a single event before we give up.
	
	//quark masses/charges used to get widths from charged radii
	double Qm_ud = 0.338; double Qm_s = 0.486; double Qm_c = 1.55; double Qm_b = 4.73;
	double chg_u = 2./3.; double chg_d = -1./3.;

	//rms charge radii
	//mesons
	double R2chg_Pi  = 0.42 ;
	double R2chg_Phi = 0.21 ;
	double R2chg_K   = 0.34 ;
	double R2chg_Jpi = 0.04 ;
	double R2chg_Ds  = 0.09 ;
	double R2chg_D   = 0.165;
	double R2chg_Ups = 0.032; //????? setting to a linear extrapolation from Jpsi -> Bc -> Ups
	double R2chg_Bc  = 0.036;
	double R2chg_B   = 0.273;
	//baryons
	double R2chg_Nuc  = 0.69 ;
	double R2chg_Omg  = 0.355;
	double R2chg_Xi   = 0.52 ;
	double R2chg_Sig  = 0.61 ;
	double R2chg_Occc = 0.179; //?????
	double R2chg_Occ  = 0.043; //?
	double R2chg_Xicc = 0.049; //?
	double R2chg_Oc   = 0.1  ; //?????
	double R2chg_Xic  = 0.24 ; //?
	double R2chg_Sigc = 0.27 ; //?
	double R2chg_Obbb = 0.001; //?????
	double R2chg_Obbc = 0.001; //?????
	double R2chg_Obb  = 0.001; //?????
	double R2chg_Xibb = 0.001; //?????
	double R2chg_Obcc = 0.001; //?????
	double R2chg_Obc  = 0.001; //?????
	double R2chg_Xibc = 0.001; //?????
	double R2chg_Ob   = 0.6  ; //?????
	double R2chg_Xib  = 0.63 ; //?????
	double R2chg_Sigb = 0.66 ; //?????

	//meson width calculations (r2) - recalc if r2chg is changed on command line...
	double SigPi2  = SigM2_calc(R2chg_Pi,  Qm_ud, Qm_ud, chg_d, chg_u);
	double SigPhi2 = SigM2_calc(R2chg_Phi, Qm_s,  Qm_s,  chg_d, -chg_d)*(2./3.); //normalizing
	double SigK2   = SigM2_calc(R2chg_K,   Qm_s,  Qm_ud, chg_d, chg_u);
	double SigJpi2 = SigM2_calc(R2chg_Jpi, Qm_c,  Qm_c,  chg_u, -chg_u)*(4./3.); //normalizing
	double SigDs2  = SigM2_calc(R2chg_Ds,  Qm_c,  Qm_s,  chg_u, chg_d);
	double SigD2   = SigM2_calc(R2chg_D,   Qm_c,  Qm_ud, chg_u, chg_d);
	double SigUps2 = SigM2_calc(R2chg_Ups, Qm_b,  Qm_b,  chg_d, -chg_d)*(2./3.); //normalizing
	double SigBc2  = SigM2_calc(R2chg_Bc,  Qm_b,  Qm_c,  chg_d, chg_u);
	double SigB2   = SigM2_calc(R2chg_B,   Qm_b,  Qm_ud, chg_d, chg_u); // (treating B_s as B)

	//baryon width calculations (r2) - recalc if r2chg is changed on command line...
	//light/strange baryons
	double SigNucR2 = SigBR2_calc(R2chg_Nuc, Qm_ud, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	double SigNucL2 = SigBL2_calc(SigNucR2,  Qm_ud, Qm_ud, Qm_ud);
	double SigOmgR2 = SigBR2_calc(R2chg_Omg, Qm_s,  Qm_s,  Qm_s,  chg_d, chg_d, chg_d);
	double SigOmgL2 = SigBL2_calc(SigOmgR2,  Qm_s,  Qm_s,  Qm_s );
	double SigXiR2  = SigBR2_calc(R2chg_Xi,  Qm_s,  Qm_s,  Qm_ud, chg_d, chg_d, chg_d);
	double SigXiL2  = SigBL2_calc(SigXiR2,   Qm_s,  Qm_s,  Qm_ud);
	double SigSigR2 = SigBR2_calc(R2chg_Sig, Qm_s,  Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	double SigSigL2 = SigBL2_calc(SigSigR2,  Qm_s,  Qm_ud, Qm_ud);

	//charm baryons
	double SigOcccR2 = SigBR2_calc(R2chg_Occc, Qm_c, Qm_c,  Qm_c,  chg_d, chg_d, chg_d); // ! maybe need to normalize? (just setting all to -1/3 for now)
	double SigOcccL2 = SigBL2_calc(SigOcccR2,  Qm_c, Qm_c,  Qm_c );
	double SigOccR2  = SigBR2_calc(R2chg_Occ,  Qm_c, Qm_c,  Qm_s,  chg_u, chg_u, chg_d);
	double SigOccL2  = SigBL2_calc(SigOccR2,   Qm_c, Qm_c,  Qm_s );
	double SigXiccR2 = SigBR2_calc(R2chg_Xicc, Qm_c, Qm_c,  Qm_ud, chg_u, chg_u, chg_d);
	double SigXiccL2 = SigBL2_calc(SigXiccR2,  Qm_c, Qm_c,  Qm_ud);
	double SigOcR2   = SigBR2_calc(R2chg_Oc,   Qm_c, Qm_s,  Qm_s,  chg_d, chg_d, chg_d); // ! setting all quark charges to -1/3
	double SigOcL2   = SigBL2_calc(SigOcR2,    Qm_c, Qm_s,  Qm_s );
	double SigXicR2  = SigBR2_calc(R2chg_Xic,  Qm_c, Qm_s,  Qm_ud, chg_u, chg_d, chg_u);
	double SigXicL2  = SigBL2_calc(SigXicR2,   Qm_c, Qm_s,  Qm_ud);
	double SigSigcR2 = SigBR2_calc(R2chg_Sigc, Qm_c, Qm_ud, Qm_ud, chg_u, chg_d, chg_u);
	double SigSigcL2 = SigBL2_calc(SigSigcR2,  Qm_c, Qm_ud, Qm_ud);

	//bottom baryons
	double SigObbbR2 = SigBR2_calc(R2chg_Obbb, Qm_b, Qm_b,  Qm_b,  chg_d, chg_d, chg_d);
	double SigObbbL2 = SigBL2_calc(SigObbbR2,  Qm_b, Qm_b,  Qm_b );
	double SigObbcR2 = SigBR2_calc(R2chg_Obbc, Qm_b, Qm_b,  Qm_c,  chg_d, chg_d, chg_d); // ! setting all quark charges to -1/3
	double SigObbcL2 = SigBL2_calc(SigObbcR2,  Qm_b, Qm_b,  Qm_c );
	double SigObbR2  = SigBR2_calc(R2chg_Obb,  Qm_b, Qm_b,  Qm_s,  chg_d, chg_d, chg_d);
	double SigObbL2  = SigBL2_calc(SigObbR2,   Qm_b, Qm_b,  Qm_s );
	double SigXibbR2 = SigBR2_calc(R2chg_Xibb, Qm_b, Qm_b,  Qm_ud, chg_d, chg_d, chg_d);
	double SigXibbL2 = SigBL2_calc(SigXibbR2,  Qm_b, Qm_b,  Qm_ud);
	double SigObccR2 = SigBR2_calc(R2chg_Obcc, Qm_b, Qm_c,  Qm_c,  chg_d, chg_u, chg_u);
	double SigObccL2 = SigBL2_calc(SigObccR2,  Qm_b, Qm_c,  Qm_c );
	double SigObcR2  = SigBR2_calc(R2chg_Obc,  Qm_b, Qm_c,  Qm_s,  chg_d, chg_d, chg_d); // ! flipping c quark charge (all to -1/3)
	double SigObcL2  = SigBL2_calc(SigObcR2,   Qm_b, Qm_c,  Qm_s );
	double SigXibcR2 = SigBR2_calc(R2chg_Xibc, Qm_b, Qm_c,  Qm_ud, chg_d, chg_u, chg_u);
	double SigXibcL2 = SigBL2_calc(SigXibcR2,  Qm_b, Qm_c,  Qm_ud);
	double SigObR2   = SigBR2_calc(R2chg_Ob,   Qm_b, Qm_s,  Qm_s,  chg_d, chg_d, chg_d);
	double SigObL2   = SigBL2_calc(SigObR2,    Qm_b, Qm_s,  Qm_s );
	double SigXibR2  = SigBR2_calc(R2chg_Xib,  Qm_b, Qm_s,  Qm_ud, chg_d, chg_d, chg_d);
	double SigXibL2  = SigBL2_calc(SigXibR2,   Qm_b, Qm_s,  Qm_ud);
	double SigSigbR2 = SigBR2_calc(R2chg_Sigb, Qm_b, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	double SigSigbL2 = SigBL2_calc(SigSigbR2,  Qm_b, Qm_ud, Qm_ud);

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
	//pythia.readString("PartonLevel:FSR=off"); //is this necessary?
  
    // Don't let pi0 decay
    //pythia.readString("111:mayDecay = off");

    // Don't let any hadron decay
    //pythia.readString("HadronLevel:Decay = off");
	
	//additional settings
	//turning off pythia checks for runtime decrease (can be turned back on if necessary, but it shouldn't make much of a difference)
	pythia.readString("Check:event = off");      // is probably a bad idea, but shouldn't really be necessary... will use a bit of runtime on event checks...
	pythia.readString("Check:history = off");    // might be a good idea to set 'off' as it saves runtime - provided we know that we've set up mother/daughter relations correctly...
	
	//making the pythia event checks a little less stringent (PYTHIA documentation already states that LHC events will occasionally violate default constraint, without concern)
	//pythia.readString("Check:epTolWarn = 1e-4");   // setting E/P conservation violation constraint somewhat weaker, just for ease
	//pythia.readString("Check:epTolErr  = 1e-2");   // setting E/P conservation violation constraint somewhat weaker, just for ease
	//pythia.readString("Check:mTolWarn  = 1e-2");   // setting EP/M conservation violation constraint somewhat weaker, just for ease
	//pythia.readString("Check:mTolErr   = 1e-1");   // setting EP/M conservation violation constraint somewhat weaker, just for ease
	
	//setting a decay threshold for subsequent hadron production
	pythia.readString("ParticleDecays:limitTau0 = on");       //When on, only particles with tau0 < tau0Max are decayed
	//pythia.readString("ParticleDecays:tau0Max = 0.000003"); //The above tau0Max, expressed in mm/c :: default = 10. :: default, mayDecay()=true for tau0 below 1000 mm
															  //set to 1E-17sec (in mm/c) to be smaller than pi0 lifetime
	pythia.readString("ParticleDecays:tau0Max = 10.0");
	
	//allowing for partonic space-time information to be used by PYTHIA
	pythia.readString("PartonVertex:setVertex = on");        //this might allow PYTHIA to keep track of partonic space-time information (default was for 'rope hadronization')
	
	//using QCD based color reconnection (original PYTHIA MPI based CR can't be used at hadron level)
	pythia.readString("ColourReconnection:reconnect = on");          //allowing color reconnections (should have been default on, but doing it here for clarity)
	pythia.readString("ColourReconnection:mode = 1");                //sets the color reconnection scheme to 'new' QCD based scheme (TODO: make sure this is better than (2)gluon move)
	pythia.readString("ColourReconnection:forceHadronLevelCR = on"); //allowing color reconnections for these constructed strings!
	//a few params for the QCD based color reconnection scheme are set below.
	pythia.readString("MultipartonInteractions:pT0Ref = 2.15");      //not sure if this is needed for this setup, but is part of the recommended 'default'
	pythia.readString("ColourReconnection:allowDoubleJunRem = off"); //default on - allows directly connected double junction systems to split into two strings
	pythia.readString("ColourReconnection:junctionCorrection = 1.15");
	pythia.readString("ColourReconnection:timeDilationMode = 3");    //allow reconnection if single pair of dipoles are in causal contact (maybe try 5 as well?)
	pythia.readString("ColourReconnection:timeDilationPar = 0.18");  //parameter used in causal interaction between strings (mode set above)(maybe try 0.073?)

    // And initialize
    pythia.init();
	
	//setting up thermal partons...
	//read in thermal partons, THEN do sibling setup...
	for(int i=0;i<HH_thermal.num();++i){
		HH_thermal[i].sibling(i); HH_thermal[i].string_id(-i);
		HH_thermal[i].is_used(true); HH_thermal[i].sibling( findthermalsibling(i, HH_thermal) ); HH_thermal[i].is_used(false);
	}
  }
}

void HybridHadronization::WriteTask(weak_ptr<JetScapeWriter> w){
  VERBOSE(8);
  auto f = w.lock();
  if ( !f ) return; 
  f->WriteComment("Hadronization Module : "+GetId());
}

//TODO: Junction Strings, Thermal Partons
void HybridHadronization::DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut){
  
  JSINFO<<"Start Hybrid Hadronization using both Recombination and PYTHIA Lund string model.";
  pythia.event.reset(); HH_shower.clear();
  
  for(unsigned int ishower=0; ishower < shower.size(); ++ishower){
	for(unsigned int ipart=0; ipart < shower.at(ishower).size(); ++ipart){
		HHparton sh_parton;
		sh_parton.is_shower(true); sh_parton.id(shower.at(ishower).at(ipart)->pid()); sh_parton.orig(0); //sh_parton.string_id(str)
		sh_parton.px(shower.at(ishower).at(ipart)->px()); sh_parton.py(shower.at(ishower).at(ipart)->py());
		sh_parton.pz(shower.at(ishower).at(ipart)->pz()); sh_parton.e(shower.at(ishower).at(ipart)->e());
		sh_parton.x( shower.at(ishower).at(ipart)->x_in().x() ); sh_parton.y( shower.at(ishower).at(ipart)->x_in().y() );
		sh_parton.z( shower.at(ishower).at(ipart)->x_in().z() ); sh_parton.x_t(shower.at(ishower).at(ipart)->x_in().t() );
		sh_parton.mass( sh_parton.e()*sh_parton.e() - sh_parton.px()*sh_parton.px() - sh_parton.py()*sh_parton.py() - sh_parton.pz()*sh_parton.pz() );
		sh_parton.mass( (sh_parton.mass() >= 0.) ? sqrt(sh_parton.mass()) : sqrt(-sh_parton.mass()) );
		sh_parton.col( shower.at(ishower).at(ipart)->color() ); sh_parton.acol( shower.at(ishower).at(ipart)->anti_color() );
		
		HH_shower.add(sh_parton);
	}
	JSDEBUG<<"Shower#"<<ishower+1 << ". Number of partons to hadronize so far: " << HH_shower.num();
  }
  JSINFO<<"# Partons to hadronize: " << HH_shower.num();
  
  int attempt_num = 0; bool run_successfully = false;
  while((attempt_num < attempts_max) && (!run_successfully)){
	HH_showerptns = HH_shower;
	//clearing hadrons and remnants collections
	HH_hadrons.clear(); HH_remnants.clear(); HH_pyremn.clear();
	
	//since we 'might' need to reset the thermal partons (if present!)... because we alter the thermal parton collection in the string repair routine
	for(int i=0; i<HH_thermal.num(); ++i){
		HH_thermal[i].is_used(false); HH_thermal[i].is_decayedglu(false); HH_thermal[i].is_remnant(false); HH_thermal[i].used_reco(false); HH_thermal[i].used_str(false);
		HH_thermal[i].orig(1); HH_thermal[i].par(-1); HH_thermal[i].status(0); HH_thermal[i].col(0); HH_thermal[i].acol(0);
		HH_thermal[i].PY_par1(-1); HH_thermal[i].PY_par2(-1); HH_thermal[i].PY_dau1(-1); HH_thermal[i].PY_dau2(-1); HH_thermal[i].PY_stat(23);
		HH_thermal[i].string_id(-i); HH_thermal[i].is_strendpt(false); HH_thermal[i].pos_str(1); HH_thermal[i].endpt_id(0);
	}
	
	//checking the shower for any color singlet particles to just dump into PYTHIA
	//also handling colored non-partonic particles
	int i_show=0;
	while(i_show<HH_showerptns.num()){
		if(HH_showerptns[i_show].id() == 21){++i_show; continue;} //is gluon
		else if(std::abs(HH_showerptns[i_show].id()) <= 6){++i_show; continue;} //is quark
		else if((std::abs(HH_showerptns[i_show].id()) >= 1103) && (std::abs(HH_showerptns[i_show].id()) <= 5503) && ((HH_showerptns[i_show].id()/10)%10 == 0)){++i_show; continue;} //is diquark
		else if(pythia.particleData.colType(HH_showerptns[i_show].id()) == 2){//this is a non-gluon color octet...
			HH_showerptns[i_show].PY_origid(HH_showerptns[i_show].id()); HH_showerptns[i_show].id(21); ++i_show; continue;
		}//this is a non-(simple)quark color triplet... (pid 42(scalar leptoquark), 4000001(1-6, excited quarks), and SUSY particles)
		else if(std::abs(pythia.particleData.colType(HH_showerptns[i_show].id())) == 1){
			HH_showerptns[i_show].PY_origid( HH_showerptns[i_show].id() );
			HH_showerptns[i_show].id( 1*(2*std::signbit(-pythia.particleData.colType(HH_showerptns[i_show].id()))-1) ); ++i_show; continue;
		} //and the two above should catch 'colored technihadrons' too!
		else if(HH_showerptns[i_show].id() == 90){HH_showerptns.remove(i_show); continue;} //PYTHIA specific pid for the 'system' as a whole, shouldn't actually be here.  Dumping it.
		//is none of the above (a colorless object of some sort (a hadron, photon, lepton...))
		HHhadron had; had.id( HH_showerptns[i_show].id() ); had.orig( HH_showerptns[i_show].orig() ); had.mass( HH_showerptns[i_show].mass() );
		had.pos(HH_showerptns[i_show].pos()); had.P(HH_showerptns[i_show].P());
		HH_hadrons.add(had); HH_showerptns.remove(i_show);
	}
/*	
	std::cout << "\n\n Start(initial shower):\n";
	for(int i=0;i<HH_showerptns.num();++i){
		std::cout << HH_showerptns[i].id() << ", " << HH_showerptns[i].is_used() << ":" << HH_showerptns[i].status() << ":" << HH_showerptns[i].col() << ":" << HH_showerptns[i].acol() << ":";
		std::cout << HH_showerptns[i].px() << ", " << HH_showerptns[i].py() << ", " << HH_showerptns[i].pz() << ", " << HH_showerptns[i].e() << ", ";
		std::cout << HH_showerptns[i].x() << ", " << HH_showerptns[i].y() << ", " << HH_showerptns[i].z() << ", " << HH_showerptns[i].x_t() << "\n";
	}
	std::cout << "End start(initial shower)\n\n"; std::flush(std::cout);
*/	
	//setting up the strings appropriately for showers - assumes that color tags are set.
	//if there are colored particles without set color tags, it will dump those partons in particular into a single string
	stringform();
	
	//running recombination
	recomb();
	
	//function to prepare strings for input into PYTHIA8
	//will need a final reindexing for py_remn (all PY_par#, PY_dau# will need to be += 1) AFTER this function (either here, or in invoke_py)
	//when recursive 'workaround' for large strings is properly handled/removed, then can put this reindexing inside the function
	if(HH_remnants.num()){bool cut = (attempt_num > -1) ? true : false; stringprep(HH_remnants, HH_pyremn, cut);}
	
	//running remaining partons through PYTHIA8 string fragmentation
	run_successfully = invoke_py();
	
	//for a successful run, go though final hadrons here and set parton parents up
	if(run_successfully){
		//partons that have parh set need to have the reco parents of that hadron
		//partons that do not have parh set need to have the parents vector reset to the original shower/thermal partons
		//until this is done, the par(i) for each parton is *wrong*
		for(int iHad=0; iHad<HH_hadrons.num(); ++iHad){
			//final state string_frag hadron
			if((HH_hadrons[iHad].is_final()) && (HH_hadrons[iHad].parh() == -1) && (HH_hadrons[iHad].parents.size() > 0)){
				std::vector<int> py_remn_parents = HH_hadrons[iHad].parents; HH_hadrons[iHad].parents.clear();
				//just making sure that we got all the partons from the parent string as parents (if there were junctions, this might not have been complete)
				for(int i=0;i<HH_pyremn.num();++i){if((HH_pyremn[i].string_id() == HH_pyremn[py_remn_parents[0]].string_id()) && (HH_pyremn[i].orig() != -1)){py_remn_parents.push_back(i);}}
				std::sort(py_remn_parents.begin(), py_remn_parents.end()); py_remn_parents.erase(std::unique(py_remn_parents.begin(), py_remn_parents.end()), py_remn_parents.end());
				std::vector<int> temp;
				for(int i=0;i<py_remn_parents.size();++i){temp.push_back(HH_pyremn[py_remn_parents[i]].par());}
				std::sort(temp.begin(), temp.end()); temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
				for(int i=0;i<temp.size();++i){if(temp[i]<-1){++temp[i];} else if(temp[i]==-1){temp[i]=-99999;}}
				HH_hadrons[iHad].parents = temp;
				
				//a last check to make sure that any strings with junctions and thermal partons are accounted for
				for(int ipar=0;ipar<py_remn_parents.size();++ipar){if(HH_pyremn[py_remn_parents[ipar]].is_thermal()){HH_hadrons[iHad].is_shth(true);HH_hadrons[iHad].is_shsh(false);break;}}
			}
			//final state reco hadron
			else if((HH_hadrons[iHad].is_final()) && (HH_hadrons[iHad].parh() >= 0) && (HH_hadrons[iHad].parh() < HH_hadrons.num()) && (HH_hadrons[iHad].parh() != iHad)){
				HH_hadrons[iHad].parents = HH_hadrons[HH_hadrons[iHad].parh()].parents;
			}
		}
	}
	
	++attempt_num;
  }
	
  if(!run_successfully){HH_hadrons.clear(); JSWARN << "This event could not be hadronized.";}
  
  for(unsigned int iHad=0; iHad<HH_hadrons.num(); ++iHad){
	if(HH_hadrons[iHad].is_final()){
		int idH = HH_hadrons[iHad].id(); double mH = HH_hadrons[iHad].mass();
		FourVector p(HH_hadrons[iHad].P()); FourVector x(HH_hadrons[iHad].pos());
		hOut.push_back(std::make_shared<Hadron> (Hadron (0,idH,1,p,x,mH)));
	}
  }
  JSINFO<<"#Showers hadronized together: " << shower.size() << ". There are " << hOut.size() << " hadrons and " << pOut.size() << " partons after Hybrid Hadronization";
  
}

//was used in the pygen code to create ordered strings (and assigned string ids) from color tags
//currently has junction functionality commented out until I know how to incorporate it within JETSCAPE
//junction functionality WILL be needed in the event that any given string configurations contain them!
void HybridHadronization::stringform(){
	
	int nstr=1;
	
	for(int i=0; i<HH_showerptns.num(); ++i){
		if(HH_showerptns[i].string_id() != 0){continue;}
		
		std::vector<int> colors;
		if(HH_showerptns[i].col()  > 0){colors.push_back(HH_showerptns[i].col()) ;}
		if(HH_showerptns[i].acol() > 0){colors.push_back(HH_showerptns[i].acol());}
		HH_showerptns[i].string_id(nstr); ++nstr;
		bool newcolor = true;
		while(newcolor){
			newcolor = false;
			//checking all final particles not in a string for new color id matches...
			int j=0; while(j<HH_showerptns.num()){
				if(HH_showerptns[j].string_id() != 0){++j; continue;}
				for(int k=0;k<colors.size();++k){
					if(     colors[k]==HH_showerptns[j].col()  && (HH_showerptns[j].acol()>0)){
						colors.push_back(HH_showerptns[j].acol()); newcolor=true; HH_showerptns[j].string_id(HH_showerptns[i].string_id()); j=-1; break;
					}
					else if(colors[k]==HH_showerptns[j].acol() && (HH_showerptns[j].col() >0)){
						colors.push_back(HH_showerptns[j].col()) ; newcolor=true; HH_showerptns[j].string_id(HH_showerptns[i].string_id()); j=-1; break;
					}
					else if(colors[k]==HH_showerptns[j].col() || colors[k]==HH_showerptns[j].acol()){HH_showerptns[j].string_id(HH_showerptns[i].string_id()); break;}
				}
				++j;
			}
			/* //this code WILL be needed in the event that any input string configurations contain junctions!
			//checking all the junctions in the event to see if there are any new color id matches...
			for(int iJun=0;iJun<pythia.event.sizeJunction();++iJun){
				bool thisjuncolors = false;
				for(int j=0;j<3;++j){for(int k=0;k<colors.size();++k){if(colors[k]==pythia.event.colJunction(iJun, j)){thisjuncolors=true;}}}
				if(!thisjuncolors){continue;}
				for(int j=0;j<3;++j){
					bool addcolor = true;
					for(int k=0;k<colors.size();++k){if(colors[k]==pythia.event.colJunction(iJun,j)){addcolor = false; break;}}
					if(addcolor){colors.push_back(pythia.event.colJunction(iJun,j)); newcolor=true;}
				}
			}
			*/
		}
	}
	
	
	//now, finding all partons in string 'i'(out of nstr-1 total strings) and reordering them
	//for "string" i=0 final partons - these are noncolored particles that should just be written out as-is (shouldn't actually be present...)
	for(int i=1;i<nstr;++i){
		
		std::vector<bool> is_used; std::vector<int> ptns_new;
		std::vector<int> ptns_strnow;
		for(int j=0; j<HH_showerptns.num(); ++j){if(HH_showerptns[j].string_id()==i){ptns_strnow.push_back(j); is_used.push_back(false);}}
		//std::vector<int> usedJuns; for(int j=0; j<NumJunctions; ++j){usedJuns.push_back(false);} int numJunsused = 0;
		
		//find 'a' quark in the event, if there are none, just grab a gluon.
		std::vector<int> stack;
		bool readcolor=true; int firstcol = 0;
		for(int j=0;j<ptns_strnow.size();++j){
			if(!(HH_showerptns[ptns_strnow[j]].id() == 21)){
				if((HH_showerptns[ptns_strnow[j]].id() > 0 && HH_showerptns[ptns_strnow[j]].id() <= 6) || (HH_showerptns[ptns_strnow[j]].id() < -6)){
					stack.push_back(ptns_strnow[j]); is_used[j]=true; firstcol=HH_showerptns[stack[0]].acol(); break;
				}
				else{stack.push_back(ptns_strnow[j]); is_used[j]=true; firstcol=HH_showerptns[stack[0]].col(); readcolor=false; break;}
			}
			//gluon loop catch
			else if(j==ptns_strnow.size()-1){stack.push_back(ptns_strnow[0]); is_used[0]=true; firstcol=HH_showerptns[stack[0]].acol();}
		}
		ptns_new.push_back(stack[0]);
		
		while(stack.size()>0){
			
			//catching the end of either a string/junction leg or we've circled back on a gluon loop.
			if(      readcolor && (HH_showerptns[stack.back()].col()  == firstcol)){stack.pop_back(); continue;}
			else if(!readcolor && (HH_showerptns[stack.back()].acol() == firstcol)){stack.pop_back(); continue;}
			
			//keeping track if we've added (a) new parton(s) to the stack
			bool found=false;
			
			for(int j=0;j<ptns_strnow.size();++j){
				if(      readcolor && !is_used[j] && (HH_showerptns[stack.back()].col()  == HH_showerptns[ptns_strnow[j]].acol())){
					stack[stack.size()-1] = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
					found=true; break;
				}
				else if(!readcolor && !is_used[j] && (HH_showerptns[stack.back()].acol() == HH_showerptns[ptns_strnow[j]].col() )){
					stack[stack.size()-1] = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
					found=true; break;
				}
			}
			
			//this code WILL be needed in the event that any input string configurations contain junctions!
/*			if(!found && ((pythia.event.sizeJunction() - numJunsused) > 0)){
				int iJun(0), iCol(0);
				if(readcolor){
					for(int iJ=0;iJ<pythia.event.sizeJunction();++iJ){
						if(usedJuns[iJ]){continue;}
						for(int iC=0;iC<3;++iC){if(pythia.event.colJunction(iJ,iC)==HH_showerptns[stack.back()].col() ){iJun=iJ;iCol=iC;found=true;break;}}if(found){break;}
					}
				}
				else{
					for(int iJ=0;iJ<pythia.event.sizeJunction();++iJ){
						if(usedJuns[iJ]){continue;}
						for(int iC=0;iC<3;++iC){if(pythia.event.colJunction(iJ,iC)==HH_showerptns[stack.back()].acol()){iJun=iJ;iCol=iC;found=true;break;}}if(found){break;}
					}
				}
				
				//find the next two partons here to add to the stack
				if(found){
					usedJuns[iJun] = true; ++numJunsused;
					readcolor = !readcolor;
					
					int cols[2]={0,0};
					int offset=0; for(int j=0;j<3;++j){if(j==iCol){offset=1;continue;} cols[j-offset]=pythia.event.colJunction(iJun,j);}
					if(cols[0]>cols[1]){int temp=cols[0]; cols[0]=cols[1]; cols[1]=temp;}
					
					if(readcolor){
						for(int j=0;j<ptns_strnow.size();++j){if(!is_used[j] && (cols[0]==HH_showerptns[ptns_strnow[j]].acol())){
							stack[stack.size()-1] = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
							break;
						}}
						for(int j=0;j<ptns_strnow.size();++j){if(!is_used[j] && (cols[1]==HH_showerptns[ptns_strnow[j]].acol())){
							stack.push_back(ptns_strnow[j]); is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
							break;
						}}
					}
					else{
						for(int j=0;j<ptns_strnow.size();++j){if(!is_used[j] && (cols[0]==HH_showerptns[ptns_strnow[j]].col())){
							stack[stack.size()-1] = ptns_strnow[j]; is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
							break;
						}}
						for(int j=0;j<ptns_strnow.size();++j){if(!is_used[j] && (cols[1]==HH_showerptns[ptns_strnow[j]].col())){
							stack.push_back(ptns_strnow[j]); is_used[j]=true; ptns_new.push_back(stack.back()); //--stack.size(); ++stack.size();
							break;
						}}
					}
				}
			}*/
			if(!found){stack.pop_back();}
		}
		
		//ordering this string
		int endpt = 0;
		for(int j=0;j<ptns_new.size();++j){
			HH_showerptns[ptns_new[j]].pos_str(j);
			HH_showerptns[ptns_new[j]].endpt_id( (HH_showerptns[ptns_new[j]].id() == 21) ? 0 : ++endpt );
			HH_showerptns[ptns_new[j]].is_strendpt( (HH_showerptns[ptns_new[j]].id() == 21) ? false : true );
		}
	}
	
	//handling any remaining colored objects that have not had color tags set
	//these are all just thrown into a single string, which may be 'cut' up later on in the string_prep function
	//these will not be 'ordered' within the string nicely...
	int numendpoint = 0; int num_str0 = -1;
	for(int i=0;i<HH_showerptns.num();++i){
		if(HH_showerptns[i].string_id() != 0){continue;}
		int iendpoint = 0;
		bool is_endpoint = false;
		if(HH_showerptns[i].id() != 21){iendpoint = ++numendpoint; is_endpoint = true;}
		HH_showerptns[i].string_id(0);
		HH_showerptns[i].pos_str(++num_str0);
		HH_showerptns[i].endpt_id(iendpoint);
		HH_showerptns[i].is_strendpt(is_endpoint);
	}
	
	//reordering HH_showerptns - based on string, and position in the string
	//isn't strictly necessary to do this here, as the critical pos_str and endpt_id variables have been set, but is 'nice'
	//std::stable_sort(&HH_showerptns[0], (&HH_showerptns[HH_showerptns.num()-1])+1, strid_compare);
	std::stable_sort(&HH_showerptns[0], (&HH_showerptns[HH_showerptns.num()-1])+1, [](const HHparton& parton1, const HHparton& parton2){return (parton1.string_id() < parton2.string_id());});
	//now that the list is sorted based on string id, going to sort the partons in each string based on the position of the partons in the string
	for(int i=0; i<HH_showerptns.num(); ++i){
		int start, prev_pos, cur_pos, lastfix;
		lastfix = 0; if(i == HH_showerptns.num()-1){lastfix = 1;}
		cur_pos = HH_showerptns[i].string_id();
		if(i==0){prev_pos = HH_showerptns[0].string_id(); start = 0;}
		if(cur_pos != prev_pos || i == HH_showerptns.num()-1){
			std::stable_sort(&HH_showerptns[start],&HH_showerptns[i]+lastfix, [](const HHparton& parton1, const HHparton& parton2){return (parton1.pos_str()   < parton2.pos_str()  );});
			start = i; prev_pos = HH_showerptns[i].pos_str();
		}
	}
}

void HybridHadronization::recomb(){
	
	//declaring a few needed values (based on values in init)
	double hbarc2 = hbarc*hbarc;
	
	//should create a new parton collection here, one with only shower quarks (inc. from gluon splitting)
	//to consider for recombination - then afterwards can reform gluons and output remnants
	parton_collection showerquarks;
	
	//clearing remnants, in case it hasn't been done before
	HH_remnants.clear();
	
	//constructing a list of all the strings in the event
	std::vector<int> list_strs;
	//adding the first string to the list
	list_strs.push_back(HH_showerptns[0].string_id());
	//looping over all the partons in the event, and writing each 'unique' string to the list
	for(int i=0;i<HH_showerptns.num();++i){
		bool str_match = false;
		for(int j=0;j<list_strs.size();++j){
			if(HH_showerptns[i].string_id() == list_strs[j]){str_match = true;}
		}
		if(!str_match){list_strs.push_back(HH_showerptns[i].string_id());}
	}
	
//*********************************************************************************************************************
//		Splitting gluons into q-qbar pairs for recombination
//*********************************************************************************************************************
	
	//starting with shower partons
	for(int i_pt=0; i_pt<HH_showerptns.num(); ++i_pt){
		//if parton is a quark, stick it into showerquarks - and set it's parent id to the quark in the original shower
		if((std::abs(HH_showerptns[i_pt].id()) <= 5 ) && (HH_showerptns[i_pt].PY_origid() == 0)){showerquarks.add(HH_showerptns.partons[i_pt]); showerquarks[showerquarks.num()-1].par(i_pt);}
		//else parton is a gluon, decay into q-qbar and stick those into the collection; set parent id appropriately for both
		else if((HH_showerptns[i_pt].id() == 21) && (HH_showerptns[i_pt].PY_origid() == 0)){
			//setting is_decayedglu; not going to set used until one of it's quarks has been used
			//will set the status to -99 here, this will be changed to -1 after first quark, and to 1 after second
			HH_showerptns[i_pt].is_decayedglu(true); HH_showerptns[i_pt].status(-99);
			
			//choosing a gluon mass - if implemented in the future, can (should) read this from the gluon entry itself (or set if necessary)
			//maybe discard gluon if it is under some threshold of mass (eg < pion?)
			//temporarily saving previously set mass here - here's a good place to check if this is even necessary?
			double temp_glumass = HH_showerptns[i_pt].mass();
			HH_showerptns[i_pt].mass( 2.*xmq + (gmax - 2.*xmq)*ran() );	// gluon virtuality
			
			//gluon decay function reads in the gluon (and the overwritten random mass), and writes the output q-qbar pair to qpair
			parton_collection qpair;
			gluon_decay(HH_showerptns[i_pt], qpair);
			
			//swapping back original gluon mass
			HH_showerptns[i_pt].mass(temp_glumass);
			//std::swap(HH_showerptns[i_pt].mass, temp_glumass);
			
			//setting the parents of the q-qbar pair to the original gluon (and other vars appropriately here)
			qpair[0].par(i_pt); qpair[1].par(i_pt);
			qpair[0].is_shower(true); qpair[1].is_shower(true);
			qpair[0].orig( HH_showerptns[i_pt].orig() ); qpair[1].orig( HH_showerptns[i_pt].orig() );
			qpair[0].string_id ( HH_showerptns[i_pt].string_id()); qpair[1].string_id ( HH_showerptns[i_pt].string_id());
			qpair[0].pos_str   ( HH_showerptns[i_pt].pos_str());   qpair[1].pos_str   ( HH_showerptns[i_pt].pos_str());
			
			//adding these partons to the collection
			showerquarks.add(qpair);
			
			//setting up the sibling relations in showerquarks
			showerquarks[showerquarks.num()-1].sibling( showerquarks.num()-2 );
			showerquarks[showerquarks.num()-2].sibling( showerquarks.num()-1 );
		}
		//else it's not a quark(u,d,s,c,b) or gluon and we're skipping it.
		else{
			//JSINFO << "\n\nThere is a parton that is not a quark(u,d,s,c,b) or gluon in the input shower.  Skipping parton in recombination.\n\n";
		}
	}

//*********************************************************************************************************************
//		Finished splitting gluons
//*********************************************************************************************************************

//*********************************************************************************************************************
//		Beginning of Recombination routine
//*********************************************************************************************************************
	
	//previously, the recombination routine ran sequentially down the list of quarks - potential bias there...
	//instead, we will randomly permute an integer array with 'n' entries from 1 to n (actually from 0 to n-1)
	//the array will be used to access the i'th element of the showerquarks collection
	//done as showerquarks[perm0_sharray[i]]
	//since we have separate shower and thermal arrays, we will force the first quark to always be from the shower
	//
	//construct perm1 array of size showerquarks.num() and perm2 array of size showerquarks.num() + HH_thermal.num()
	//to access an element from either shower or thermal partons
	//there will be showerquarks.num() positive entries and HH_thermal.num() negative entries (going to exclude 0 in these, just -1 i's)
	//thus, perm2 = { 3, -42, -33,...} will access i=2 from showerquarks, then 41 from thermal, then 32 from thermal...
	//
	//lastly, need to bypass scenarios where the quark has been used, or we're trying to use the same quark twice
	//to determine if we're trying to use the same quark twice, check if perm[i]=perm[j]
	//if used=true OR p[i]=p[j]?, then skip this attempt (continue works well for these cases?)
	
	//constructing permutation arrays; shower quarks are > 0, thermal are < 0
	//perm2 will always access element [std::abs(perm2[i]) - 1]
	int perm1[showerquarks.num()], perm2[showerquarks.num()+HH_thermal.num()];
	for(int i=0;i<showerquarks.num();++i){perm1[i]=i;}
	for(int i=0;i<HH_thermal.num();++i){perm2[i]=(i-HH_thermal.num());}
	for(int i=HH_thermal.num();i<showerquarks.num()+HH_thermal.num();++i){perm2[i]=(i-HH_thermal.num()+1);}
	
	//permuting these arrays using Fisher-Yates algorithm
	//-- To shuffle an array a of n elements (indices 0..n-1):
	//for i from 0 to n−2 do
	//  j ← random integer such that i ≤ j < n
	//  exchange a[i] and a[j]
	for(int i=0;i<showerquarks.num()-1;++i){
		int ranelement = i+floor((showerquarks.num()-i)*ran());
		int temp = perm1[i]; perm1[i] = perm1[ranelement]; perm1[ranelement] = temp;
	}
	for(int i=0;i<showerquarks.num()+HH_thermal.num()-1;++i){
		int ranelement = i+floor((showerquarks.num()+HH_thermal.num()-i)*ran());
		int temp = perm2[i]; perm2[i] = perm2[ranelement]; perm2[ranelement] = temp;
	}
	
	//looping over all the quarks that we have...
	//'q1' loops over all quarks in the shower
	//'q2' loops over all quarks in the event
	//'q3' loops over all quarks in the event, starting from 'q2' and ending at the last quark
	//when 'q2' is at the last quark, we will not consider quark 'q3' - can only make a meson at that point...
	
	parton_collection considering; int element[3];
	
	for(int q1=0;q1<showerquarks.num();++q1){
		//accessing first considered quark
		//set q1 variables here
		
		//continue; //uncomment to turn off recombination in a computationally friendly way...
		
		//resetting madehadron flag, now that we're going to be considering a new hadron
		bool madehadron = false;
		
		//taking the id from the permutation array, and turning it into quark array element
		element[0] = perm1[q1];
		
		//skipping if current quark has been used or isn't a u,d,s quark or antiquark
		if(showerquarks[element[0]].status() != 0 || showerquarks[element[0]].is_used()){continue;}
		else if(std::abs(showerquarks[element[0]].id()) > 5){JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for recombination, THIS SHOULD NOT HAPPEN!"; continue;}
		
		//assigning quark values to considering - and setting the status to -991
		//if at the end of the check, we make a hadron, we will set all -99* entries to 1, else we'll set back to 0
//		if(considering.num() > 0){JSWARN << "SOMETHING HAS GONE VERY WRONG WITH COLLECTION; HAS " << considering.num() << " ENTRIES (SHOULD BE 0)"; considering.clear();}
		
		considering.add(showerquarks[element[0]]);
		showerquarks[element[0]].status(-991);
		
		for(int q2=0;q2<showerquarks.num()+HH_thermal.num();++q2){
			//set q2 variables here - if we can form a meson, then skip q3 loop
			//also skip q3 loop if q2 is at last quark
			
			double recofactor2 = 0.;
			
			//accessing the second considered quark
			//this will skip over non-quark entries in HH_thermal
			if(perm2[q2]>0){/*is shower quark*/ element[1] = perm2[q2] - 1;}
			else{/*is thermal quark*/ element[1] = perm2[q2] + 1; if(std::abs(HH_thermal[-element[1]].id()) > 5){continue;}}
			
			//checking to see if quark is in shower or thermal
			//only need to check if quark is the same as the previous quark IFF it is in the shower
			//only need to check if quark is from the same gluon if it is from a gluon decay in the shower...
			//want to check if we've messed up and are accessing a gluon if it is in original shower/thermal
			//need to check if used for all cases...
			if(perm2[q2]>0){
				//skipping if current quark has been used
				if(showerquarks[element[1]].status() != 0 || showerquarks[element[1]].is_used()){continue;}
				//skipping if current quark is the same as q1
				else if(element[0]==element[1]){continue;}
				//skipping if the current quark is from the same gluon as q1
				else if((showerquarks[element[0]].par() != -1) && (showerquarks[element[0]].par() == showerquarks[element[1]].par())){continue;}
				//skipping if the current quark is not in the same string as q1 (both are shower partons)
				else if(showerquarks[element[1]].string_id() != showerquarks[element[0]].string_id()){continue;}
				//skipping if current quark is not a u,d,s quark
				else if(std::abs(showerquarks[element[1]].id()) > 5){JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for recombination, THIS SHOULD NOT HAPPEN!"; continue;}
				
				considering.add(showerquarks[element[1]]);
				showerquarks[element[1]].status(-992);
				recofactor2 = sh_recofactor*sh_recofactor;
			}
			else if(perm2[q2]<0){
				//skipping if current quark has been used
				if(HH_thermal[-element[1]].status() != 0 || HH_thermal[-element[1]].is_used()){continue;}
				//skipping if current quark is not a u,d,s quark
				else if(std::abs(HH_thermal[-element[1]].id()) > 5){JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for recombination, THIS SHOULD NOT HAPPEN!"; continue;}
				
				//checking distance cut ONLY if this is a thermal parton - skip if dist2 > dist2cut
				FourVector pos_ptn1 = considering[0].pos(); FourVector pos_ptn2 = HH_thermal[-element[1]].pos();
				double dt = pos_ptn1.t() - pos_ptn2.t();
				if(dt > 0.){
					double dt_E = dt/HH_thermal[-element[1]].e(); // P/E * dT = dist = P*(dT/E)
					pos_ptn2.Set(pos_ptn2.x()+HH_thermal[-element[1]].px()*dt_E,pos_ptn2.y()+HH_thermal[-element[1]].py()*dt_E,pos_ptn2.z()+HH_thermal[-element[1]].pz()*dt_E,0.);
				}
				else{
					double dt_E = -dt/considering[0].e();
					pos_ptn1.Set(pos_ptn1.x()+considering[0].px()*dt_E,pos_ptn1.y()+considering[0].py()*dt_E,pos_ptn1.z()+considering[0].pz()*dt_E, 0.);
				}
				if(dif2(pos_ptn1,pos_ptn2) > dist2cut){continue;}
				
				considering.add(HH_thermal[-element[1]]);
				HH_thermal[-element[1]].status(-992);
				recofactor2 = th_recofactor*sh_recofactor;
			}
			else{JSWARN << "SOMETHING WENT HORRIBLY WRONG - DO NOT KNOW WHERE CURRENT QUARK CAME FROM?!";}
			
			//now that we have two 'acceptable' quarks to check hadron formation against,
			//there is no reason to bother checking if we can make a baryon if we have a q-qbar at this point...
			//will skip third loop in this case - otherwise we will check if we can make a baryon...
			if((considering[0].id()*considering[1].id() > 0) && (q2 < showerquarks.num()+HH_thermal.num()-1)){for(int q3=q2+1;q3<showerquarks.num()+HH_thermal.num();++q3){
				
				double recofactor3 = recofactor2;
				
				//removing all but the first two entries in the considering collection...
				//this should have been done before, but have this here just in case - remove if this doesn't ever trigger
				//while(considering.num() > 2){considering--;}
				while(considering.num() > 2){considering.partons.pop_back();}
				
				//accessing the third considered quark
				//this will skip over non-quark entries in HH_thermal
				if(perm2[q3]>0){/*is shower quark*/ element[2] = perm2[q3] - 1;}
				else{/*is thermal quark*/ element[2] = perm2[q3] + 1; if(std::abs(HH_thermal[-element[2]].id()) > 5){continue;}}
				
				//now that we have q3, we need to check if it is valid:
				//q3 needs to be checked if used (all cases)
				//q3 needs to be checked if it is a erroneously accessed parton (not u,d,s quark)
				//q3 needs to be checked against q1 to see if it is the same, or from the same gluon
				//q3 does not need to be checked against q2 to see if it is the same
				//q3 does need to be checked to see if it is from the same gluon as q2
				//q3 needs to be checked to make sure that it can form a baryon with q1 and q2 (checking against either is ok)
				//check:  used, sameq1, sameg_q1, sameg_q2, isglu
				if(perm2[q3]>0){
					//skipping if the current quark cannot make a baryon with other two quarks
					if(showerquarks[element[2]].id()*considering[0].id() < 0){continue;}
					//skipping if current quark has been used
					if(showerquarks[element[2]].status() != 0 || showerquarks[element[2]].is_used()){continue;}
					//skipping if current quark is the same as q1
					else if(element[0]==element[2]){continue;}
					//skipping if the current quark is from the same gluon as q1
					else if((showerquarks[element[0]].par() != -1) && (showerquarks[element[0]].par() == showerquarks[element[2]].par())){continue;}
					//skipping if the current quark is from the same gluon as q2
					else if((perm2[q2] > 0) && (showerquarks[element[1]].par() != -1) && (showerquarks[element[1]].par() == showerquarks[element[2]].par())){continue;}
					//skipping if the current quark is not in the same string as q1
					//(q2 MUST be in the same string as q1 if it's in the shower, and doesn't need to be checked if thermal)
					else if(showerquarks[element[2]].string_id() != showerquarks[element[0]].string_id()){continue;}
					//skipping if current quark is not a u,d,s quark
					else if(std::abs(showerquarks[element[2]].id()) > 5){JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for recombination, THIS SHOULD NOT HAPPEN!"; continue;}
					
					considering.add(showerquarks[element[2]]);
					showerquarks[element[2]].status(-993);
					recofactor3 = sh_recofactor*recofactor2;
				}
				else if(perm2[q3]<0){
					//skipping if the current quark cannot make a baryon with other two quarks
					if(HH_thermal[-element[2]].id()*considering[0].id() < 0){continue;}
					//skipping if current quark has been used
					if(HH_thermal[-element[2]].status() != 0 || HH_thermal[-element[2]].is_used()){continue;}
					//skipping if the current quark is from the same gluon as q2
					else if((perm2[q2] < 0) && (HH_thermal[-element[1]].par() != -1) && (HH_thermal[-element[1]].par() == HH_thermal[-element[2]].par())){continue;}
					//skipping if current quark is not a u,d,s quark
					else if(std::abs(HH_thermal[-element[2]].id()) > 5){JSWARN << "SOMETHING OTHER THAN u,d,s,c,b WAS considered for recombination, THIS SHOULD NOT HAPPEN!"; continue;}
					
					//checking distance cut ONLY if this is a thermal parton - skip if dist2 > dist2cut
					FourVector pos_ptn1 = considering[0].pos(); FourVector pos_ptn2 = considering[1].pos(); FourVector pos_ptn3 = HH_thermal[-element[2]].pos();
					if(     (pos_ptn1.t() > pos_ptn2.t()) && (pos_ptn1.t() > pos_ptn3.t())){
						double dt_E2 = (pos_ptn1.t() - pos_ptn2.t())/considering[1].e(); double dt_E3 = (pos_ptn1.t() - pos_ptn3.t())/HH_thermal[-element[2]].e();
						pos_ptn2.Set(pos_ptn2.x()+considering[1].px()*dt_E2,pos_ptn2.y()+considering[1].py()*dt_E2,pos_ptn2.z()+considering[1].pz()*dt_E2,0.);
						pos_ptn3.Set(pos_ptn3.x()+HH_thermal[-element[2]].px()*dt_E3,pos_ptn3.y()+HH_thermal[-element[2]].py()*dt_E3,pos_ptn3.z()+HH_thermal[-element[2]].pz()*dt_E3,0.);
					}
					else if((pos_ptn2.t() > pos_ptn1.t()) && (pos_ptn2.t() > pos_ptn3.t())){
						double dt_E1 = (pos_ptn2.t() - pos_ptn1.t())/considering[0].e(); double dt_E3 = (pos_ptn2.t() - pos_ptn3.t())/HH_thermal[-element[2]].e();
						pos_ptn1.Set(pos_ptn1.x()+considering[0].px()*dt_E1,pos_ptn1.y()+considering[0].py()*dt_E1,pos_ptn1.z()+considering[0].pz()*dt_E1,0.);
						pos_ptn3.Set(pos_ptn3.x()+HH_thermal[-element[2]].px()*dt_E3,pos_ptn3.y()+HH_thermal[-element[2]].py()*dt_E3,pos_ptn3.z()+HH_thermal[-element[2]].pz()*dt_E3,0.);
					}
					else{
						double dt_E1 = (pos_ptn3.t() - pos_ptn1.t())/considering[0].e(); double dt_E2 = (pos_ptn3.t() - pos_ptn2.t())/considering[1].e();
						pos_ptn1.Set(pos_ptn1.x()+considering[0].px()*dt_E1,pos_ptn1.y()+considering[0].py()*dt_E1,pos_ptn1.z()+considering[0].pz()*dt_E1,0.);
						pos_ptn2.Set(pos_ptn2.x()+considering[1].px()*dt_E2,pos_ptn2.y()+considering[1].py()*dt_E2,pos_ptn2.z()+considering[1].pz()*dt_E2,0.);
					}
					if((dif2(pos_ptn3,pos_ptn1) > dist2cut) || (dif2(pos_ptn3,pos_ptn2) > dist2cut) || (dif2(pos_ptn1,pos_ptn2) > dist2cut)){continue;}
					
					considering.add(HH_thermal[-element[2]]);
					HH_thermal[-element[2]].status(-993);
					recofactor3 = th_recofactor*recofactor2;
				}
				else{JSWARN << "SOMETHING WENT HORRIBLY WRONG - DO NOT KNOW WHERE CURRENT QUARK CAME FROM?!";}
				
				//now that we *could* form a baryon, now we check if we actually do form one
				//baryon momentum
				FourVector Pbaryon;
				Pbaryon.Set(considering[0].px()+considering[1].px()+considering[2].px(),considering[0].py()+considering[1].py()+considering[2].py(),considering[0].pz()+considering[1].pz()+considering[2].pz(),0.);
				
				//baryon(CM) velocity
				FourVector betaB; //really p[i]/e below
				betaB.Set(Pbaryon.x()/(considering[0].e()+considering[1].e()+considering[2].e()),Pbaryon.y()/(considering[0].e()+considering[1].e()+considering[2].e()),Pbaryon.z()/(considering[0].e()+considering[1].e()+considering[2].e()),0.);
				betaB.Set(betaB.x(),betaB.y(),betaB.z(),1./(sqrt(1. - (betaB.x()*betaB.x() + betaB.y()*betaB.y() + betaB.z()*betaB.z()))));
				
				//boosting into CM frame
				FourVector pos_BCM[3], p_BCM[3];
				pos_BCM[0] = considering[0].boost_pos(betaB); pos_BCM[1] = considering[1].boost_pos(betaB); pos_BCM[2] = considering[2].boost_pos(betaB);
				  p_BCM[0] = considering[0].boost_P(betaB);     p_BCM[1] = considering[1].boost_P(betaB);     p_BCM[2] = considering[2].boost_P(betaB);
				
				//velocities in CM frame
				FourVector v_BCM[3];
				v_BCM[0].Set(p_BCM[0].x()/p_BCM[0].t(),p_BCM[0].y()/p_BCM[0].t(),p_BCM[0].z()/p_BCM[0].t(),0.); //these are really p[i]/e
				v_BCM[1].Set(p_BCM[1].x()/p_BCM[1].t(),p_BCM[1].y()/p_BCM[1].t(),p_BCM[1].z()/p_BCM[1].t(),0.);
				v_BCM[2].Set(p_BCM[2].x()/p_BCM[2].t(),p_BCM[2].y()/p_BCM[2].t(),p_BCM[2].z()/p_BCM[2].t(),0.);
				
				//propagating quarks until time of youngest quark
				double curtime = std::max(std::max(pos_BCM[0].t(), pos_BCM[1].t()), pos_BCM[2].t());
				FourVector cur_pos[3];
				cur_pos[0].Set(pos_BCM[0].x()+v_BCM[0].x()*(curtime-pos_BCM[0].t()),pos_BCM[0].y()+v_BCM[0].y()*(curtime-pos_BCM[0].t()),pos_BCM[0].z()+v_BCM[0].z()*(curtime-pos_BCM[0].t()),curtime);
				cur_pos[1].Set(pos_BCM[1].x()+v_BCM[1].x()*(curtime-pos_BCM[1].t()),pos_BCM[1].y()+v_BCM[1].y()*(curtime-pos_BCM[1].t()),pos_BCM[1].z()+v_BCM[1].z()*(curtime-pos_BCM[1].t()),curtime);
				cur_pos[2].Set(pos_BCM[2].x()+v_BCM[2].x()*(curtime-pos_BCM[2].t()),pos_BCM[2].y()+v_BCM[2].y()*(curtime-pos_BCM[2].t()),pos_BCM[2].z()+v_BCM[2].z()*(curtime-pos_BCM[2].t()),curtime);
				
				//finding position of CM at curtime
				FourVector pos_CM;
				pos_CM.Set(
				(cur_pos[0].x()*considering[0].mass()+cur_pos[1].x()*considering[1].mass()+cur_pos[2].x()*considering[2].mass())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				(cur_pos[0].y()*considering[0].mass()+cur_pos[1].y()*considering[1].mass()+cur_pos[2].y()*considering[2].mass())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				(cur_pos[0].z()*considering[0].mass()+cur_pos[1].z()*considering[1].mass()+cur_pos[2].z()*considering[2].mass())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				curtime
				);
				
				//finding position of baryon in lab frame
				betaB.Set(-betaB.x(),-betaB.y(),-betaB.z(),betaB.t());
				FourVector pos_lab = HHboost(betaB, pos_CM);
				
				
				//finding relative momenta of partons in CM frame
				FourVector k_rel[2];
				k_rel[0].Set(
				(considering[1].mass()*p_BCM[0].x()-considering[0].mass()*p_BCM[1].x())/(considering[0].mass()+considering[1].mass()),
				(considering[1].mass()*p_BCM[0].y()-considering[0].mass()*p_BCM[1].y())/(considering[0].mass()+considering[1].mass()),
				(considering[1].mass()*p_BCM[0].z()-considering[0].mass()*p_BCM[1].z())/(considering[0].mass()+considering[1].mass()),
				0.);
				k_rel[1].Set(
				(considering[2].mass()*(p_BCM[0].x()+p_BCM[1].x())-(considering[0].mass()+considering[1].mass())*p_BCM[2].x())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				(considering[2].mass()*(p_BCM[0].y()+p_BCM[1].y())-(considering[0].mass()+considering[1].mass())*p_BCM[2].y())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				(considering[2].mass()*(p_BCM[0].z()+p_BCM[1].z())-(considering[0].mass()+considering[1].mass())*p_BCM[2].z())/(considering[0].mass()+considering[1].mass()+considering[2].mass()),
				0.);
				
				//finding relative positions of partons in CM frame
				FourVector pos_rel[2];
				pos_rel[0].Set((cur_pos[0].x()-cur_pos[1].x())/sqrt(2.),(cur_pos[0].y()-cur_pos[1].y())/sqrt(2.),(cur_pos[0].z()-cur_pos[1].z())/sqrt(2.),0.);
				pos_rel[1].Set(
				((cur_pos[0].x()*considering[0].mass()+cur_pos[1].x()*considering[1].mass())/(considering[0].mass()+considering[1].mass())-cur_pos[2].x())*sqrt(2./3.),
				((cur_pos[0].y()*considering[0].mass()+cur_pos[1].y()*considering[1].mass())/(considering[0].mass()+considering[1].mass())-cur_pos[2].y())*sqrt(2./3.),
				((cur_pos[0].z()*considering[0].mass()+cur_pos[1].z()*considering[1].mass())/(considering[0].mass()+considering[1].mass())-cur_pos[2].z())*sqrt(2./3.),
				0.);
				
				
				double SigRB2 = SigNucR2; double SigLB2 = SigNucL2;
				int sortid[3] = {std::abs(considering[0].id()),std::abs(considering[1].id()),std::abs(considering[2].id())};
				std::stable_sort(std::begin(sortid), std::end(sortid), std::greater<int>());
				
				//for particles we don't want to form, setting recofactor3 to 0
				if(     sortid[0] == 3){
					if(     sortid[1] == 3){
						if(     sortid[2] == 3){SigRB2 = SigOmgR2;  SigLB2 = SigOmgL2;}
						else{                   SigRB2 = SigXiR2;   SigLB2 = SigXiL2;}
					}
					else{                       SigRB2 = SigSigR2;  SigLB2 = SigSigL2;}
				}
				else if(sortid[0] == 4){
					if(     sortid[1] == 4){
						if(     sortid[2] == 4){SigRB2 = SigOcccR2; SigLB2 = SigOcccL2;}
						else if(sortid[2] == 3){SigRB2 = SigOccR2;  SigLB2 = SigOccL2;}
						else{                   SigRB2 = SigXiccR2; SigLB2 = SigXiccL2;}
					}
					else if(sortid[1] == 3){
						if(     sortid[2] == 3){SigRB2 = SigOcR2;   SigLB2 = SigOcL2;}
						else{                   SigRB2 = SigXicR2;  SigLB2 = SigXicL2;}
					}
					else{                       SigRB2 = SigSigcR2; SigLB2 = SigSigcL2;}
				}
				else if(sortid[0] == 5){
					if(     sortid[1] == 5){
						if(     sortid[2] == 5){SigRB2 = SigObbbR2; SigLB2 = SigObbbL2; recofactor3=0.;}
						else if(sortid[2] == 4){SigRB2 = SigObbcR2; SigLB2 = SigObbcL2; recofactor3=0.;}
						else if(sortid[2] == 3){SigRB2 = SigObbR2;  SigLB2 = SigObbL2;  recofactor3=0.;}
						else{                   SigRB2 = SigXibbR2; SigLB2 = SigXibbL2; recofactor3=0.;}
					}
					else if(sortid[1] == 4){
						if(     sortid[2] == 4){SigRB2 = SigObccR2; SigLB2 = SigObccL2; recofactor3=0.;}
						else if(sortid[2] == 3){SigRB2 = SigObcR2;  SigLB2 = SigObcL2;  recofactor3=0.;}
						else{                   SigRB2 = SigXibcR2; SigLB2 = SigXibcL2; recofactor3=0.;}
					}
					else if(sortid[1] == 3){
						if(     sortid[2] == 3){SigRB2 = SigObR2;   SigLB2 = SigObL2;}
						else{                   SigRB2 = SigXibR2;  SigLB2 = SigXibL2;}
					}
					else{                       SigRB2 = SigSigbR2; SigLB2 = SigSigbL2;}
				}
				
				//precalc's for Wigner Wavefunction
				//0:x, 1:y, 2:z ::: urho:(rel. between partons 1,2), ulamb:(rel. between partons (1,2),3)
				double urho[3], ulamb[3];
				urho[0] =  0.5*(pos_rel[0].x()*pos_rel[0].x()/SigRB2 + k_rel[0].x()*k_rel[0].x()*SigRB2/hbarc2);
				urho[1] =  0.5*(pos_rel[0].y()*pos_rel[0].y()/SigRB2 + k_rel[0].y()*k_rel[0].y()*SigRB2/hbarc2);
				urho[2] =  0.5*(pos_rel[0].z()*pos_rel[0].z()/SigRB2 + k_rel[0].z()*k_rel[0].z()*SigRB2/hbarc2);
				ulamb[0] = 0.5*(pos_rel[1].x()*pos_rel[1].x()/SigLB2 + k_rel[1].x()*k_rel[1].x()*SigLB2/hbarc2);
				ulamb[1] = 0.5*(pos_rel[1].y()*pos_rel[1].y()/SigLB2 + k_rel[1].y()*k_rel[1].y()*SigLB2/hbarc2);
				ulamb[2] = 0.5*(pos_rel[1].z()*pos_rel[1].z()/SigLB2 + k_rel[1].z()*k_rel[1].z()*SigLB2/hbarc2);
				
				//1D GS Wig. wavefunction
				double wig0[2][3];
				wig0[0][0] = std::exp(-urho[0]);  wig0[0][1] = std::exp(-urho[1]);  wig0[0][2] = std::exp(-urho[2]);
				wig0[1][0] = std::exp(-ulamb[0]); wig0[1][1] = std::exp(-ulamb[1]); wig0[1][2] = std::exp(-ulamb[2]);
				//3D GS Wig. wavefunction
				double WigB[2];
				WigB[0] = wig0[0][0]*wig0[0][1]*wig0[0][2]*wig0[1][0]*wig0[1][1]*wig0[1][2];
				
				//summing up 3D Wig. wavefunctions over nlev excited states
				WigB[1] = 0.; //sumWigB;
				
				double wigE[2][3];
				for(int i=0;i<2;++i){for(int j=0;j<2;++j){wigE[i][j]=0.;}}
				
				wigE[0][0] = wig0[0][0];
				for(int iRx=0; iRx<=maxE_level; ++iRx){
					wigE[0][1] = wig0[0][1];
					for(int iRy=0; iRy<=maxE_level-iRx; ++iRy){
						wigE[0][2] = wig0[0][2];
						for(int iRz=0; iRz<=maxE_level-iRx-iRy; ++iRz){
							wigE[1][0] = wig0[1][0];
							for(int iLx=0; iLx<=maxE_level-iRx-iRy-iRz; ++iLx){
								wigE[1][1] = wig0[1][1];
								for(int iLy=0; iLy<=maxE_level-iRx-iRy-iRz-iLx; ++iLy){
									wigE[1][2] = wig0[1][2];
									for(int iLz=0; iLz<=maxE_level-iRx-iRy-iRz-iLx-iLy; ++iLz){
										WigB[1] += wigE[0][0]*wigE[0][1]*wigE[0][2]*wigE[1][0]*wigE[1][1]*wigE[1][2];
										wigE[1][2] *= ulamb[2]/((double(iLz))+1.);
									}
									wigE[1][1] *= ulamb[1]/((double(iLy))+1.);
								}
								wigE[1][0] *= ulamb[0]/((double(iLx))+1.);
							}
							wigE[0][2] *= urho[2]/((double(iRz))+1.);
						}
						wigE[0][1] *= urho[1]/((double(iRy))+1.);
					}
					wigE[0][0] *= urho[0]/((double(iRx))+1.);
				}
				
				//checking if we form a baryon - there were variables suggesting differing recomb. probabilities could
				//have been included if baryon was strange or light - but they were identical in the fortran code
				//this will have to be altered if we want to include that here.
				
				//Checking if baryon is formed (either ground or excited state)
				double rndbaryon = ran();
				
				if(WigB[1]*recofactor3 >= rndbaryon){
					//*******************************************string repair functionality below*******************************************
					//determining which partons, of the 3 being considered, are endpoints
					int numendpoints = 0; bool is_endpoint[3];  is_endpoint[0]=false;  is_endpoint[1]=false;  is_endpoint[2]=false;
									if(showerquarks[element[0]].is_strendpt()){++numendpoints; is_endpoint[0] = true;}
					if(perm2[q2]>0){if(showerquarks[element[1]].is_strendpt()){++numendpoints; is_endpoint[1] = true;}}
					else{            if(HH_thermal[-element[1]].is_strendpt()){++numendpoints; is_endpoint[1] = true;}}
					if(perm2[q3]>0){if(showerquarks[element[2]].is_strendpt()){++numendpoints; is_endpoint[2] = true;}}
					else{            if(HH_thermal[-element[2]].is_strendpt()){++numendpoints; is_endpoint[2] = true;}}
					
					//setting the current string value to the value of the first parton (since it can't be a thermal parton)
					//this might have to change if partons are allowed to recombine from multiple strings...
					int cur_str = showerquarks[element[0]].string_id();
					
					if(numendpoints == 0){//no parton is an endpoint
						//finding the string id for this new string...
						//start by choosing an appropriate string 'index'
						int new_str = 100*cur_str + 1;
						//looping over all the current string indices to determine if the first guess is unique; if not, then increment it until it is unique
						//this might loop multiple times if the chosen index is not unique
						{int i=0; while(i<list_strs.size()){if(new_str == list_strs[i]){++new_str; i=-1;} ++i;}}
						//need to add the new string to the list of strings, unless for some reason we want to keep adding the 'new' strings together?
						list_strs.push_back(new_str);
						
						showerquarks[showerquarks[element[0]].sibling()].string_id(new_str);
						showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
						showerquarks[showerquarks[element[0]].sibling()].pos_str(0);
						showerquarks[showerquarks[element[0]].sibling()].endpt_id(1);
						if(perm2[q2]>0){
							showerquarks[showerquarks[element[1]].sibling()].string_id(new_str);
							showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
							showerquarks[showerquarks[element[1]].sibling()].pos_str(0);
							showerquarks[showerquarks[element[1]].sibling()].endpt_id(1);
						}
						else{
							int thermsib = findthermalsibling(-element[1], HH_thermal);
							HH_thermal[thermsib].string_id(new_str);
							HH_thermal[thermsib].is_strendpt(true);
							HH_thermal[thermsib].pos_str(0);
							HH_thermal[thermsib].endpt_id(1);
						}
						if(perm2[q3]>0){
							showerquarks[showerquarks[element[2]].sibling()].string_id(new_str);
							showerquarks[showerquarks[element[2]].sibling()].is_strendpt(true);
							showerquarks[showerquarks[element[2]].sibling()].pos_str(0);
							showerquarks[showerquarks[element[2]].sibling()].endpt_id(1);
						}
						else{
							int thermsib = findthermalsibling(-element[2], HH_thermal);
							HH_thermal[thermsib].string_id(new_str);
							HH_thermal[thermsib].is_strendpt(true);
							HH_thermal[thermsib].pos_str(0);
							HH_thermal[thermsib].endpt_id(1);
						}
					}
					else if(numendpoints == 1){//only one parton is an endpoint
						if(is_endpoint[0]){
							if(perm2[q2]>0){
								showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[1]].sibling()].pos_str(showerquarks[element[0]].pos_str());
								showerquarks[showerquarks[element[1]].sibling()].endpt_id(showerquarks[element[0]].endpt_id());
							}
							else{
								int thermsib = findthermalsibling(-element[1], HH_thermal);
								HH_thermal[thermsib].string_id(cur_str);
								HH_thermal[thermsib].is_strendpt(true);
								HH_thermal[thermsib].pos_str(showerquarks[element[0]].pos_str());
								HH_thermal[thermsib].endpt_id(showerquarks[element[0]].endpt_id());
							}
							if(perm2[q3]>0){
								showerquarks[showerquarks[element[2]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[2]].sibling()].pos_str(showerquarks[element[0]].pos_str());
								showerquarks[showerquarks[element[2]].sibling()].endpt_id(showerquarks[element[0]].endpt_id());
							}
							else{
								int thermsib = findthermalsibling(-element[2], HH_thermal);
								HH_thermal[thermsib].string_id(cur_str);
								HH_thermal[thermsib].is_strendpt(true);
								HH_thermal[thermsib].pos_str(showerquarks[element[0]].pos_str());
								HH_thermal[thermsib].endpt_id(showerquarks[element[0]].endpt_id());
							}
						}
						else if(is_endpoint[1]){
							if(perm2[q2]>0){
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(showerquarks[element[1]].pos_str());
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(showerquarks[element[1]].endpt_id());
								if(perm2[q3]>0){
									showerquarks[showerquarks[element[2]].sibling()].is_strendpt(true);
									showerquarks[showerquarks[element[2]].sibling()].pos_str(showerquarks[element[1]].pos_str());
									showerquarks[showerquarks[element[2]].sibling()].endpt_id(showerquarks[element[1]].endpt_id());
								}
								else{
									int thermsib = findthermalsibling(-element[2], HH_thermal);
									HH_thermal[thermsib].string_id(cur_str);
									HH_thermal[thermsib].is_strendpt(true);
									HH_thermal[thermsib].pos_str(showerquarks[element[1]].pos_str());
									HH_thermal[thermsib].endpt_id(showerquarks[element[1]].endpt_id());
								}
							}
							else{//these shouldn't actually ever trigger, as thermal partons should never 'be' an endpoint...
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(HH_thermal[-element[1]].pos_str());
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(HH_thermal[-element[1]].endpt_id());
								if(perm2[q3]>0){
									showerquarks[showerquarks[element[2]].sibling()].is_strendpt(true);
									showerquarks[showerquarks[element[2]].sibling()].pos_str(HH_thermal[-element[1]].pos_str());
									showerquarks[showerquarks[element[2]].sibling()].endpt_id(HH_thermal[-element[1]].endpt_id());
								}
								else{
									int thermsib = findthermalsibling(-element[2], HH_thermal);
									HH_thermal[thermsib].string_id(cur_str);
									HH_thermal[thermsib].is_strendpt(true);
									HH_thermal[thermsib].pos_str(HH_thermal[-element[1]].pos_str());
									HH_thermal[thermsib].endpt_id(HH_thermal[-element[1]].endpt_id());
								}
							}
						}
						else{//is_endpoint[2]
							if(perm2[q3]>0){
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(showerquarks[element[2]].pos_str());
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(showerquarks[element[2]].endpt_id());
								if(perm2[q2]>0){
									showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
									showerquarks[showerquarks[element[1]].sibling()].pos_str(showerquarks[element[2]].pos_str());
									showerquarks[showerquarks[element[1]].sibling()].endpt_id(showerquarks[element[2]].endpt_id());
								}
								else{
									int thermsib = findthermalsibling(-element[1], HH_thermal);
									HH_thermal[thermsib].string_id(cur_str);
									HH_thermal[thermsib].is_strendpt(true);
									HH_thermal[thermsib].pos_str(showerquarks[element[2]].pos_str());
									HH_thermal[thermsib].endpt_id(showerquarks[element[2]].endpt_id());
								}
							}
							else{//these shouldn't actually ever trigger, as thermal partons should never 'be' an endpoint...
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(HH_thermal[-element[2]].pos_str());
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(HH_thermal[-element[2]].endpt_id());
								if(perm2[q2]>0){
									showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
									showerquarks[showerquarks[element[1]].sibling()].pos_str(HH_thermal[-element[2]].pos_str());
									showerquarks[showerquarks[element[1]].sibling()].endpt_id(HH_thermal[-element[2]].endpt_id());
								}
								else{
									int thermsib = findthermalsibling(-element[1], HH_thermal);
									HH_thermal[thermsib].string_id(cur_str);
									HH_thermal[thermsib].is_strendpt(true);
									HH_thermal[thermsib].pos_str(HH_thermal[-element[2]].pos_str());
									HH_thermal[thermsib].endpt_id(HH_thermal[-element[2]].endpt_id());
								}
							}
						}
					}
					else if(numendpoints == 2){//there are two partons that are endpoints
						if(!is_endpoint[0]){
							//replacing the closest endpoint with the sibling of the non-endpoint parton
							//since neither of the endpoints can be thermal partons, don't need to check
							int set_pos, set_endptid;
							if(std::abs(showerquarks[element[0]].pos_str()-showerquarks[element[1]].pos_str()) <= std::abs(showerquarks[element[0]].pos_str()-showerquarks[element[2]].pos_str())){
								set_pos = showerquarks[element[1]].pos_str(); set_endptid = showerquarks[element[1]].endpt_id();
							}
							else{set_pos = showerquarks[element[2]].pos_str(); set_endptid = showerquarks[element[2]].endpt_id();}
							showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
							showerquarks[showerquarks[element[0]].sibling()].pos_str(set_pos);
							showerquarks[showerquarks[element[0]].sibling()].endpt_id(set_endptid);
						}
						else if(!is_endpoint[1]){
							int set_pos, set_endptid;
							if(std::abs(showerquarks[element[1]].pos_str()-showerquarks[element[0]].pos_str()) <= std::abs(showerquarks[element[1]].pos_str()-showerquarks[element[2]].pos_str())){
								set_pos = showerquarks[element[0]].pos_str(); set_endptid = showerquarks[element[0]].endpt_id();
							}
							else{set_pos = showerquarks[element[2]].pos_str(); set_endptid = showerquarks[element[2]].endpt_id();}
							if(perm2[q2]>0){
								showerquarks[showerquarks[element[1]].sibling()].string_id(cur_str);
								showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[1]].sibling()].pos_str(set_pos);
								showerquarks[showerquarks[element[1]].sibling()].endpt_id(set_endptid);
							}
							else{
								int thermsib = findthermalsibling(-element[1], HH_thermal);
								HH_thermal[thermsib].string_id(cur_str);
								HH_thermal[thermsib].is_strendpt(true);
								HH_thermal[thermsib].pos_str(set_pos);
								HH_thermal[thermsib].endpt_id(set_endptid);
							}
						}
						else{//is_endpoint[2]
							int set_pos, set_endptid;
							if(std::abs(showerquarks[element[2]].pos_str()-showerquarks[element[0]].pos_str()) <= std::abs(showerquarks[element[2]].pos_str()-showerquarks[element[1]].pos_str())){
								set_pos = showerquarks[element[0]].pos_str(); set_endptid = showerquarks[element[0]].endpt_id();
							}
							else{set_pos = showerquarks[element[1]].pos_str(); set_endptid = showerquarks[element[1]].endpt_id();}
							if(perm2[q3]>0){
								showerquarks[showerquarks[element[2]].sibling()].string_id(cur_str);
								showerquarks[showerquarks[element[2]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[2]].sibling()].pos_str(set_pos);
								showerquarks[showerquarks[element[2]].sibling()].endpt_id(set_endptid);
							}
							else{
								int thermsib = findthermalsibling(-element[2], HH_thermal);
								HH_thermal[thermsib].string_id(cur_str);
								HH_thermal[thermsib].is_strendpt(true);
								HH_thermal[thermsib].pos_str(set_pos);
								HH_thermal[thermsib].endpt_id(set_endptid);
							}
						}
					}
					//else{//all three partons are endpoints
					//	//whelp, this string is now a gluon loop?
					//}
					//*******************************************string repair functionality above*******************************************
					
					//now we're forming the hadron
					HHhadron formedhadron;
					//setting the hadron values: is a recombined hadron, mass, and parents
					formedhadron.is_recohad(true); formedhadron.mass( p_BCM[0].t() + p_BCM[1].t() + p_BCM[2].t() );
					formedhadron.add_par(showerquarks[element[0]].par());
					if(perm2[q2]>0){formedhadron.add_par(showerquarks[element[1]].par());}else{formedhadron.add_par(element[1]);}
					if(perm2[q3]>0){formedhadron.add_par(showerquarks[element[2]].par());}else{formedhadron.add_par(element[2]);}
					
					//now setting if baryon is in excited state
					if(WigB[0]*recofactor3 < rndbaryon){formedhadron.is_excited(true);}
					
					//setting if there are any thermal partons used to make the hadron (with the '0' thermal parent as -99999 so that it doesn't conflict with '0' shower parton)
					if(     perm2[q2]>0 && perm2[q3]>0){/*is sh-sh-sh*/ formedhadron.is_shsh(true);}
					else if(perm2[q2]>0 && perm2[q3]<0){/*is sh-sh-th*/ formedhadron.is_shth(true); if(element[2] == 0){formedhadron.parents[2] = -99999;}}
					else if(perm2[q2]<0 && perm2[q3]>0){/*is sh-th-sh*/ formedhadron.is_shth(true); if(element[1] == 0){formedhadron.parents[1] = -99999;}}
					else if(perm2[q2]<0 && perm2[q3]<0){/*is sh-th-th*/ formedhadron.is_shth(true); if(element[1] == 0){formedhadron.parents[1] = -99999;}
						                                                                            if(element[2] == 0){formedhadron.parents[2] = -99999;}
					}
					
					//setting hadron position and momentum vectors
					Pbaryon.Set(Pbaryon.x(),Pbaryon.y(),Pbaryon.z(),sqrt(Pbaryon.x()*Pbaryon.x() + Pbaryon.y()*Pbaryon.y() + Pbaryon.z()*Pbaryon.z() + formedhadron.mass()*formedhadron.mass()));
					formedhadron.pos(pos_lab); formedhadron.P(Pbaryon);
					
					//need to choose *what* hadron we've formed... base this on the parton id's, mass, & if excited
					//might want to do this differently? void f'n(partoncollection, formedhadron)?
					set_baryon_id(considering, formedhadron);
					
					//need to add the hadron to the collection
					HH_hadrons.add(formedhadron);
					
					//now that we've formed the hadron, need to set ALL (3) the 'considering' flags to used
					                showerquarks[element[0]].status(1); showerquarks[element[0]].is_used(true);
					if(perm2[q2]>0){showerquarks[element[1]].status(1); showerquarks[element[1]].is_used(true);}
					else{            HH_thermal[-element[1]].status(1);  HH_thermal[-element[1]].is_used(true);}
					if(perm2[q3]>0){showerquarks[element[2]].status(1); showerquarks[element[2]].is_used(true);}
					else{            HH_thermal[-element[2]].status(1);  HH_thermal[-element[2]].is_used(true);}
					
					madehadron = true; considering.clear();
					break;
				}
				else{
					//since we've not formed a baryon on this try, need to revert the third quark 'considering' used flag
					if(perm2[q3]>0){showerquarks[element[2]].status(0);}
					else{            HH_thermal[-element[2]].status(0);}
					
					//and remove the third entry in considering
					considering.partons.pop_back(); //considering--;
				}
				
			}}
			else if(considering[0].id()*considering[1].id() < 0){
				//now that we *could* form a meson, now we check if we actually do form one
				//meson momentum
				FourVector Pmeson;
				Pmeson.Set(considering[0].px()+considering[1].px(),considering[0].py()+considering[1].py(),considering[0].pz()+considering[1].pz(),0.);
				
				//meson(CM) velocity
				FourVector betaM;
				betaM.Set(Pmeson.x()/(considering[0].e()+considering[1].e()),Pmeson.y()/(considering[0].e()+considering[1].e()),Pmeson.z()/(considering[0].e()+considering[1].e()),0.);
				betaM.Set(betaM.x(),betaM.y(),betaM.z(),1./(sqrt(1.-(betaM.x()*betaM.x() + betaM.y()*betaM.y() + betaM.z()*betaM.z()))));
				
				//boosting into CM frame
				FourVector pos_MCM[2], p_MCM[2];
				pos_MCM[0] = considering[0].boost_pos(betaM); pos_MCM[1] = considering[1].boost_pos(betaM);
				  p_MCM[0] = considering[0].boost_P(betaM);     p_MCM[1] = considering[1].boost_P(betaM);
				
				//velocities in CM frame
				FourVector v_MCM[2];
				v_MCM[0].Set(p_MCM[0].x()/p_MCM[0].t(),p_MCM[0].y()/p_MCM[0].t(),p_MCM[0].z()/p_MCM[0].t(),0.);
				v_MCM[1].Set(p_MCM[1].x()/p_MCM[1].t(),p_MCM[1].y()/p_MCM[1].t(),p_MCM[1].z()/p_MCM[1].t(),0.);
				
				//propagating quarks until time of youngest quark
				//is just max(pos_MCM[0].t(), pos_MCM[1].t());
				double curtime = (pos_MCM[0].t() > pos_MCM[1].t()) ? pos_MCM[0].t() : pos_MCM[1].t();
				FourVector cur_pos[2];
				cur_pos[0].Set(pos_MCM[0].x()+v_MCM[0].x()*(curtime-pos_MCM[0].t()),pos_MCM[0].y()+v_MCM[0].y()*(curtime-pos_MCM[0].t()),pos_MCM[0].z()+v_MCM[0].z()*(curtime-pos_MCM[0].t()),curtime);
				cur_pos[0].Set(pos_MCM[1].x()+v_MCM[1].x()*(curtime-pos_MCM[1].t()),pos_MCM[1].y()+v_MCM[1].y()*(curtime-pos_MCM[1].t()),pos_MCM[1].z()+v_MCM[1].z()*(curtime-pos_MCM[1].t()),curtime);
				
				//finding position of CM at curtime
				FourVector pos_CM;
				pos_CM.Set(
				(cur_pos[0].x()*considering[0].mass()+cur_pos[1].x()*considering[1].mass())/(considering[0].mass()+considering[1].mass()),
				(cur_pos[0].y()*considering[0].mass()+cur_pos[1].y()*considering[1].mass())/(considering[0].mass()+considering[1].mass()),
				(cur_pos[0].z()*considering[0].mass()+cur_pos[1].z()*considering[1].mass())/(considering[0].mass()+considering[1].mass()),
				curtime);
				
				//finding position of meson in lab frame
				betaM.Set(-betaM.x(),-betaM.y(),-betaM.z(),betaM.t());
				FourVector pos_lab = HHboost(betaM, pos_CM);
				
				//finding relative momenta of partons in CM frame
				FourVector k_rel;
				k_rel.Set((p_MCM[1].x()-p_MCM[0].x())*(p_MCM[1].x()-p_MCM[0].x())/4.,(p_MCM[1].y()-p_MCM[0].y())*(p_MCM[1].y()-p_MCM[0].y())/4.,(p_MCM[1].z()-p_MCM[0].z())*(p_MCM[1].z()-p_MCM[0].z())/4.,0.);
				k_rel.Set(k_rel.x(),k_rel.y(),k_rel.z(),k_rel.x()+k_rel.y()+k_rel.z());
				
				//finding relative positions of partons in CM frame
				FourVector pos_rel;
				pos_rel.Set((cur_pos[0].x()-cur_pos[1].x())*(cur_pos[0].x()-cur_pos[1].x()),(cur_pos[0].y()-cur_pos[1].y())*(cur_pos[0].y()-cur_pos[1].y()),(cur_pos[0].z()-cur_pos[1].z())*(cur_pos[0].z()-cur_pos[1].z()),0.);
				pos_rel.Set(pos_rel.x(),pos_rel.y(),pos_rel.z(),pos_rel.x()+pos_rel.y()+pos_rel.z());
				
				//setting appropriate sigma...
				double SigM2 = SigPi2;
				int sortid[2] = {0,0};
				if(std::abs(considering[0].id()) >= std::abs(considering[1].id())){sortid[0] = std::abs(considering[0].id()); sortid[1] = std::abs(considering[1].id());}
				else{sortid[0] = std::abs(considering[1].id()); sortid[1] = std::abs(considering[0].id());}
				
				if(     sortid[0] == 3){
					if(     sortid[1] == 3){SigM2 = SigPhi2;}
					else{                   SigM2 = SigK2;}
				}
				else if(sortid[0] == 4){
					if(     sortid[1] == 4){SigM2 = SigJpi2;}
					else if(sortid[1] == 3){SigM2 = SigDs2;}
					else{                   SigM2 = SigD2;}
				}
				else if(sortid[0] == 5){
					if(     sortid[1] == 5){SigM2 = SigUps2;}
					else if(sortid[1] == 4){SigM2 = SigBc2;}
					else if(sortid[1] == 3){SigM2 = SigB2;}
					else{                   SigM2 = SigB2;}
				}
				
				//Calc'ing Wig. wavefunction
				double WigM[2]; WigM[1] = 0.;//sumWigM;
				//3D GS Wig. wavefunction
				WigM[0] = std::exp(-pos_rel.t()/(2.*SigM2) - k_rel.t()*SigM2/hbarc2);
				
				//summing up 3D Wig. wavefunctions over maxE_level excited states
				double u[4];
				u[1] = 0.5*(pos_rel.x()/SigM2 + k_rel.x()*SigM2/hbarc2);
				u[2] = 0.5*(pos_rel.y()/SigM2 + k_rel.y()*SigM2/hbarc2);
				u[3] = 0.5*(pos_rel.z()/SigM2 + k_rel.z()*SigM2/hbarc2);
				u[0] = u[1] + u[2] + u[3];
				
				//this will fail if iME is too large (but why would you want anything quite that big?)
				for(int iME=0; iME<=maxE_level; ++iME){WigM[1] += std::pow(u[0],iME)*std::exp(-u[0])/double(std::tgamma(iME+1));} //std::tgamma(iME+1)==factorial(N)
				
				//Checking if meson is formed (either ground or excited state)
				double rndmeson = ran();
				
				if(WigM[1]*recofactor2 >= rndmeson){
					//*******************************************string repair functionality below*******************************************
					//determining which partons, of the 2 being considered, are endpoints
					int numendpoints = 0; bool is_endpoint[2];  is_endpoint[0]=false;  is_endpoint[1]=false;
									if(showerquarks[element[0]].is_strendpt()){++numendpoints; is_endpoint[0] = true;}
					if(perm2[q2]>0){if(showerquarks[element[1]].is_strendpt()){++numendpoints; is_endpoint[1] = true;}}
					else{            if(HH_thermal[-element[1]].is_strendpt()){++numendpoints; is_endpoint[1] = true;}}
					
					//setting the current string value to the value of the first parton (since it can't be a thermal parton)
					//this might have to change if partons are allowed to recombine from multiple strings...
					int cur_str = showerquarks[element[0]].string_id();
					
					if(numendpoints == 0){//no parton is an endpoint
						//finding the string id for this new string...
						//start by choosing an appropriate string 'index'
						int new_str = 100*cur_str + 1;
						//looping over all the current string indices to determine if the first guess is unique; if not, then increment it until it is unique
						//this might loop multiple times if the chosen index is not unique
						{int i=0; while(i<list_strs.size()){if(new_str == list_strs[i]){++new_str; i=-1;} ++i;}}
						//need to add the new string to the list of strings, unless for some reason we want to keep adding the 'new' strings together?
						list_strs.push_back(new_str);
						
						showerquarks[showerquarks[element[0]].sibling()].string_id(new_str);
						showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
						showerquarks[showerquarks[element[0]].sibling()].pos_str(0);
						showerquarks[showerquarks[element[0]].sibling()].endpt_id(1);
						if(perm2[q2]>0){
							showerquarks[showerquarks[element[1]].sibling()].string_id(new_str);
							showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
							showerquarks[showerquarks[element[1]].sibling()].pos_str(0);
							showerquarks[showerquarks[element[1]].sibling()].endpt_id(1);
						}
						else{
							int thermsib = findthermalsibling(-element[1], HH_thermal);
							HH_thermal[thermsib].string_id(new_str);
							HH_thermal[thermsib].is_strendpt(true);
							HH_thermal[thermsib].pos_str(0);
							HH_thermal[thermsib].endpt_id(1);
						}
					}
					else if(numendpoints == 1){//only one parton is an endpoint
						if(is_endpoint[0]){
							if(perm2[q2]>0){
								showerquarks[showerquarks[element[1]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[1]].sibling()].pos_str(showerquarks[element[0]].pos_str());
								showerquarks[showerquarks[element[1]].sibling()].endpt_id(showerquarks[element[0]].endpt_id());
							}
							else{
								int thermsib = findthermalsibling(-element[1], HH_thermal);
								HH_thermal[thermsib].string_id(cur_str);
								HH_thermal[thermsib].is_strendpt(true);
								HH_thermal[thermsib].pos_str(showerquarks[element[0]].pos_str());
								HH_thermal[thermsib].endpt_id(showerquarks[element[0]].endpt_id());
							}
						}
						else if(is_endpoint[1]){
							if(perm2[q2]>0){
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(showerquarks[element[1]].pos_str());
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(showerquarks[element[1]].endpt_id());
							}
							else{//these shouldn't actually ever trigger, as thermal partons should never 'be' an endpoint...
								showerquarks[showerquarks[element[0]].sibling()].is_strendpt(true);
								showerquarks[showerquarks[element[0]].sibling()].pos_str(0); //HH_thermal[-element[1]].pos_str;
								showerquarks[showerquarks[element[0]].sibling()].endpt_id(1); //HH_thermal[-element[1]].endpt_id;
							}
						}
					}
					//else{//both endpoints are endpoints
					//	//whelp, this string is now a gluon loop?
					//}
					//*******************************************string repair functionality above*******************************************
					
					//now we're forming the hadron
					HHhadron formedhadron;
					//setting the hadron values: is a recombined hadron, mass, and parents - setting the par3 flag to 999999 to denote the lack of a 3rd parent parton
					formedhadron.is_recohad(true); formedhadron.mass( p_MCM[0].t() + p_MCM[1].t() );
					//formedhadron.par1 = element[0]; formedhadron.par2 = element[1]; formedhadron.par3 = 999999;
					formedhadron.add_par(showerquarks[element[0]].par());
					if(perm2[q2]>0){formedhadron.add_par(showerquarks[element[1]].par());}else{formedhadron.add_par(element[1]);}
					
					//now setting if baryon is in excited state
					if(WigM[0]*recofactor2 < rndmeson){formedhadron.is_excited(true);}
					
					//setting if there are any thermal partons used to make the hadron: (setting the '0' thermal parent to -99999 so that it doesn't conflict with '0' shower parton)
					if(     perm2[q2]>0){/*is sh-sh*/ formedhadron.is_shsh(true);}
					else if(perm2[q2]<0){/*is sh-th*/ formedhadron.is_shth(true); if(element[1] == 0){formedhadron.parents[1] = -99999;}}
					
					//setting hadron position and momentum vectors
					Pmeson.Set(Pmeson.x(),Pmeson.y(),Pmeson.z(),sqrt(Pmeson.x()*Pmeson.x() + Pmeson.y()*Pmeson.y() + Pmeson.z()*Pmeson.z() + formedhadron.mass()*formedhadron.mass()));
					formedhadron.pos(pos_lab); formedhadron.P(Pmeson);
					
					//need to choose *what* hadron we've formed... base this on the parton id's, mass, & if excited
					set_meson_id(considering, formedhadron);
					
					//need to add the hadron to the collection
					HH_hadrons.add(formedhadron);
					
					//now that we've formed the hadron, need to set ALL (both) the 'considering' flags to used
					                showerquarks[element[0]].status(1); showerquarks[element[0]].is_used(true);
					if(perm2[q2]>0){showerquarks[element[1]].status(1); showerquarks[element[1]].is_used(true);}
					else{            HH_thermal[-element[1]].status(1);  HH_thermal[-element[1]].is_used(true);}
					
					//now that we've formed the hadron, break to first loop here!
					madehadron = true; considering.clear();
					break;
				}
			}
			
			//if we've formed a hadron - break to first loop
			if(madehadron){break;}
			
			//since we CAN'T form a baryon on this try, need to revert the second quark 'considering' used flag
			if(perm2[q2]>0){showerquarks[element[1]].status(0);}
			else{            HH_thermal[-element[1]].status(0);}
			
			//and remove the second entry in considering (putting in if statement in case we tried and failed to make a hadron;
			considering.partons.pop_back(); //considering--;
		}
		//if we've formed a hadron - continue to next parton in first loop...
		if(madehadron){continue;}
		
		//since we've not formed a hadron with the first quark, need to revert the first quark 'considering' used flag
		//only need to reset these for shower quarks as the first quark cannot be a thermal quark
		//and remove the first(only) entry in considering
		showerquarks[element[0]].status(0); considering.partons.pop_back(); //considering--;
	}
	
	//all possibilities have been considered, all used quark flags for showerquarks are set appropriately
	//time for cleanup
	
	//set the used quarks in the original shower to reflect that they were used in reco module and that they were actually used
	//set the fully used gluons in the original shower to reflect that they were used in reco module and that they were actually used
	//set the partially used gluons in the original shower to reflect that they were used in reco module and that they were actually used
	//give used quarks and fully used gluons a status of '1'; give partially used gluons a status of '-1'
	//stick all unused quarks and 'completely' unused gluons into remnants (make sure to set the parents to the original shower partons appropriately)
	//write updated string information into shower, so that it is set properly in remnants
	//for the thermal array, all the used flags should already be set (and have a status of '1')
	//if this is used in a loop (as the original version should be doing) - those will have to be reset before reuse
	//otherwise, this is perfect for medium feedback
	
	//using quarks in showerquarks to set the flags appropriately for partons in shower
	for(int i=0; i<showerquarks.num(); ++i){
		//if we have a quark in the original shower
		if((std::abs(HH_showerptns[showerquarks[i].par()].id()) <= 5) && (showerquarks[i].is_used())){
			HH_showerptns[showerquarks[i].par()].is_used(true); HH_showerptns[showerquarks[i].par()].status(1); HH_showerptns[showerquarks[i].par()].used_reco(true);
		}
		//if this quark is from a split gluon in the original shower
		else if(std::abs(HH_showerptns[showerquarks[i].par()].id()) == 21 && showerquarks[i].is_used()){
			//if this is the first used quark in a splitting, set the parent gluon to used, and status of -1
			if(HH_showerptns[showerquarks[i].par()].status() == -99){
				HH_showerptns[showerquarks[i].par()].status(-1); HH_showerptns[showerquarks[i].par()].is_used(true); HH_showerptns[showerquarks[i].par()].used_reco(true);
			}
			//if this is the second (last) used quark in a splitting, set the status to 1
			else if(HH_showerptns[showerquarks[i].par()].status() == -1){HH_showerptns[showerquarks[i].par()].status(1);}
			//remove this check if it never throws.
			else{JSWARN << "SOMETHING HAS GONE VERY WRONG WITH REFORMING GLUON IN POS: " << showerquarks[i].par(); int val; std::cin >> val; showerquarks[i].par(val);}
		}
	}
	//need to run back through shower; if there are any gluons that didn't get used at all in the shower - restore them (and output to remnants below)
	for(int i=0; i<HH_showerptns.num(); ++i){if(HH_showerptns[i].status() == -99){HH_showerptns[i].is_decayedglu(false); HH_showerptns[i].status(0);}}
	
	//need to update string information in shower from showerquarks
	for(int i=0; i<HH_showerptns.num(); ++i){if(!HH_showerptns[i].is_used()){
		for(int j=0; j<showerquarks.num(); ++j){if(showerquarks[j].par() == i && !showerquarks[j].is_used()){
			HH_showerptns[i].string_id(   showerquarks[j].string_id());
			HH_showerptns[i].is_strendpt( showerquarks[j].is_strendpt());
			HH_showerptns[i].pos_str(     showerquarks[j].pos_str());
			HH_showerptns[i].endpt_id(    showerquarks[j].endpt_id());
			break;
		}}
	}}
	//now all the partons in shower have had flags appropriately set; the 'partially' used gluons have a status of -1
	
	//sticking all unused partons into remnants; keeping order intact (a partially used gluon is replaced with it's unused quark)
	for(int i=0; i<HH_showerptns.num(); ++i){
		//if unused parton, write into remnants
		if(HH_showerptns[i].status() == 0){HH_remnants.add(HH_showerptns[i]); HH_remnants[HH_remnants.num() - 1].par(i); HH_showerptns[i].is_remnant(true);}
		//if 'partially' used gluon, write unused daughter quark into remnants
		else if(HH_showerptns[i].status() == -1){
			//finding the unused quark for this gluon and adding it to remnants (have to loop over as we only keep track of parents, not daughters)
			for(int j=0; j<showerquarks.num(); ++j){if(showerquarks[j].par() == i && !showerquarks[j].is_used()){HH_remnants.add(showerquarks[j]); break;}}
			HH_showerptns[i].is_remnant(true);
		}
	}
	
	//appending the thermal partons used in the string repair functionality into remnants - order is NOT preserved...
	//is later sorted based not only on the string, but also on the position of the parton IN the string
	//can use the fact that the thermal partons needed either will be endpoints, or will be in a string with id > 0
	for(int i=0; i<HH_thermal.num(); ++i){
		//if this thermal parton is ... then add it to the remnants collection
		if(HH_thermal[i].is_used()){HH_thermal[i].status(1); HH_thermal[i].used_reco(true);}
		if(HH_thermal[i].is_strendpt()){HH_remnants.add(HH_thermal[i]); HH_remnants[HH_remnants.num() - 1].par(-i-1); HH_thermal[i].is_remnant(true);}
	}
	
	//hadrons have been recombined, and output to hadron collection
	//remnants have been collected, and output to remnant collection
	//shower partons have all been updated appropriately
	//thermal partons have all been updated appropriately - make sure that the thermal partons are reset before the recomb module is called again...
	
	//Future plans:
	//since we've already read in the thermal parton array, can use that to get the 'fake' parton for the necessary string, if present?
	//include color after completion?
	
//end of recombination routine
}

//sets id of formed baryon based on quark content, mass of quark system, and if the baryon formed into an excited state
void HybridHadronization::set_baryon_id(parton_collection& qrks, HHhadron& had){
	
	//assigning quark_ids in descending order to construct baryon id
	int id[3] = {qrks[0].id(), qrks[1].id(), qrks[2].id()};
	if (std::abs(id[0]) < std::abs(id[1])){std::swap(id[0],id[1]);}
	if (std::abs(id[1]) < std::abs(id[2])){std::swap(id[1],id[2]);}
	if (std::abs(id[0]) < std::abs(id[1])){std::swap(id[0],id[1]);}
	
/*	//http://pdg.lbl.gov/2017/listings/contents_listings.html
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
	//id and mass list:
	//ground state: (111, 222, 333 configurations prohibited...)
	//	2212, 2112 - .938  p,n
	//	3222, 3212, 3112 - 1.190   sigma;  3122 - 1.115 lambda
	//	3322, 3312 - 1.315  xi
	//
	//excited:
	//	2224, 2214, 2114, 1114 - 1.232  delta
	//	3224, 3214, 3114 - 1.190   sigma;
	//	3324, 3314 - 1.315  xi
	//	3334 - 1.672 omega
//	^^^^^ MIGHT NOT NEED BUT KEEP FOR NOW ^^^^^^
	
	if(id[0]==id[1] && id[0]==id[2]){had.id( 1000*std::abs(id[0]) + 100*std::abs(id[1]) + 10*std::abs(id[2])+4 );}  // J=3/2 only
	else if(id[0]==id[1] || id[0]==id[2] || id[1]==id[2]){
		if(ran()>0.333){had.id( 1000*std::abs(id[0]) + 100*std::abs(id[1]) + 10*std::abs(id[2])+4 );}  // J=3/2 
		else{           had.id( 1000*std::abs(id[0]) + 100*std::abs(id[1]) + 10*std::abs(id[2])+2 );}
	}
	else{
		double prb = ran();
		if(     prb>0.333){had.id( 1000*std::abs(id[0]) + 100*std::abs(id[1]) + 10*std::abs(id[2])+4 );}  // J=3/2
		else if(prb>0.166){had.id( 1000*std::abs(id[0]) + 100*std::abs(id[1]) + 10*std::abs(id[2])+2 );}  // J=1/2 higher mass
		else{              had.id( 1000*std::abs(id[0]) + 10*std::abs(id[1]) + 100*std::abs(id[2])+2 );}  // J=1/2 lower mass
	}	// note the swap of quark index in last line
	
	// This would be how to make excited N, Delta and Lambda. RIGHT NOW DISABLED IN FIRST IF STATEMENT
	//there are no excited baryon codes of this form; there are only spin excited states in PYTHIA's ParticleData class
//	if(false && had.is_excited){
//		if(had.id == 2212 || had.id == 2112 || had.id == 2214 || had.id == 2114 || had.id == 2224 || had.id ==1114 || had.id == 3122){
//			had.id += 100000;  // This needs to be adjusted to make the actual excited baryon codes
//		}
//	}
	
	had.id( had.id()*(2*std::signbit(-id[0])-1) );
	return;
}

//sets id of formed meson based on quark content, mass of quark system, and if the meson formed into an excited state
void HybridHadronization::set_meson_id(parton_collection& qrks, HHhadron& had){
	
	//assigning quark_ids in descending order to construct meson id
	int id[2] = {qrks[0].id(), qrks[1].id()};
	if(std::abs(qrks[1].id()) > std::abs(qrks[0].id())){id[0] = qrks[1].id(); id[1] = qrks[0].id();}
	
//this whole thing needs to be updated when spin/excited state n is included - will make for better physics
//it probably still needs to be updated!
//	double mass_pi0, mass_eta, mass_omegam, mass_etap, mass_phi; //mass_rho;
//	mass_pi0    = 0.1349770;
//	mass_eta    = 0.547862;
//	mass_rho    = 0.7690;
//	mass_omegam = 0.78265;
//	mass_etap   = 0.95778;
//	mass_phi    = 1.019460;
//// MASSES ABOVE NOT REALLY NEEDED IN SIMPLE QUARK MODEL BASED APPROACH	
	
	int baseid = 0;
	//if isospin I3=0
	if(id[0] == -id[1]){
		if(ran() > 0.25){	//spin triplet
			if(std::abs(id[0])==5){baseid = 553;}						// Upsilon
			if(std::abs(id[0])==4){baseid = 443;}					// J/psi
			if(std::abs(id[0])==3){baseid = 333;}					// phi
			if((std::abs(id[0])==1) || (std::abs(id[0]) == 2)){
				if(ran()>0.5){baseid = 223;} else{baseid = 113;}
			}	// omega and rho
		}	
		else{				// spin singlet
			if(std::abs(id[0])==5){baseid = 551; }						// etaB
			if(std::abs(id[0])==4){baseid = 441; }					// etac
			if(std::abs(id[0])==3){if(ran()>0.666){baseid = 331;} else{baseid = 221;}}	// eta' and eta
			if(std::abs(id[0])<3){
				double prb = ran();
				if(prb>0.5){baseid = 111;}							// pi0
				else if(prb>0.333){baseid = 221;}					// eta	
				else{baseid = 331;}									// eta'
			}
		}
		//if(had.is_excited){baseid += 100000;}
	}
	// if isospin I3 not 0
	else{
		int basesign = 1;
		baseid = 100*std::abs(id[0])+10*std::abs(id[1]);
		if(id[0]%2 == 0){basesign = 2*std::signbit(-id[0])-1;}
		else{            basesign = 2*std::signbit( id[0])-1;}
		
		//if(had.is_excited){baseid += 100000;}
		if(ran() > 0.25){baseid += 3;}
		else{baseid += 1;}
		
		baseid *= basesign;
	}
	
	had.id( baseid );
	return;
}

//gluon decay function
void HybridHadronization::gluon_decay(HHparton& glu, parton_collection& qrks){
	
	HHparton q1, q2;
	//these are placeholders - might want to instead use values directly from partons...
	double qmass, glu_e;
	
	//if set to already be not on-shell, but not initially set!
	//glu.mass = sqrt(glu.e()*glu.e() - glu.px()*glu.px() - glu.py()*glu.py() - glu.pz()*glu.pz());
	//glu_e = glu.e();
	glu_e = sqrt(glu.mass()*glu.mass() + glu.px()*glu.px() + glu.py()*glu.py() + glu.pz()*glu.pz());
	
	//choosing qqbar ids (u, d, or s)
	//assuming that xms >= xmq (bad things *could* happen if not...)
	if(glu.mass() > 2.*xms){
		//******** ratio = Gamma(g->ssbar)/Gamma(g->uubar, ddbar) ******
		double ratio = 0.5*sqrt((glu.mass()*glu.mass()-4.*xms*xms)/(glu.mass()*glu.mass()-4.*xmq*xmq))*((glu.mass()*glu.mass()+2.*xms*xms)/(glu.mass()*glu.mass()+2.*xmq*xmq));
		double prob = ran();
		if(prob <= ratio/(1.+ratio)){qmass = xms; q1.id(3); q2.id(-3); q1.mass(xms); q2.mass(xms);}
		else if((prob > ratio/(1.+ratio)) && (prob <= (0.5+ratio)/(1.+ratio))){qmass = xmq; q1.id(1); q2.id(-1); q1.mass(xmq); q2.mass(xmq);}
		else{ /*if (prob > (0.5+ratio)/(1.+ratio))*/ qmass = xmq; q1.id(2); q2.id(-2); q1.mass(xmq); q2.mass(xmq);}
	}
	else{
		double prob = ran();
		if(prob <= 0.5){qmass = xmq; q1.id(1); q2.id(-1); q1.mass(xmq); q2.mass(xmq);}
		else{           qmass = xmq; q1.id(2); q2.id(-2); q1.mass(xmq); q2.mass(xmq);}
	}
	
	//gluon velocity
	FourVector Betag;
	Betag.Set(-glu.px()/glu_e,-glu.py()/glu_e,-glu.pz()/glu_e,0.);
	double sum2 = Betag.x()*Betag.x() + Betag.y()*Betag.y() + Betag.z()*Betag.z();
	if(sum2 > 0.){Betag.Set(Betag.x(),Betag.y(),Betag.z(),1./sqrt(sum2));}
	
	//setting the q-qbar position (done as it was in FORTRAN code - might be better to use inbuilt boost somehow?)
	FourVector position;
	double tau = hbarc*Betag.t()/glu.mass();
	position.Set(glu.x()-Betag.x()*tau,glu.y()-Betag.y()*tau,glu.z()-Betag.z()*tau,glu.x_t()+tau);
	q1.pos(position); q2.pos(position);
	
	//setting the q-qbar momenta, starting in gluon rest frame
	FourVector Pq_CM, Pq1, Pq2;
	double pq = (glu.mass() > 2.*qmass) ? sqrt(glu.mass()*glu.mass()/4. - qmass*qmass) : 0.;
	double theta = acos(1.-2.*ran()); double phi = 2*pi*ran();
	Pq_CM.Set(pq*sin(theta)*cos(phi),pq*sin(theta)*sin(phi),pq*cos(theta),sqrt(qmass*qmass+pq*pq));
	Pq1 = HHboost(Betag,Pq_CM);
	Pq_CM.Set(-Pq_CM.x(),-Pq_CM.y(),-Pq_CM.z(),Pq_CM.t());
	Pq2 = HHboost(Betag,Pq_CM);
	q1.P(Pq1); q2.P(Pq2);
	
	qrks.add(q1); qrks.add(q2);
}

//finding a thermal sibling for a thermal parton in therm
int HybridHadronization::findthermalsibling(int ithm, parton_collection& therm){
	if(!therm[therm[ithm].sibling()].is_used() && (therm[therm[ithm].sibling()].string_id() < 0) && (ithm != therm[ithm].sibling())){return therm[ithm].sibling();}
	int qrk_close = -1; double dist2min = 999999999999.;
	for(int i=0;i<therm.num();++i){
		if((therm[ithm].id() * therm[i].id() > 0) || therm[i].is_used()){continue;}
		double distnow = therm[ithm].posDif2(therm[i]) + (therm[ithm].x_t()-therm[i].x_t())*(therm[ithm].x_t()-therm[i].x_t());
		if(distnow < dist2min){qrk_close = i; dist2min = distnow;}
	}
	if(qrk_close == -1){qrk_close = ithm;}
	return qrk_close;
}

//prepares remnant partons/strings for PYTHIA string hadronization
//sorts strings, ensures strings are in 'valid' configurations, assigns color/anticolor tags
//TODO: this might be where to use thermal partons to enforce color neutrality
void HybridHadronization::stringprep(parton_collection& SP_remnants, parton_collection& SP_prepremn, bool cutstr){
	
	//declaring a parton collection to hold a single string in the event, to be 'worked on'
	parton_collection current_str;
	
	//saving the last used color tag id - just increment this up as needed.
	int lastused_tag = 0;
	
	//sort remnants based on string, and position in the string
	std::stable_sort(&SP_remnants[0], (&SP_remnants[SP_remnants.num()-1])+1, [](const HHparton& parton1, const HHparton& parton2){return (parton1.string_id() < parton2.string_id());});
	//now that the list is sorted based on string id, going to sort the partons in each string based on the position of the partons in the string
	for(int i=0; i<SP_remnants.num(); ++i){
		int start, prev_pos, cur_pos, lastfix;
		lastfix = 0; if(i == SP_remnants.num()-1){lastfix = 1;}
		cur_pos = SP_remnants[i].string_id();
		if(i==0){prev_pos = SP_remnants[0].string_id(); start = 0;}
		if(cur_pos != prev_pos || i == SP_remnants.num()-1){
			std::stable_sort(&SP_remnants[start],&SP_remnants[i]+lastfix, [](const HHparton& parton1, const HHparton& parton2){return (parton1.pos_str() < parton2.pos_str()  );});
			start = i;prev_pos = SP_remnants[i].pos_str();
		}
	}
	
	//now that remnants has been sorted properly, need to check the strings to ensure that they're valid.
	//check the number of quarks/antiquarks/(diquarks - though not quite implemented yet)
		//make sure that the numbers of these can form a color singlet (which shouldn't actually break)
		//then make sure that the gluons are in the appropriate position (actually *internal* to the string)
	//lastly, run back through remnants and collect any colorless partons?
	
	int start, prev_str, upd_str;
	for(int par_i=0; par_i<=SP_remnants.num(); ++par_i){
		
		int str_len;
		//need a 'fix' for the first string
		if(par_i==0){prev_str = SP_remnants[0].string_id(); start = 0;}
		//updating the string of the most recent parton
		if(par_i<SP_remnants.num()){upd_str = SP_remnants[par_i].string_id();}else{upd_str = prev_str + 1;}
		//if the current 'parton' is not actually a parton (lepton, hadron, etc.), just dump it in hadrons..?
		//if(!(SP_remnants[par_i].id == 21) && !(std::abs(SP_remnants[par_i].id) <= 6) &&
		// !((std::abs(SP_remnants[par_i].id) >= 1103) && (std::abs(SP_remnants[par_i].id) <= 5503) && ((SP_remnants[par_i].id/10)%10 == 0))){SP_prepremn.add(SP_remnants[par_i]); continue;}
		//if(!(SP_remnants[par_i].id == 21) && !(std::abs(SP_remnants[par_i].id) <= 6) &&
		// !((std::abs(SP_remnants[par_i].id) >= 1103) && (std::abs(SP_remnants[par_i].id) <= 5503) && ((SP_remnants[par_i].id/10)%10 == 0))){continue;}
		//if the current parton is in the same event as the previous, get the next parton - otherwise repair (if necessary) the completed string (start:par_i-1)
		if(upd_str == prev_str){current_str.add(SP_remnants[par_i]); continue;}
		
		//now that we have this complete string, going to check it over to make sure that there's at least one nonthermal parton - otherwise we dump it.
		bool dumpstring = true;
		for(int i=0; i<SP_remnants.num(); ++i){if(!SP_remnants[i].is_thermal()){dumpstring = false; break;}}
		if(dumpstring){
			start = par_i; prev_str = upd_str;
			current_str.clear();
			if(par_i<SP_remnants.num()){current_str.add(SP_remnants[par_i]);}
			continue;
		}
		
		str_len = current_str.num();
		
		//moving on to repair this string...
		//counting up the total number and net number of quarks & diquarks (and the number of gluons)
		int numqrk, netqrk, numdiqrk, netdiqrk, numglu; //numjun; numjun = 0;
		numqrk = 0; netqrk = 0; numdiqrk = 0; netdiqrk = 0; numglu = 0;
		for(int i=0; i<current_str.num(); ++i){
			if(  std::abs(current_str[i].id()) <= 6 && (current_str[i].id() !=0)){++numqrk; netqrk += (2*std::signbit(-current_str[i].id())-1);}
			else if( current_str[i].id()  == 21){++numglu;}
			else if((std::abs(current_str[i].id())>=1103)&&(std::abs(current_str[i].id())<=5503)&&((current_str[i].id()/10)%10==0)){++numdiqrk;netdiqrk+=(4*std::signbit(-current_str[i].id())-2);}
		}
		
		//checking/enforcing color neutrality (will attempt to use thermal partons, if any are present, to accomplish this)
		bool fakepartonadded = false; HHparton fakeparton;
		if((netqrk+netdiqrk)%3 != 0){
			fakepartonadded = true;
			if(    std::abs((netqrk+netdiqrk)%3) == 1){fakeparton.id( -(2*std::signbit(-netqrk-netdiqrk)-1)*(int(1+2*ran())) );}
			else{/*std::abs((netqrk+netdiqrk)%3) == 2*/fakeparton.id(  (2*std::signbit(-netqrk-netdiqrk)-1)*(int(1+2*ran())) );}
			
			int ilast = 0;
			if(std::abs(current_str[current_str.num()-1].id()) == 21){ilast = current_str.num()-1;}
			else if(std::abs(current_str[0].id()) == 21){ilast = 0;}
			else{ilast = (current_str.num()-1)/2;}
			
			if(false){//HH_thermal.num() > 0
				
				int qrk_close = -1; double dist2min = 999999999999.;
				for(int i=0;i<HH_thermal.num();++i){
					if((fakeparton.id() * HH_thermal[i].id() < 0) || HH_thermal[i].is_used()){continue;}
					double distnow = current_str[ilast].posDif2(HH_thermal[i]) + (current_str[ilast].x_t()-HH_thermal[i].x_t())*(current_str[ilast].x_t()-HH_thermal[i].x_t());
					if(distnow < dist2min){qrk_close = i; dist2min = distnow;}
				}
				
				if(qrk_close > -1){
					HH_thermal[qrk_close].string_id(prev_str); HH_thermal[qrk_close].is_remnant(true); HH_thermal[qrk_close].is_used(true); HH_thermal[qrk_close].used_str(true);
					fakeparton.orig(1); fakeparton.is_remnant(true); fakeparton.is_strendpt(true); fakeparton.string_id(prev_str); fakeparton.pos_str(-1);
					fakeparton.px(HH_thermal[qrk_close].px()); fakeparton.py(HH_thermal[qrk_close].py()); fakeparton.pz(HH_thermal[qrk_close].pz()); fakeparton.e(HH_thermal[qrk_close].e());
					fakeparton.x( HH_thermal[qrk_close].x() ); fakeparton.y( HH_thermal[qrk_close].y() ); fakeparton.z( HH_thermal[qrk_close].z() ); fakeparton.x_t(HH_thermal[qrk_close].x_t());
					//will reposition fake parton within the string in the next repair section...
					++str_len; ++numqrk; netqrk += (2*std::signbit(-fakeparton.id())-1);
				}
				else{
					fakeparton.orig(-1); fakeparton.is_remnant(true); fakeparton.is_strendpt(true); fakeparton.string_id(prev_str); fakeparton.pos_str(-1);
					fakeparton.px(0.); fakeparton.py(0.); fakeparton.pz(0.); fakeparton.e(0.); fakeparton.x(0.); fakeparton.y(0.); fakeparton.z(0.); fakeparton.x_t(0.);
					fakeparton.mass(xmq); fakeparton.e(xmq);
					//will reposition fake parton within the string in the next repair section...
					++str_len; ++numqrk; netqrk += (2*std::signbit(-fakeparton.id())-1);
				}
			}
			else{
				fakeparton.orig(-1); fakeparton.is_remnant(true); fakeparton.is_strendpt(true); fakeparton.string_id(prev_str); fakeparton.pos_str(-1);
				fakeparton.px(0.); fakeparton.py(0.); fakeparton.pz(0.); fakeparton.e(0.); fakeparton.x(0.); fakeparton.y(0.); fakeparton.z(0.); fakeparton.x_t(0.);
				fakeparton.mass(xmq); fakeparton.e(xmq);
				//will reposition fake parton within the string in the next repair section...
				++str_len; ++numqrk; netqrk += (2*std::signbit(-fakeparton.id())-1);
			}
		}
		
		//need to ensure that gluons do NOT form 'ends' of this string (move leading/trailing gluons to internal of the string)
		//also need to ensure that the 'fake' parton is stuck onto the appropriate end of the string (if added)
		//UNLESS it's a gluon loop - in that case, make sure that the color tags are set correctly in the next section
		//OR UNLESS it is completely fixed by reversing the order of the string (likely a junction type string with a 'trailing' segment)
		
		//moving on to color tag assignments
		//figure out how to 'properly' discern junction systems, and where the junction is in such a system...
			//for these, assign colors by "depth first search", need to denote a 'structure' for this system somehow...
			//for a first attempt, will treat a 'junction' type string by connecting a leading (or trailing) connection?
			//ex: q1 - g1 - g2 - q2 - g3 - g4 - g5 - g6 - q3, where g2 connects between g4 and g5
		
		//structure-wise, creating a 2d vector(n x n) to denote which partons in the string are next to other partons
		//allows for a unique definition of the string, and allows for the depth-search algorithm to assign color indices correctly
		std::vector<std::vector<bool>> connections; connections.resize(str_len, std::vector<bool>(str_len, false));
		
		//sorting the string into appropriate procedure based on number of quarks (& diquarks)
		//assuming that any 'simple' string (number of quarks <= 3) has some sort of 'physical' position ordering (with 'minor' deviations)
		//not bothering to check numq=1 as that is explicitly forbidden by the above check for color neutrality
		if(numqrk + 2*numdiqrk == 0){ //no quarks, this is a gluon loop
			if(numglu > 1){ //just chain the gluons into a loop
				for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;}
				connections[0][1]=true; connections[0][str_len-1]=true;
				connections[str_len-1][0]=true; connections[str_len-1][str_len-2]=true;
				
			}
			else{ //a single gluon can't be a color singlet, but we're going to 'force' this string to be one by splitting it into a q,qbar string
				//choosing a gluon mass - if implemented in the future, can (should) read this from the gluon entry itself (or set if necessary)
				//temporarily saving previously set mass here - here's a good place to check if this is even necessary?
				//maybe discard gluon if it is under some threshold of mass (eg < pion?)
				double temp_glumass = current_str[0].mass();
				current_str[0].mass( 2.*xmq + (gmax - 2.*xmq)*ran() );	// gluon virtuality
				
				//gluon decay function reads in the gluon (and the overwritten random mass), and writes the output q-qbar pair to qpair
				parton_collection qpair;
				gluon_decay(current_str[0], qpair);
				
				//swapping back original gluon mass
				//std::swap(current_str[0].mass, temp_glumass);
				current_str[0].mass(temp_glumass);
				
				//setting the vars of the q-qbar pair to the original gluon (and the flags of the original gluon)
				qpair[0].par(        current_str[0].par());        qpair[1].par(        current_str[0].par());
				qpair[0].is_shower(  current_str[0].is_shower());  qpair[1].is_shower(  current_str[0].is_shower());
				qpair[0].is_thermal( current_str[0].is_thermal()); qpair[1].is_thermal( current_str[0].is_thermal());
				qpair[0].orig(       current_str[0].orig());       qpair[1].orig(       current_str[0].orig());
				qpair[0].string_id(  current_str[0].string_id());  qpair[1].string_id(  current_str[0].string_id());
				qpair[0].pos_str(    0);                           qpair[1].pos_str(    0);
				qpair[0].is_remnant( true);                        qpair[1].is_remnant( true);
				if(qpair[0].par() >= 0){HH_shower[qpair[0].par()].is_used(true); HH_shower[qpair[0].par()].used_str(true); HH_shower[qpair[0].par()].is_decayedglu(true);}
				connections.resize(2, std::vector<bool>(2, false)); connections[0][1] = true; connections[1][0] = true;
				
				//since we work on 'current_str' for color assignments, we need to overwrite it with qpair...
				current_str.clear(); current_str.add(qpair);
				
			}
		}
		else if(numqrk + 2*numdiqrk == 2){ //two quarks, since color neutrality is already enforced, numqrk=2 and numdiqrk=0
			if((std::abs(current_str[0].id()) <= 6) && (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //GOOD string!  Don't need to re-sort.
				for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
				
			}
			else if((std::abs(current_str[0].id()) <= 6) != (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //There's a quark on one end of the string, but not the other...
				//going to make this easy, if the quark doesn't start the string (and so must end it), then we're going to reverse the order
				//this lets us treat both scenarios identically
				//UNLESS there's a fake parton, in which case we just stick that at the end after we reverse it so that the original quark is at the start
				if(std::abs(current_str[current_str.num()-1].id()) <= 6){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
				
				//now we find the quark that's not at the beginning at the string, and put it on the end by reversing the order from it to the end.
				//ex: q1-g1-g2-g3-q2-g4-g5  ==>  q1-g1-g2-g3-g5-g4-q2
				if(!fakepartonadded){for(int i=1;i<current_str.num();++i){if(std::abs(current_str[i].id())<=6){std::reverse(&current_str[i],(&current_str[current_str.num()-1])+1);break;}}}
				
				//and now the string is good!
				for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
				
			}
			else{ //Worst case scenario - the quarks are both somewhere internal to the string...
				//to remain consistent with the previous, we'll reverse the string from the start to the first quark, then from the second quark to the end
				//ex: g1-g2-q1-g3-q2-g4-g5  ==>  q1-g2-g1-g3-g5-g4-q2
				//UNLESS we have a fake parton, then we'll reverse from the quark to the closest end
					//then we'll reverse the whole string to put the quark in the first position, if necessary.
				//then we can just add the fake parton at the end.
				if(!fakepartonadded){
					for(int i=0; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[0], &current_str[i+1]); break;}}
					for(int i=1; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
				}
				else{ //fakepartonadded
					for(int i=0; i<current_str.num(); ++i){
						if(std::abs(current_str[i].id()) <= 6){
							if(2*i <= current_str.num()){std::reverse(&current_str[0], &current_str[i+1]);}
							else{std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1);}
							break;
						}
					}
					//now need to ensure that this string starts with the quark
					if(std::abs(current_str[current_str.num()-1].id()) <= 6){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
				}
				
				//and now the string is good!
				for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
				
			}
		}
		else if(numqrk + 2*numdiqrk == 3){ //can get this either by 3quarks (or antiquarks), or by diquark-quark string
			if(numqrk == 3){ //this is a 3 quark junction system
				//if two quarks are next to each other (first two?)(based on same position), then pretend that they're (+ the junction) are a diquark, and make a "diq" - - q string
				//if none of the quarks are next to each other, make sure that the ends are quarks (reversing from the 'outermost' quarks to the 'quarkless' ends)
					//then find the longest chain of gluons, split it in the middle for the junction, then attach the shortest gluon chain (only two chains!) to the junction
						//terminate this on the remaining quark.
				bool semi_diqrk = false; int firstq_diq = 0;
				for(int i=0; i<current_str.num()-1; ++i){if((std::abs(current_str[i].id()) <= 6) && (std::abs(current_str[i+1].id()) <= 6)){semi_diqrk = true; firstq_diq = i; break;}}
				
				if(semi_diqrk){ //we're treating this string configuration as a pseudo-diquark - quark string
					if((firstq_diq == 0) && (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //starts on the (p)diq, ends on the q; good!
						//connecting the two quarks of the p-diq together, and to the next link in the string
						connections[0][1] = true; connections[0][2] = true; connections[1][0] = true; connections[1][2] = true;
						//connecting the first link in the string to the first two quarks, and then to the next link
						connections[2][0] = true; connections[2][1] = true; connections[2][3] = true;
						for(int i=3; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[str_len-1][str_len-2]=true;
						
					}
					else if((std::abs(current_str[0].id()) <= 6) && (firstq_diq == current_str.num()-2)){ //starts on the q, ends on the (p)diq; reverse, then good!
						std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
						//connecting the two quarks of the p-diq together, and to the next link in the string
						connections[0][1] = true; connections[0][2] = true; connections[1][0] = true; connections[1][2] = true;
						//connecting the first link in the string to the first two quarks, and then to the next link
						connections[2][0] = true; connections[2][1] = true; connections[2][3] = true;
						for(int i=3; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[str_len-1][str_len-2]=true;
						
					}
					else if((firstq_diq == 0) || (firstq_diq == current_str.num()-2)){ //(p)diq terminates the string, but the other end isn't properly terminated; fix!
						if(firstq_diq == current_str.num()-2){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
						
						//now that the string starts with the (p)diq, make it so that the other end is terminated with the quark
						if(!fakepartonadded){
							for(int i=2; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
						}
						
						//and now the string is good!
						//connecting the two quarks of the p-diq together, and to the next link in the string
						connections[0][1] = true; connections[0][2] = true; connections[1][0] = true; connections[1][2] = true;
						//connecting the first link in the string to the first two quarks, and then to the next link
						connections[2][0] = true; connections[2][1] = true; connections[2][3] = true;
						for(int i=3; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[str_len-1][str_len-2]=true;
						
					}
					else if((std::abs(current_str[0].id()) <= 6) || (std::abs(current_str[current_str.num()-1].id()) <= 6)){//q terminates the string, but (p)diq doesn't terminate other end; fix!
						//first start by putting the terminating quark at the start of the string
						if(std::abs(current_str[current_str.num()-1].id()) <= 6){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
						
						//now we reverse from the FIRST of the p-diq to the end of the string, to terminate the string
						if(!fakepartonadded){
							for(int i=1; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
						}
						
						//now the string ends on the p-diq - reverse the string to put it at the front
						std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
						
						//and now the string is good!
						//connecting the two quarks of the p-diq together, and to the next link in the string
						connections[0][1] = true; connections[0][2] = true; connections[1][0] = true; connections[1][2] = true;
						//connecting the first link in the string to the first two quarks, and then to the next link
						connections[2][0] = true; connections[2][1] = true; connections[2][3] = true;
						for(int i=3; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[str_len-1][str_len-2]=true;
						
					}
					else{ //neither the (p)diq nor the q terminate the string; fix!
						//reversing from start of string to first q(or p-diq), then from second q (or 3rd if the first was a p-diq) to the end of the string
						if(!fakepartonadded){
							int k=0; bool diqstartsstr=false;
							for(int i=k; i<current_str.num(); ++i){
								if(std::abs(current_str[i].id()) <= 6){k=i; if(i == firstq_diq){k=i+1; diqstartsstr=true;} std::reverse(&current_str[0], (&current_str[k])+1); break;}
							}
							for(int i=k+1; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
							if(!diqstartsstr){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);} //reversing the string to put the p-diq at the start
						}
						else{ //fake parton has been added, the two actual quarks are must be next to each other
							for(int i=0; i<current_str.num(); ++i){
								if(std::abs(current_str[i].id()) <= 6){
									if(2*i <= current_str.num()){std::reverse(&current_str[0], (&current_str[i+1])+1);}
									else{std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
									break;
								}
							}
						}
						
						//and now the string is good!
						//connecting the two quarks of the p-diq together, and to the next link in the string
						connections[0][1] = true; connections[0][2] = true; connections[1][0] = true; connections[1][2] = true;
						//connecting the first link in the string to the first two quarks, and then to the next link
						connections[2][0] = true; connections[2][1] = true; connections[2][3] = true;
						for(int i=3; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[str_len-1][str_len-2]=true;
						
					}
				}
				else{ //this string configuration is a 'true' junction system (recombination did NOT add this junction!)
					//we want this string to look like q1-gN1-q2-gN2-q3
					//this string will then have a junction added in the middle of the longest gluon chain (either gN1 or gN2 will be broken in half)
					//the other gluon chain (gN2, or gN1) will be connected to the junction, and will terminate with the final quark
					
					if(fakepartonadded){
						//need to consider 3 cases - the 2 quarks can form strings as:  [1] g(n1)-q1-g(n2)-q2-g(n3); [2] q1-g(n1)-q2-g(n2); or [3] q1-g(n1)-q2
						//case [1] can be 'reduced' to the second case by reversing from whichever quark is closest to an end of the string, to that end
						//reverse [2(1)] if necessary to force the 'end' quark to occupy position 1 - now we can just tack the fake parton on the end
						//case [3] will need the fake parton added in the middle - will just attach it at the end (eg. g(n2) -> n2=0)
						if(     (std::abs(current_str[0].id()) <= 6) && (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //case [3]
							//not sure that we really need to do anything here to prep the string for color tagging (since we're just sticking the fake quark on the end)
							//the fake quark will then just get stuck internal to the gluon chain
						}
						else if((std::abs(current_str[0].id()) <= 6) || (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //case [2]
							if(!(std::abs(current_str[0].id()) <= 6)){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
						}
						else{ //case [1]
							int n,m; n=0; m=current_str.num()-1;
							for(int i=0;   i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){n=i; break;}}
							for(int i=n+1; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){m=i; break;}}
							if(n < current_str.num()-m){std::reverse(&current_str[0], (&current_str[n])+1);} //reverse to 1st quark (from start), else from 2nd quark (to end)
							else{std::reverse(&current_str[m], (&current_str[current_str.num()-1])+1); std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
						}
					}
					else{
						//also need to consider 3 cases here: can form strings as:  [1] q1-gN1-q2-gN2-q3; [2] q1-gN1-q2-gN2-q3-gN3; [3] gN1-q1-gN2-q2-gN3-q3-gN4
							//eg. both 'ends' properly terminate, only one 'end' terminates, or no end terminates
						//for case [2], reverse the end of the string, up to the last quark - now the string terminates on 'both' sides correctly, giving case [1]
						//for case [3], reverse from the start of the string to the first quark, and from the last quark to the end of the string - now have case [1]
						if((std::abs(current_str[0].id()) <= 6) && (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //case [1]
							//all is good, just set up the connections matrix
						}
						else if((std::abs(current_str[0].id()) <= 6) || (std::abs(current_str[current_str.num()-1].id()) <= 6)){ //case [2]
							//ensuring that the first parton is a quark
							if(std::abs(current_str[current_str.num()-1].id()) <= 6){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
							//finding last quark, reversing to the end to terminate the string
							for(int i=current_str.num()-1; i>0; --i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
						}
						else{ //case [3]
							//starting string with first quark
							for(int i=1; i<current_str.num(); ++i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[0], &current_str[i+1]); break;}}
							//ending string with last quark
							for(int i=current_str.num()-1; i>0; --i){if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}}
						}
					}
					
					//now the string looks like q1-gN1-q2-gN2-(q3) (q3 may be a fake parton)
					//count the gluons in N1 and in N2 - split the larger gluon chain in half, attach the shorter gluon chain to that break - now junction
					//finding the length of the 2 gluon chains
					int pos=0;
					int gluchainlen[2]; gluchainlen[0] = 0; gluchainlen[1] = 0;
					for(int i=pos+1; i<current_str.num(); ++i){if(current_str[i].id() == 21){++gluchainlen[0];}else{pos=i; break;}}
					for(int i=pos+1; i<current_str.num(); ++i){if(current_str[i].id() == 21){++gluchainlen[1];}else{       break;}}
					//breaking the longer chain in half, and attaching the shorter chain to it
					if(gluchainlen[0] >= gluchainlen[1]){ //breaking the first chain
						pos = gluchainlen[0]/2;
						int k = 0;
						connections[0][1] = true;
						for(int i=k+1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true; if(current_str[i+1].id() != 21){k=i+2; connections[i+1][i]=true; break;}}
						if(k<str_len-1){connections[k][k+1]=true;}
						connections[k][pos]=true; connections[k][pos+1]=true; connections[pos][k]=true; connections[pos+1][k]=true;
						for(int i=k+1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;}
						if(k != str_len-1){connections[str_len-1][str_len-2]=true;}
					}
					else{ //breaking the second chain
						pos = gluchainlen[1]/2;
						int k = str_len-1;
						connections[str_len-1][str_len-2] = true;
						for(int i=k-1; i>0; --i){connections[i][i-1]=true; connections[i][i+1]=true; if(current_str[i-1].id() != 21){k=i-2; connections[i-1][i]=true; break;}}
						connections[k][str_len-1-pos]=true; connections[k][str_len-1-pos-1]=true; connections[str_len-1-pos][k]=true; connections[str_len-1-pos-1][k]=true;
						if(k>0){connections[k][k-1]=true;}
						for(int i=k-1; i>0; --i){connections[i][i-1]=true; connections[i][i+1]=true;}
						if(k != 0){connections[0][1]=true;}
					}
					
				}
				
			}
			else{ //this is a q--diq system
				if((std::abs(current_str[0].id()) >= 1103) && (std::abs(current_str[0].id()) <= 5503) && ((current_str[0].id()/10)%10 == 0) && (std::abs(current_str[current_str.num()-1].id()) <= 6)){
					//GOOD string!  Don't need to re-sort.
					for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
				}
				else if((std::abs(current_str[0].id()) <= 6) &&
				  (std::abs(current_str[current_str.num()-1].id()) >= 1103) && (std::abs(current_str[current_str.num()-1].id()) <= 5503) && ((current_str[current_str.num()-1].id()/10)%10 == 0)){
					//Just going to reverse string to put the diquark first... now good
					std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
					for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
					
				}
				else if(((std::abs(current_str[0].id()) >= 1103) && (std::abs(current_str[0].id()) <= 5503) && ((current_str[0].id()/10)%10 == 0)) ||
						((std::abs(current_str[current_str.num()-1].id())>=1103)&&(std::abs(current_str[current_str.num()-1].id())<=5503)&&((current_str[current_str.num()-1].id()/10)%10==0))){
					//The diquark terminates this string, but the quark needs to be fixed
					//ensuring that the diquark starts the string
					if(!((std::abs(current_str[0].id()) >= 1103) && (std::abs(current_str[0].id()) <= 5503) && ((current_str[0].id()/10)%10 == 0))){
						std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
					}
					
					//find the quark that's not at the beginning at the string, and put it on the end by reversing the order from it to the end.
					//ex: q1-g1-g2-g3-q2-g4-g5  ==>  q1-g1-g2-g3-g5-g4-q2
					if(!fakepartonadded){for(int i=1; i<current_str.num(); ++i){
						if(std::abs(current_str[i].id()) <= 6){std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;}
					}}
					
					//and now the string is good!
					for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
					
				}
				else if((std::abs(current_str[0].id()) <= 6) || (std::abs(current_str[current_str.num()-1].id()) <= 6)){
					//The quark terminates one end of the string, but the diquark does not. (no fake parton here...)
					//going to start by reversing the string so that the quark starts it, if necessary
					if(std::abs(current_str[current_str.num()-1].id()) <= 6){std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);}
					
					//now we reverse the order from the diquark to the other end of the string, to properly terminate the string.
					if(!fakepartonadded){for(int i=1;i<current_str.num();++i){if((std::abs(current_str[i].id())>=1103) && (std::abs(current_str[i].id())<=5503) && ((current_str[i].id()/10)%10==0)){
						std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;
					}}}
					
					//now we reverse the entire string's order to put the diquark at the beginning.
					std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
					
					//and now the string is good!
					for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
					
				}
				else{ //Worst case scenario - the quarks are both somewhere internal to the string...
					//to remain consistent with the previous, we'll reverse the string from the start to the first (di)quark, then from the second (di)quark to the end
					//then we'll reverse the string if necessary to put the diquark at the beginning of the string
					//ex: g1-g2-q-g3-diq-g4-g5  ==>  q-g2-g1-g3-g5-g4-diq  ==>  diq-g4-g5-g3-g1-g2-q
					//UNLESS we have a fake parton, then we'll reverse from the diquark to the closest end
						//then we'll reverse the whole string to put the diquark in the first position, if necessary.
					//then we can just add the fake parton at the end.
					if(!fakepartonadded){
						for(int i=0; i<current_str.num(); ++i){
							if((std::abs(current_str[i].id()) <= 6) || ((std::abs(current_str[i].id()) >= 1103) && (std::abs(current_str[i].id()) <= 5503) && ((current_str[i].id()/10)%10 == 0))){
								std::reverse(&current_str[0], &current_str[i]); break;
						}}
						for(int i=1; i<current_str.num(); ++i){
							if((std::abs(current_str[i].id()) <= 6) || ((std::abs(current_str[i].id()) >= 1103) && (std::abs(current_str[i].id()) <= 5503) && ((current_str[i].id()/10)%10 == 0))){
								std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1); break;
						}}
						if(!((std::abs(current_str[0].id()) >= 1103) && (std::abs(current_str[0].id()) <= 5503) && ((current_str[0].id()/10)%10 == 0))){
							std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
						}
					}
					else{ //fakepartonadded
						for(int i=0; i<current_str.num(); ++i){
							if((std::abs(current_str[i].id()) >= 1103) && (std::abs(current_str[i].id()) <= 5503) && ((current_str[i].id()/10)%10 == 0)){
								if(2*i <= current_str.num()){std::reverse(&current_str[0], &current_str[i]);}
								else{std::reverse(&current_str[i], (&current_str[current_str.num()-1])+1);}
								break;
							}
						}
						//now need to ensure that this string starts with the quark
						if((std::abs(current_str[current_str.num()-1].id())>=1103) && (std::abs(current_str[current_str.num()-1].id())<=5503) && ((current_str[current_str.num()-1].id()/10)%10==0)){
							std::reverse(&current_str[0], (&current_str[current_str.num()-1])+1);
						}
					}
					
					//and now the string is good!
					for(int i=1; i<str_len-1; ++i){connections[i][i+1]=true; connections[i][i-1]=true;} connections[0][1]=true; connections[str_len-1][str_len-2]=true;
					
				}
			}
		}
		else if(numqrk+2*numdiqrk >= 4){
			//all of these are not 'simple' structures, and must contain at least one junction (except for a diq - - diq string, but that shouldn't really happen...
			
//**********************************************************************************************************************************************************************************
//	Beginning of string cutting into 'stringlets' string handling for >=4 quark strings
//**********************************************************************************************************************************************************************************
			if(cutstr){
				//maybe?, as a first 'attempt' that will conserve baryon number, collapse ALL quark systems with the same position (that were added to the endpoint as pairs) into diquarks
					//what to do if there are more than 2 quarks at the same position?
				//construct diquark/quark strings, then form a string from any leftovers (if any are left).  There MAY be a junction here, but only one
					//if there is more than one junction needed, cut it up into n single junction strings
				//q-qbar strings SHOULD NOT BE MADE until all of the 'diquark'-ish systems are cut
				
				//TODO: if this setup is left in place, then try to accomplish this without recursion, for ease of compiler optimization
				//the easiest way at the moment to accomplish the above is to split all the partons in the current strings into separate strings, each with its own string label
				//then call this function recursively (this will only have a depth of 1; we shouldn't ever call this function inside itself again) with all the current_str partons
				//all these partons will then be shovelled into the output remnant collection, and we'll need to continue to the next string...
				
				//need to break up current_str s.t. the different strings are given *new* string values!
				//to ensure that we have a unique value without knowing how many times we've recursively called the function, scan over remnants, find the max string_id, then add 1
				int new_string_id = 0;
				for(int i=0;i<SP_remnants.num();++i){if(SP_remnants[i].string_id() > new_string_id){new_string_id=SP_remnants[i].string_id();}}
				for(int i=0;i<SP_prepremn.num();++i){if(SP_prepremn[i].string_id() > new_string_id){new_string_id=SP_prepremn[i].string_id();}}
				++new_string_id;
				
				//now that we know where to start the new string_ids, we need to do the string cutting...
				//doing this differently than in the FORTRAN version
					//that version, in order, did: (1)pair off di(anti)quarks with (anti)quarks, (2)pair off q-qbar, (3)make q-q-q(bar) junction systems with the remainder
					//this one will, instead, form 'legs' with (anti)(di)quarks ending gluon chains, which will be paired off together based on distance/quarktype
					//TODO?: RNG to form q-q-q instead of q-qbar if multiple junctions (& antijunctions) present, to simulate junction-antijunction annihilaton
				if(fakepartonadded){current_str.add(fakeparton);}
				std::vector<std::vector<int>> ptn_legs;  std::vector<bool> valid_leg; std::vector<int> gluons; std::vector<bool> qtype_leg;
				std::vector<std::vector<int>> cutstrings;
				//if 2 valid legs closer than mindist_glu - remove BOTH of those legs from 'valid_leg' and find next closest gluon to any valid leg
				//if 0 valid legs are left, make a 'new' leg to form a gluon loop?
				for(int i=0; i<current_str.num(); ++i){
					if(current_str[i].id() == 21){gluons.push_back(i);}
					else{std::vector<int> tempptn; tempptn.push_back(i); ptn_legs.push_back(tempptn); tempptn.clear();}
				}
				for(int i=0; i<ptn_legs.size(); ++i){valid_leg.push_back(true);}
				for(int i=0; i<ptn_legs.size(); ++i){
					if((current_str[ptn_legs[i][0]].id() > 0 && current_str[ptn_legs[i][0]].id() <= 6) || (current_str[ptn_legs[i][0]].id() < -6)){qtype_leg.push_back(true);}
					else{qtype_leg.push_back(false);}
				}
				bool lastleg_gluloop = false;
				while(gluons.size() > 0){
					double mindist_glu  = 999999999999.; int glu_min=0 ; int leg_glu_min=-1;
					double mindist_legs = 999999999999.; int leg1_min=-1; int leg2_min=-1;
					for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){leg_glu_min = i; break;}}
					for(int i=0; i<ptn_legs.size(); ++i){
						if(!valid_leg[i]){continue;}
						for(int j=0; j<ptn_legs.size(); ++j){if(i==j){continue;} if(!valid_leg[j]){continue;}
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[j][ptn_legs[j].size()-1]]);
							if(dist_now < mindist_legs){mindist_legs = dist_now; leg1_min = i;   leg2_min = j;}
						}
						for(int j=0; j<gluons.size();   ++j){
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[gluons[j]]);
							if(dist_now < mindist_glu ){mindist_glu  = dist_now; leg_glu_min = i; glu_min = j;}
						}
					}
					int num_legs = 0; for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){++num_legs;}}
					if((num_legs > 0) && ((mindist_glu <= mindist_legs) || (num_legs < 2))){
						if(leg_glu_min>=0){ptn_legs[leg_glu_min].push_back(gluons[glu_min]); gluons.erase(gluons.begin()+glu_min);}
					}
					else if(num_legs > 0){
						if((leg1_min>=0)&&(leg2_min>=0)){valid_leg[leg1_min] = false; valid_leg[leg2_min] = false;}
					}
					else{
						std::vector<int> tempptn; tempptn.push_back(gluons[glu_min]); ptn_legs.push_back(tempptn); tempptn.clear();
						gluons.erase(gluons.begin()+glu_min); valid_leg.push_back(true); lastleg_gluloop = true;
					}
				}
				
				//now need to merge the legs
				//will use the distances between the leg ends to determine which legs to merge
				int num_legs = ptn_legs.size();
				for(int i=0; i<valid_leg.size(); ++i){valid_leg[i] = true;}
				//if we had a gluon loop for one of the legs, we write it out here
				if(lastleg_gluloop){
					std::vector<int> tempstr;
					//inserting odd entries at the beginning to make a 'zippered up' gluon loop
					for(int i=0; i<ptn_legs[ptn_legs.size()-1].size(); ++i){
						if(i%2==0){tempstr.push_back(ptn_legs[ptn_legs.size()-1][i]);}
						else{tempstr.insert(tempstr.begin(),ptn_legs[ptn_legs.size()-1][i]);}
					}
					cutstrings.push_back(tempstr);
					valid_leg[valid_leg.size()-1] = false; --num_legs;
				}
				//go over all (anti)(di)quark legs and merge them into q-*-q *-q or q-qbar strings based on momentum distance
				while(num_legs > 0){
					int leg1_sameq = -1; int leg2_sameq = -1; int leg3_sameq = -1; int leg1_oppq = -1; int leg2_oppq = -1;
					double min_dist12_sameq = -1.; double min_dist23_sameq = -1.; double min_dist13_sameq = -1.; double min_dist_oppq = -1.;
					bool checkJ = false; bool checkqq = false;
					
					for(int i=0; i<ptn_legs.size(); ++i){
						if(!valid_leg[i]){continue;}
						for(int j=0; j<ptn_legs.size(); ++j){if(i==j){continue;} if(!valid_leg[j]){continue;}
							if(qtype_leg[i] == qtype_leg[j]){
								for(int k=0; k<ptn_legs.size(); ++k){if(i==k || j==k){continue;} if(!valid_leg[k]){continue;}
									if(qtype_leg[i] == qtype_leg[k]){
										double dist1_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[j][ptn_legs[j].size()-1]]);
										double dist2_now = current_str[ptn_legs[j][ptn_legs[j].size()-1]].pDif2(current_str[ptn_legs[k][ptn_legs[k].size()-1]]);
										double dist3_now = current_str[ptn_legs[k][ptn_legs[k].size()-1]].pDif2(current_str[ptn_legs[i][ptn_legs[i].size()-1]]);
										double distsum_now = dist1_now + dist2_now + dist3_now;
										double distsum_prev = min_dist12_sameq + min_dist23_sameq + min_dist13_sameq;
										if((distsum_prev < 0) || (distsum_now < distsum_prev)){
											checkJ = true; leg1_sameq = i; leg2_sameq = j; leg3_sameq = k;
											min_dist12_sameq = dist1_now; min_dist23_sameq = dist2_now; min_dist13_sameq = dist3_now;
										}
									}
								}
							}
							else{
								double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[j][ptn_legs[j].size()-1]]);
								if((min_dist_oppq < 0) || (dist_now < min_dist_oppq)){checkqq = true; leg1_oppq = i; leg2_oppq = j; min_dist_oppq = dist_now;}
							}
						}
					}
					
					if(checkJ && ((min_dist_oppq < 0) || (min_dist12_sameq + min_dist23_sameq + min_dist13_sameq) < (3.*min_dist_oppq))){
						std::vector<int> tempstr;
						int sorted_legs[3]; sorted_legs[0]=leg1_sameq; sorted_legs[1]=leg2_sameq; sorted_legs[2]=leg3_sameq;
						if(ptn_legs[sorted_legs[1]].size() > ptn_legs[sorted_legs[0]].size()){std::swap(sorted_legs[0],sorted_legs[1]);}
						if(ptn_legs[sorted_legs[2]].size() > ptn_legs[sorted_legs[1]].size()){std::swap(sorted_legs[1],sorted_legs[2]);}
						if(ptn_legs[sorted_legs[1]].size() > ptn_legs[sorted_legs[0]].size()){std::swap(sorted_legs[0],sorted_legs[1]);}
						for(int i=0;i<ptn_legs[sorted_legs[0]].size();++i){tempstr.push_back(  ptn_legs[ sorted_legs[0] ][i]  );}
						for(int i=ptn_legs[sorted_legs[1]].size()-1;i>=0;--i){tempstr.push_back(ptn_legs[sorted_legs[1]][i]);}
						for(int i=ptn_legs[sorted_legs[2]].size()-1;i>=0;--i){tempstr.push_back(ptn_legs[sorted_legs[2]][i]);}
						cutstrings.push_back(tempstr);
						valid_leg[leg1_sameq] = false; valid_leg[leg2_sameq] = false; valid_leg[leg3_sameq] = false; num_legs -= 3;
					}
					else if(checkqq){
						std::vector<int> tempstr;
						for(int i=0;i<ptn_legs[leg1_oppq].size();++i){tempstr.push_back(ptn_legs[leg1_oppq][i]);}
						for(int i=ptn_legs[leg2_oppq].size()-1;i>=0;--i){tempstr.push_back(ptn_legs[leg2_oppq][i]);}
						cutstrings.push_back(tempstr);
						valid_leg[leg1_oppq] = false; valid_leg[leg2_oppq] = false; num_legs -= 2;
					}
				}
				
				//now cutstrings[i][j] has cut the 'current_str' into 'i' strings of 'j' partons
				//relabel all partons into corresponding string_id 'i' and pos_str 'j'
				for(int i=0;i<cutstrings.size();++i){
					for(int j=0;j<cutstrings[i].size();++j){current_str[cutstrings[i][j]].string_id(i+new_string_id); current_str[cutstrings[i][j]].pos_str(j);}
				}
				
				//recursive call, to deal with the current 'stringlets'.
				//if this is commented out, then we will skip this string wholesale
				parton_collection str_completed;
				stringprep(current_str, str_completed, false);
				
				//now need to reindex str_completed to prepare to add to final output
				//need to increment up the lastused_tag(color) so that this new string doesn't inherit anything from the previous string...
				++lastused_tag; int tagmax = 0;
				for(int i=0; i<str_completed.num(); ++i){
					if(str_completed[i].PY_par1() >= 0){str_completed[i].PY_par1( str_completed[i].PY_par1()+SP_prepremn.num() );}
					if(str_completed[i].PY_par2() >= 0){str_completed[i].PY_par2( str_completed[i].PY_par2()+SP_prepremn.num() );}
					if(str_completed[i].PY_dau1() >= 0){str_completed[i].PY_dau1( str_completed[i].PY_dau1()+SP_prepremn.num() );}
					if(str_completed[i].PY_dau2() >= 0){str_completed[i].PY_dau2( str_completed[i].PY_dau2()+SP_prepremn.num() );}
					if(str_completed[i].col()     >  0){str_completed[i].col(     str_completed[i].col()    +lastused_tag      );}
					if(str_completed[i].acol()    >  0){str_completed[i].acol(    str_completed[i].acol()   +lastused_tag      );}
					
					if(str_completed[i].col()  > tagmax){tagmax = str_completed[i].col();}
					if(str_completed[i].acol() > tagmax){tagmax = str_completed[i].acol();}
				}
				//set lastused_tag to be the last used color tag of these new strings...
				if(tagmax > 0){lastused_tag = tagmax;}
				
				//before we finally finish with this string, we need to fix the color tags...
				//these will start at '0' when they should have started at lastused_tag
				//easiest way to fix this will be to add lastused_tag to every color tag in str_completed
				
				//adding the current string to the output parton collection.
				if(str_completed.num() > 0){SP_prepremn.add(str_completed);}
				
				//now that this string is good, fill the next string
				start = par_i; prev_str = upd_str;
				current_str.clear();
				if(par_i<SP_remnants.num()){current_str.add(SP_remnants[par_i]);}
				continue;
			}
			
//**********************************************************************************************************************************************************************************
//	End of string cutting into 'stringlets' string handling for >=4 quark strings
//
//	!!!!!!!!!!!!ONLY ONE OF THESE SECTIONS SHOULD BE ACTIVE AT ANY GIVEN TIME!!!!!!!!!!!!
//	Just set the bool 'cutstr' to choose which...
//
//	Beginning of distance-measure parton rearrangement into junction leg string handling for >=4 quark strings (creating a single 'superstring')
//**********************************************************************************************************************************************************************************
			
			else{
				if(fakepartonadded){current_str.add(fakeparton);}
				std::vector<std::vector<int>> pseudojuns;
				std::vector<std::vector<int>> ptn_legs;  std::vector<bool> valid_leg; std::vector<int> gluons;
				std::vector<bool> qtype_leg;
				for(int i=0; i<current_str.num(); ++i){
					if(current_str[i].id() == 21){gluons.push_back(i);}
					else{std::vector<int> tempptn; tempptn.push_back(i); ptn_legs.push_back(tempptn); tempptn.clear();}
				}
				for(int i=0; i<ptn_legs.size(); ++i){valid_leg.push_back(true);}
				for(int i=0; i<ptn_legs.size(); ++i){
					if((current_str[ptn_legs[i][0]].id() > 0 && current_str[ptn_legs[i][0]].id() <= 6) || (current_str[ptn_legs[i][0]].id() < -6)){qtype_leg.push_back(true);}
					else{qtype_leg.push_back(false);}
				}
				
				while(gluons.size() > 0){
					double mindist_glu  = 9999999999999.; int glu_min=0 ; int leg_glu_min=-1;
					double mindist_legs = 999999999999.; int leg1_min; int leg2_min;
					for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){leg_glu_min = i; break;}}
					for(int i=0; i<ptn_legs.size(); ++i){
						if(!valid_leg[i]){continue;}
						for(int j=0; j<ptn_legs.size(); ++j){if(i==j){continue;} if(!valid_leg[j]){continue;} if(qtype_leg[i] != qtype_leg[j]){continue;}
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[j][ptn_legs[j].size()-1]]);
							if(dist_now < mindist_legs){mindist_legs = dist_now; leg1_min = i;   leg2_min = j;}
						}
						for(int j=0; j<gluons.size();   ++j){
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[gluons[j]]);
							if(dist_now < mindist_glu ){mindist_glu  = dist_now; leg_glu_min = i; glu_min = j;}
						}
					}
					int num_legs = 0; for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){++num_legs;}}
					if((mindist_glu <= mindist_legs) || (num_legs < 3)){ptn_legs[leg_glu_min].push_back(gluons[glu_min]); gluons.erase(gluons.begin()+glu_min);}
					else{
						mindist_glu  = 999999999999.;
						for(int j=0; j<gluons.size(); ++j){
							double dist_now = current_str[ptn_legs[leg1_min][ptn_legs[leg1_min].size()-1]].pDif2(current_str[gluons[j]]);
							if(dist_now < mindist_glu ){mindist_glu = dist_now; leg_glu_min = leg1_min; glu_min = j;}
								   dist_now = current_str[ptn_legs[leg2_min][ptn_legs[leg2_min].size()-1]].pDif2(current_str[gluons[j]]);
							if(dist_now < mindist_glu ){mindist_glu = dist_now; leg_glu_min = leg2_min; glu_min = j;}
						}
						std::vector<int> temp; temp.push_back(gluons[glu_min]); ptn_legs.push_back(temp); temp.clear();
						temp.push_back(leg1_min); temp.push_back(leg2_min); temp.push_back(ptn_legs.size()-1); pseudojuns.push_back(temp); temp.clear();
						valid_leg[leg1_min] = false; valid_leg[leg2_min] = false; valid_leg.push_back(true); gluons.erase(gluons.begin()+glu_min);
						qtype_leg.push_back(!qtype_leg[leg1_min]);
					}
				}
				
				//at this point, we have 2, 3, or more!
					//need to fix this so that the 'more' case is handled correctly...
					//do this by finding the two closest legs and merge them with a fake gluon.
					//repeat this procedure until there are 3 legs left (using legs, not fake gluon, to define distance). This forms the final junction
				//either way, write these remaining legs to psuedojunctions
				//TODO: a TOTAL distance minimization routine here
				int num_legs = 0; for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){++num_legs;}}
				int num_fakeglu = 0;
				double prev_leg = -1.;
				while(num_legs >= 3){
					if(prev_leg < 0){//this is the first merger
						int leg1_min = 0; int leg2_min = 0; double mindist = 999999999999.;
						//for(int i=0; i<ptn_legs.size()-num_fakeglu; ++i){for(int j=i+1; j<ptn_legs.size()-num_fakeglu; ++j){
						for(int i=0; i<ptn_legs.size(); ++i){for(int j=i+1; j<ptn_legs.size(); ++j){
							if(!valid_leg[i]){continue;} if(!valid_leg[j]){continue;} if(qtype_leg[i] != qtype_leg[j]){continue;}
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[j][ptn_legs[j].size()-1]]);
							if(dist_now < mindist){mindist = dist_now; leg1_min = i; leg2_min = j;}
						}}
						
						HHparton fakeglu;
						fakeglu.id(21); fakeglu.orig(-1); fakeglu.is_remnant(true); fakeglu.string_id(prev_str); fakeglu.pos_str(-1);
						fakeglu.px(0.); fakeglu.py(0.); fakeglu.pz(0.); fakeglu.e(0.); fakeglu.x(0.); fakeglu.y(0.); fakeglu.z(0.); fakeglu.x_t(0.);
						fakeglu.mass(0.); fakeglu.PY_stat(-99);
						current_str.add(fakeglu);
						
						std::vector<int> temp; temp.push_back(current_str.num()-1); ptn_legs.push_back(temp); temp.clear();
						temp.push_back(leg1_min); temp.push_back(leg2_min); temp.push_back(ptn_legs.size()-1); pseudojuns.push_back(temp); temp.clear();
						valid_leg[leg1_min] = false; valid_leg[leg2_min] = false; valid_leg.push_back(true); ++num_fakeglu; num_legs -= 2;
						qtype_leg.push_back(!qtype_leg[leg1_min]);
						
						//since for the first merger, we chose two legs, we need to find which leg we want to compare to, to find the next leg
						//this is done by going through *both* and finding which is closer to any other leg in the string
						//this will technically repeat for the next leg, but it's a small price to pay.
						mindist = 999999999999.;
						for(int i=0; i<ptn_legs.size()-num_fakeglu; ++i){
							if(!valid_leg[i]){continue;} if(qtype_leg[i] != qtype_leg[qtype_leg.size()-1]){continue;}
							double dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[leg1_min][ptn_legs[leg1_min].size()-1]]);
							if(dist_now < mindist){prev_leg = leg1_min;}
								   dist_now = current_str[ptn_legs[i][ptn_legs[i].size()-1]].pDif2(current_str[ptn_legs[leg2_min][ptn_legs[leg2_min].size()-1]]);
							if(dist_now < mindist){prev_leg = leg2_min;}
						}
					}
					else{//we have a prev parton to compare against
						int leg_min = 0; double mindist = 999999999999.;
						for(int i=0; i<ptn_legs.size()-num_fakeglu; ++i){
							if(!valid_leg[i]){continue;}  if(qtype_leg[i] == qtype_leg[prev_leg]){continue;}
							double dist_now = current_str[ptn_legs[prev_leg][ptn_legs[prev_leg].size()-1]].pDif2(current_str[ptn_legs[i][ptn_legs[i].size()-1]]);
							if(dist_now < mindist){mindist = dist_now; leg_min = i;}
						}
						
						HHparton fakeglu;
						fakeglu.id(21); fakeglu.orig(-1); fakeglu.is_remnant(true); fakeglu.string_id(prev_str); fakeglu.pos_str(-1);
						fakeglu.px(0.); fakeglu.py(0.); fakeglu.pz(0.); fakeglu.e(0.); fakeglu.x(0.); fakeglu.y(0.); fakeglu.z(0.); fakeglu.x_t(0.);
						fakeglu.mass(0.); fakeglu.PY_stat(-99);
						current_str.add(fakeglu);
						
						std::vector<int> temp; temp.push_back(current_str.num()-1); ptn_legs.push_back(temp); temp.clear();
						//int prev_glu = pseudojuns[pseudojuns.size()-1][pseudojuns[pseudojuns.size()-1].size()-1];
						int prev_glu = ptn_legs.size()-2;
						temp.push_back(prev_glu); temp.push_back(leg_min); temp.push_back(ptn_legs.size()-1); pseudojuns.push_back(temp); temp.clear();
						valid_leg[prev_glu] = false; valid_leg[leg_min] = false; valid_leg.push_back(true); ++num_fakeglu; --num_legs;
						qtype_leg.push_back(!qtype_leg[qtype_leg.size()-1]);
						
						prev_leg = leg_min;
					}
				}
				
				//now, we've added num_fakeglu fake gluons to current_str - we need to resize connections for the additional partons.
				connections.clear(); connections.resize(current_str.num(), std::vector<bool>(current_str.num(), false));
				
				//now there's either 2 or 3 legs left... this ends the string.
				std::vector<int> tempptn;
				for(int i=0; i<valid_leg.size(); ++i){if(valid_leg[i]){tempptn.push_back(i);}}
				pseudojuns.push_back(tempptn); tempptn.clear();
				
				//now, we have a nicely set up configuration, with partons in 'current_str' referenced to by ptn_legs, with 'junctions' determined by pseudojuns
				//to set up the connections 2d vector, we go through each of the ptn_legs entries, and connect each parton to the next (and to the one before)
				//then we set up the connections between legs via the pseudojuns array
					//for each 'junction', the last parton of the first n-1 legs connects to the first parton of the last leg
					//EXCEPT for the last 'junction', where the last elements of all connect to each other
				
				//starting by setting up connections within each leg, for all legs
					//unless the leg only has one parton(quark), in which case it gets handled only in the leg connection step after
				for(int i=0; i<ptn_legs.size(); ++i){
					if(ptn_legs[i].size() <= 1){continue;}
					connections[ptn_legs[i][0]][ptn_legs[i][1]]=true;
					for(int j=1; j<ptn_legs[i].size()-1; ++j){connections[ptn_legs[i][j]][ptn_legs[i][j+1]]=true; connections[ptn_legs[i][j]][ptn_legs[i][j-1]]=true;}
					connections[ptn_legs[i][ptn_legs[i].size()-1]][ptn_legs[i][ptn_legs[i].size()-2]]=true;
				}
				//now using pseudojuns to connect legs together
				for(int ipJn=0; ipJn<pseudojuns.size()-1; ++ipJn){
					std::vector<int> to_connect;
					for(int ipJL=0; ipJL<pseudojuns[ipJn].size()-1; ++ipJL){
						to_connect.push_back(ptn_legs[pseudojuns[ipJn][ipJL]][ptn_legs[pseudojuns[ipJn][ipJL]].size()-1]);
					}
					to_connect.push_back(ptn_legs[pseudojuns[ipJn][pseudojuns[ipJn].size()-1]][0]);
					
					for(int i=0; i<to_connect.size(); ++i){
						for(int j=0; j<to_connect.size(); ++j){
							if(i==j){continue;}
							connections[to_connect[i]][to_connect[j]]=true;
						}
					}
				}
				//connect the last legs in pseudojun[last]
				std::vector<int> to_connect;
				for(int ipJL=0; ipJL<pseudojuns[pseudojuns.size()-1].size(); ++ipJL){
					to_connect.push_back(ptn_legs[pseudojuns[pseudojuns.size()-1][ipJL]][ptn_legs[pseudojuns[pseudojuns.size()-1][ipJL]].size()-1]);
				}
				for(int i=0; i<to_connect.size(); ++i){
					for(int j=0; j<to_connect.size(); ++j){
						if(i==j){continue;}
						connections[to_connect[i]][to_connect[j]]=true;
					}
				}
				//and now we've set up the string & connections like we want!
				//we didn't rearrange the 'current_str' collection, but we shouldn't have needed to...
				//the subsequent code *should* be able to handle 'weird' cases (like starting from an intermediate gluon...)
				
	//******************************************************************************************************************************************************************************
	//	End of distance-measure parton rearrangement into junction leg string handling for >=4 quark strings
	//******************************************************************************************************************************************************************************
			}
		}
		
		//the strings are 'repaired' and the connection matrix is set up for the current string
		//now we're going to set up the color tags
		//first, we finally add the fake parton, if present, to the end of the current_str collection
			//note: connections should already include this; just hadn't placed it into current_str
		if(fakepartonadded && !(numqrk+2*numdiqrk >= 4)){current_str.add(fakeparton);}
		
		//setting up a 'labelled' list, used to determine which partons have been accounted for
		std::vector<bool> labelled; for(int i=0; i<current_str.num(); ++i){labelled.push_back(false);}
		//setting up a stack for depth-search algorithm for the string
		std::vector<int> stack;
		//setting up a list of junction partons - making a point to ensure that the color tags of these change!
		std::vector<std::vector<int>> junctions;
		//keeping track of what to set the next partons color and anticolor tags to (if it needs to be inherited or set to a new value)
		std::vector<std::vector<int>> settag; for(int i=0; i<current_str.num(); ++i){settag.push_back({0,0});}
		//need to increment up the lastused_tag so that this new string doesn't inherit anything from the previous string...
		++lastused_tag;
		//finding the first parton (would rather not start on a gluon, can be problematic, if the algorithm below isn't fixed to account)
		int first_ptn = 0;
		if((numqrk+2*numdiqrk > 0) && (current_str[0].id() == 21)){while(current_str[first_ptn].id() == 21){++first_ptn;}}
		//setting the flag for the first parton
		if((current_str[first_ptn].id() > 0 && current_str[first_ptn].id() <= 6) || (current_str[first_ptn].id() < -6)){settag[first_ptn] = {lastused_tag,0};}
		else if(current_str[first_ptn].id() == 21){                                                                     settag[first_ptn] = {lastused_tag,++lastused_tag};}
		else{                                                                                                           settag[first_ptn] = {0,lastused_tag};}
		
		//traversing the string, placing partons into the stack as encountered, then labelling the last parton entered (first-in, last-out).
		stack.push_back(first_ptn);
		while(stack.size() > 0){
			int current = stack.back(); stack.pop_back();
			if(!labelled[current]){
				
				//while we know that the connections vector stores the connections between ALL the partons, we want to know which partons THIS parton is connected to.
				//and also any partons that this parton is connected to that haven't been labelled yet...
				//we'll also load any connected partons that haven't been labelled into the stack
				std::vector<int> connected_all; std::vector<int> connected_new; int num_pre = 0;
				for(int i=0; i<connections.size(); ++i){if(connections[current][i]){
					connected_all.push_back(i); if(!labelled[i]){stack.push_back(i); connected_new.push_back(i);} if(i<current){++num_pre;}
				}}
				
				//now we'll set color/anticolor for quarks and figure out the tags for the next (unlabelled) parton(s)
				if((std::abs(current_str[current].id())<=6) || ((std::abs(current_str[current].id())>=1103) && (std::abs(current_str[current].id())<=5503) && ((current_str[current].id()/10)%10 == 0))){
					//setting tag
					if((current_str[current].id() > 0 && current_str[current].id() <= 6) || (current_str[current].id() < -6)){
						if(     settag[current][0] >   0){current_str[current].col( settag[current][0] );}
						else if(settag[current][0] == -2){current_str[current].col( ++lastused_tag );}
						else{
							//this junction leg needs to be 'repaired' by having the col and anticolor tags swapped for all partons back to the junction...
							current_str[current].col( (settag[current][1] > 0) ? settag[current][1] : ++lastused_tag );
							int ptn_now = current; std::vector<int> prev_ptns; prev_ptns.push_back(ptn_now);
							while(true){
								int conn_now = 0; int ptn_new = -1; bool foundnew = false;
								for(int i=0; i<connections.size(); ++i){
									if(connections[ptn_now][i]){
										bool ptn_isprev = false;
										for(int j=0; j<prev_ptns.size(); ++j){if(i == prev_ptns[j]){ptn_isprev = true;}}
										if(!ptn_isprev){++conn_now; foundnew = true; ptn_new = i; prev_ptns.push_back(i);}
									}
								}
								if((conn_now == 1) && foundnew){
									int temp = current_str[ptn_new].col();
									current_str[ptn_new].col(current_str[ptn_new].acol());
									current_str[ptn_new].acol(temp); ptn_now = ptn_new;
								}
								else{break;}
							}
						}
						//setting settag for next parton(s)
						if(connected_all.size() == 1){if(connected_new.size()){
							settag[connected_new[0]][1] = current_str[current].col();
							if(current_str[connected_new[0]].id() == 21){settag[connected_new[0]][0] = -1;}
						}}
						else{for(int i=0; i<connected_new.size(); ++i){  settag[connected_new[i]][0] = -2; settag[connected_new[i]][1] = -1;}}
						
					}
					else{
						if(     settag[current][1] >   0){current_str[current].acol( settag[current][1] );}
						else if(settag[current][1] == -2){current_str[current].acol( ++lastused_tag );}
						else{
							//this junction leg needs to be 'repaired' by having the col and anticolor tags swapped for all partons back to the junction...
							current_str[current].acol( (settag[current][0] > 0) ? settag[current][0] : ++lastused_tag );
							int ptn_now = current; std::vector<int> prev_ptns; prev_ptns.push_back(ptn_now);
							while(true){
								int conn_now = 0; int ptn_new = -1; bool foundnew = false;
								for(int i=0; i<connections.size(); ++i){
									if(connections[ptn_now][i]){
										bool ptn_isprev = false;
										for(int j=0; j<prev_ptns.size(); ++j){if(i == prev_ptns[j]){ptn_isprev = true;}}
										if(!ptn_isprev){++conn_now; foundnew = true; ptn_new = i; prev_ptns.push_back(i);}
									}
								}
								if((conn_now == 1) && foundnew){
									int temp = current_str[ptn_new].col();
									current_str[ptn_new].col(current_str[ptn_new].acol());
									current_str[ptn_new].acol(temp); ptn_now = ptn_new;
								}
								else{break;}
							}
						}
						//setting settag for next parton(s)
						if(connected_all.size() == 1){if(connected_new.size()){
							settag[connected_new[0]][0] = current_str[current].acol();
							if(current_str[connected_new[0]].id() == 21){settag[connected_new[0]][1] = -1;}
						}}
						else{for(int i=0; i<connected_new.size(); ++i){  settag[connected_new[i]][0] = -1; settag[connected_new[i]][1] = -2;}}
						
					}
					if(connected_all.size() > 1){
						std::vector<int> temp; temp.push_back(current); for(int i=0; i<connected_all.size(); ++i){temp.push_back(connected_all[i]);} junctions.push_back(temp);
					}
				}
				else{ //is gluon: set both tags and figure out tags for next (unlabelled) partons
				
					//setting color tag for gluons
					if(     settag[current][0] >   0){current_str[current].col( settag[current][0] );}
					else if(settag[current][0] == -1){current_str[current].col( ++lastused_tag );}
					else if(settag[current][0] == -2){current_str[current].col( ++lastused_tag );}
					//else if(!settag[0][current]){ERROR}
					
					//setting anticolor tag for gluons
					if(     settag[current][1] >   0){current_str[current].acol( settag[current][1] );}
					else if(settag[current][1] == -1){current_str[current].acol( ++lastused_tag );}
					else if(settag[current][1] == -2){current_str[current].acol( ++lastused_tag );}
					//else if(!settag[1][current]){ERROR}
					
					//figuring out tags for next (unlabelled) parton(s)
					//cases to catch:
						//0 - trivial case, current gluon is only connected to 1 unlabelled parton
						//1 - gluon loop case, 1st gluon of gluon loop (or other loop structure with cyclic color tags - gluon connected to 2 !unlabelled! gluons not in a junction)
							//warning, if the system *somehow* wound up with a nonjunction gluon connected to an unlabelled quark and another gluon, set tags 'appropriately'
								//though this should never actually happen...
						//2 - junction case, gluon connects to 'n' partons as part of a junction system - all partons must get new color tags
					if(     connected_all.size() == 1){
						//this is a special case - a gluon loop composed of only two gluons...
						//I'm going to assume 'hopefully' that there are only 2 gluons in the string (how else is this case going to happen?!) and set both tags in one go!
						current_str[0].col( settag[0][0] ); current_str[0].acol( settag[0][1] );
						current_str[1].col( settag[0][1] ); current_str[1].acol( settag[0][0] );
						break;
					}
					else if(connected_all.size() == 2){
						//check here for the gluon loop (>2 gluons) case - make sure not to overwrite any preexisting > 0 entries in settag!
						if((connected_new.size() == 2) && (current_str[connected_new[0]].id() == 21) && (current_str[connected_new[1]].id() == 21)){
							//one of these need to inherit the color, the other, anticolor...
							settag[connected_new[0]][0] = 0;                           settag[connected_new[0]][1] = current_str[current].col();
							settag[connected_new[1]][0] = current_str[current].acol(); settag[connected_new[1]][1] = -1;
						}
						//this case technically shouldn't happen... but we'll catch and treat, just in case...
						else if(connected_new.size() == 2){
							if((current_str[connected_new[0]].id() == 21) && (current_str[connected_new[1]].id() == 21)){
								settag[connected_new[0]][0] = current_str[current].acol(); settag[connected_new[0]][1] = -1;
								settag[connected_new[1]][0] = -1;                          settag[connected_new[1]][1] = current_str[current].col();
							}
							else if((std::abs(current_str[connected_new[0]].id()) <= 6) ||
							((std::abs(current_str[connected_new[0]].id())>=1103) && (std::abs(current_str[connected_new[0]].id())<=5503) && ((current_str[connected_new[0]].id()/10)%10==0))){
								settag[connected_new[0]][0] = current_str[current].acol(); settag[connected_new[0]][1] = -1;
								settag[connected_new[1]][0] = -1;                          settag[connected_new[1]][1] = current_str[current].col();
							}
							else{
								settag[connected_new[0]][0] = -1;                          settag[connected_new[0]][1] = current_str[current].col();
								settag[connected_new[1]][0] = current_str[current].acol(); settag[connected_new[1]][1] = -1;
							}
						}
						//this case is the 'expected' behavior with gluons internal to a string
						else if(connected_new.size() == 1){
							if(     (settag[current][0] >   0) && (settag[connected_new[0]][0] < 1)){settag[connected_new[0]][0] = current_str[current].acol();}
							else if((settag[current][0] == -1) && (settag[connected_new[0]][0] < 1)){settag[connected_new[0]][0] = -1;}
							else if((settag[current][0] == -2) && (settag[connected_new[0]][0] < 1)){settag[connected_new[0]][0] = current_str[current].acol();}
							                                                                    
							if(     (settag[current][1] >   0) && (settag[connected_new[0]][1] < 1)){settag[connected_new[0]][1] = current_str[current].col();}
							else if((settag[current][1] == -1) && (settag[connected_new[0]][1] < 1)){settag[connected_new[0]][1] = -1;}
							else if((settag[current][1] == -2) && (settag[connected_new[0]][1] < 1)){settag[connected_new[0]][1] = current_str[current].col();}
						}
						//else{}//shouldn't happen unless this is the final gluon of a gluon loop...
					}
					else{ //connected_all.size() >= 3
						//this is a gluon in a junction - and there might be already labelled partons connected...
						//already labelled ones do not need new tags
						//unlabelled partons that are also part of the junction need new tags - and inheritance set appropriately (likely opposite of this parton)
						//unlabelled partons (should only actually be 1 parton - not dealing with 'gluon loops' with junctions) that are NOT part of this junction need tags set
						//need to also catch configurations with a gluon being connected to multiple junctions (or adjacent to)
						
						//first, we construct a list of all the junctions that this parton is part of (it just might be multiple!)
							//this is done by assuming at least one junction exists
							//we then construct a vector of junction(s), with one junction initially, and the first entry is the current parton
							//then we run over the list of connected partons
							//for each of those, we will run over each list - if this parton is connected to all partons in the junction, it is taken as part of the junction
							//if the parton is not a part of any preexisting junctions, we create a new junction with the 'current' parton and this parton
							//once we've run through all connected partons, we look at the vector of junctions
							//any 'junction' composed of more than 2 partons is a 'true' junction and needs to be added to the master junction vector
							//any 'junction' composed of 2 partons is NOT a true junction, and the 2nd parton (not the 'current') is not part of any junction
						
						//setting up a few vectors
						std::vector<std::vector<int>> pseudoJunctions; std::vector<int> temppartons;
						temppartons.push_back(current); pseudoJunctions.push_back(temppartons); temppartons.clear();
						
						//this is running over ALL connected partons just in case we have a new parton that is part of a junction with an unlabelled one?
						for(int i=0; i<connected_all.size(); ++i){
							bool added = false;
							for(int j=0; j<pseudoJunctions.size(); ++j){
								bool conn_toall = true;
								for(int k=0; k<pseudoJunctions[j].size(); ++k){
									if(!connections[connected_all[i]][pseudoJunctions[j][k]]){conn_toall = false;}
								}
								if(conn_toall){pseudoJunctions[j].push_back(connected_all[i]); added = true;}
							}
							if(!added){
								temppartons.push_back(current); temppartons.push_back(connected_all[i]);
								pseudoJunctions.push_back(temppartons); temppartons.clear();
							}
						}
						
						//now we go ahead and take all the vectors in pseudoJunctions that are >2 partons ('true' junctions) and add them to junctions
						for(int i=0; i<pseudoJunctions.size(); ++i){if(pseudoJunctions[i].size() > 2){junctions.push_back(pseudoJunctions[i]);}}
						
						//from the above, we know which partons are part of which junctions, and which are not part of any junctions
						//to make what I'm doing next obvious, I'm going to split the connected partons into two vectors - ones part of a(n) junction(s), and ones not
						//any unlabelled partons in a junction get inheritance tags opposite of 'current'
						//any unlabelled partons not in a junction will continue the inheritance tags of 'current'
						//going to take care not to overwrite any preexisting ( != 0 ) settag entries
						
						//splitting all connected_new entries into "junction partons" and "not-junction partons"
						std::vector<int> ptns_inJun; std::vector<int> ptns_notJ;
						for(int i=0; i<connected_new.size(); ++i){
							bool found_parton = false;
							for(int j=0; j<pseudoJunctions.size(); ++j){
								for(int k=1; k<pseudoJunctions[j].size(); ++k){
									if(connected_new[i] == pseudoJunctions[j][k]){
										if(pseudoJunctions[j].size() > 2){ptns_inJun.push_back(connected_new[i]);}
										else{ptns_notJ.push_back(connected_new[i]);}
										found_parton = true; break;
									}
								}
								if(found_parton){break;}
							}
							//shouldn't ever happen, but we'll catch it anyways...
							//assuming that a not-found parton is not in a junction...
							if(!found_parton){ptns_notJ.push_back(connected_new[i]);}
						}
						
						//now, we run over both the inJun and notJ parton lists and set any unset settag entries
						//notJ partons get the inheritance tags of this parton (as if it is just continuing the string, which it is!)
						//inJ partons get flipped inheritance tags w.r.t 'current'
						for(int i=0; i<ptns_notJ.size();  ++i){
							if(     (settag[current][0] >   0) && (settag[ptns_notJ[i]][0] <  1)){settag[ptns_notJ[i]][0] = current_str[current].acol();}
							else if((settag[current][0] == -1) && (settag[ptns_notJ[i]][0] == 0)){settag[ptns_notJ[i]][0] = -1;}
							else if((settag[current][0] == -2) && (settag[ptns_notJ[i]][0] == 0)){settag[ptns_notJ[i]][0] = current_str[current].acol();}
							                                                                
							if(     (settag[current][1] >   0) && (settag[ptns_notJ[i]][1] <  1)){settag[ptns_notJ[i]][1] = current_str[current].col();}
							else if((settag[current][1] == -1) && (settag[ptns_notJ[i]][1] == 0)){settag[ptns_notJ[i]][1] = -1;}
							else if((settag[current][1] == -2) && (settag[ptns_notJ[i]][1] == 0)){settag[ptns_notJ[i]][1] = current_str[current].col();}
						}
						for(int i=0; i<ptns_inJun.size(); ++i){
							if(     (settag[current][0] >   0) && (settag[ptns_inJun[i]][0] <  1)){settag[ptns_inJun[i]][0] = -1;}
							else if((settag[current][0] == -1) && (settag[ptns_inJun[i]][0] == 0)){settag[ptns_inJun[i]][0] = -2;}
							else if((settag[current][0] == -2) && (settag[ptns_inJun[i]][0] == 0)){settag[ptns_inJun[i]][0] = -1;}
							                                                                 
							if(     (settag[current][1] >   0) && (settag[ptns_inJun[i]][1] <  1)){settag[ptns_inJun[i]][1] = -1;}
							else if((settag[current][1] == -1) && (settag[ptns_inJun[i]][1] == 0)){settag[ptns_inJun[i]][1] = -2;}
							else if((settag[current][1] == -2) && (settag[ptns_inJun[i]][1] == 0)){settag[ptns_inJun[i]][1] = -1;}
						}
					}
				}
				
				//now that we've labelled the color and prepared to set up any connected unlabelled partons, we're setting the fact that this parton is now labelled
				labelled[current] = true;
			}
		}
		
		//may need to clean up junctions vector from duplicate entries
		//any junction that is connected to another junction needs to have fake quarks added on the legs to use for pythia history
		//id's of the 3 quarks (inc. 'fake' ones if needed) that a junction should have will determine which particle serves as the mother for the junction
		//if the mother needs to have color, try to determine which pythia particle (BSM?) to set for the mother !!THIS SHOULDN'T HAPPEN CURRENTLY!!
		
		//finishing up settings for current string (namely, adding fake particles for pythia history) - need to trace legs of string/junction to assign the parents to endpoint quarks
		//first, sort all entries in the junctions vector, then remove duplicate entries (can arise from quarks being directly attached to junction)
		for(int i=0; i<junctions.size(); ++i){std::stable_sort(junctions[i].begin(),junctions[i].end());}
		std::stable_sort(junctions.begin(),junctions.end());
		std::vector<std::vector<int>>::iterator iter;
		iter = std::unique(junctions.begin(), junctions.end());
		junctions.resize(std::distance(junctions.begin(),iter));
		
		//now, for each junction, trace out all 3 legs to the endpoints
		//if it hits another junction before that, we'll need to add 2 fake partons (one for each junction) - these will need to be color correlated to the internal gluon link
			//these fake partons do NOT need to have daughters set
		//then the three partons will need to have their parent id's set to a fake parent added (colorless, unless we have a 'colored' junction?)
			//these fake particles WILL need to have daughters set?
		std::vector<std::vector<int>> py_siblings;
		for(int iJun=0; iJun<junctions.size(); ++iJun){
			//for this junction, take each parton connected - run along chain back to endpoint or next junction;
			//done by finding a parton that this parton is connected to, that is NOT a parton in the junction, UNLESS this parton is a quark (means this parton is the end of the leg)
			//move along this gluon chain until a quark is reached, OR a gluon that connects to more than 2 gluons (a junction)
			//if we have an internal junction link, just leave the original gluon in this list (we'll use it later to assign the color tag of the fake parton correctly)
			std::vector<int> temp;
			for(int iPar=0; iPar<junctions[iJun].size(); ++iPar){
				if(current_str[junctions[iJun][iPar]].id() != 21){temp.push_back(junctions[iJun][iPar]);}
				else{
					//find the next parton in this leg
					//start by looking for all connected partons
					std::vector<int> connected;
					for(int i=0; i<connections.size(); ++i){if(connections[junctions[iJun][iPar]][i]){connected.push_back(i);}}
					
					//now we go over the list of connected partons, and find the one that's not in the junction
					//then we clear the 'connected' list and add the original junction parton; then we add this parton we just found (that's not connected to the junction)
					//if we don't find a parton, then this gluon is the 'whole' chain (it connects to 2 junctions itself)
					bool findnext = false;
					for(int i=0; i<connected.size(); ++i){
						int nconn=0; for(int j=0; j<connections.size(); ++j){if(connections[connected[i]][j]){++nconn;}}
						if((nconn <=2 && current_str[connected[i]].id() == 21) || (nconn == 1)){
							int tempconn = connected[i]; connected.clear(); connected.push_back(junctions[iJun][iPar]); connected.push_back(tempconn);
							if(current_str[tempconn].id() == 21){findnext=true;}
							break;
						} //the below condition only triggers IF we don't find a valid candidate
						else if(i == connected.size()-1){connected.clear(); connected.push_back(junctions[iJun][iPar]);}
					}
					
					//now we move along the chain as long as we have a gluon that's not connected to a junction
					//if we reach a quark that's only connected to 1 parton, we've found the endpoint of the leg
					//if we reach a gluon that's connected to >= 3 partons, or a quark that's connected to >= 2 partons, we've hit the next junction
						//in this case, we just leave the original gluon
					if(findnext){
						int nconn=0; for(int i=0; i<connections.size(); ++i){if(connections[connected[1]][i]){++nconn;}}
						while((current_str[connected[connected.size()-1]].id() == 21) && (nconn <= 2)){
							for(int i=0; i<connections.size(); ++i){
								if(connections[connected[connected.size()-1]][i]){
									bool skip = false;
									for(int j=0; j<connected.size(); ++j){if(i == connected[j]){skip = true;}}
									if(!skip){connected.push_back(i); break;}
								}
							}
							//now need to find nconn
							nconn=0; for(int i=0; i<connections.size(); ++i){if(connections[connected[connected.size()-1]][i]){++nconn;}}
						}
						
						//we've put the last 'parton' into connected - if it's a quark that's only connected to a single parton(gluon), it's the endpoint
						//else we've either got a gluon connected to >2 partons, or a quark connected to >=2 partons, so we've hit a junction
						//if(current_str[connected[connected.size()-1]].id() != 21 && nconn < 2){temp.push_back(connected[connected.size()-1]);}
						if(current_str[connected[connected.size()-1]].id() != 21){temp.push_back(connected[connected.size()-1]);}
						else{temp.push_back(connected[0]);}
					}
					else{
						if(current_str[connected[connected.size()-1]].id() != 21){temp.push_back(connected[connected.size()-1]);}
						else{temp.push_back(connected[0]);}
					}
				}
			}
			py_siblings.push_back(temp); temp.clear();
		}
		
		//now, for each sibling 'group' in siblings, need to add a fake mother and set the daughters field of it appropriately
		//need to break sorting of the string to put siblings together (set daughters of the mother as [first,last] - if 3,4,5 then set 3,5 and 4 MUST be the remaining sibling)
		//and make sure to set the mother of the sibling partons too!
		//also, if this parton is a gluon, then it means that we need to add in a fake q,qbar pair on the opposing ends of the internal gluon chain
		//nicely though, we just need to add in a fake q (or qbar); use inheritflags vector to pick and tag for this gluon appropriately
			//HOWEVER, the internal gluon chain will all have to have mothers set to BOTH of these fake partons
		//and if there are more than 3 partons, this must have been a multiple junction system that 'collapsed'.
			//treating this by 'throwing a Hail Mary' and attaching them ALL to a single mother; hope that PYTHIA can figure out whats going on... NOPE DOESN'T WORK.
		//declaring a collection to be added to final output
		parton_collection str_completed;
		std::vector<int> fakequarktype; for(int i=0; i<current_str.num(); ++i){fakequarktype.push_back(0);}
		for(int iSib=0; iSib<py_siblings.size(); ++iSib){
			int nq=0; int nqbar=0;
			std::vector<std::vector<int>> gluchains;
			std::vector<int> fakeptns; std::vector<int> fakeptns_colid;
			//first pass, run through and catch any gluons (from internal gluon chains that connect 2 junctions)
			//toss those into a gluon chain, put that chain into a vector of multiple chains (for just this junction)
			//toss in the fake parton for this side of each of those chains
			for(int iPtn=0; iPtn<py_siblings[iSib].size(); ++iPtn){
				if(current_str[py_siblings[iSib][iPtn]].id() == 21){
					//this is a gluon that's part of a chain connecting two junctions
					//find the next gluon in the chain (if present!)
					//we're going to make sure to start the chain with this gluon - if it doesn't start with this gluon, it WILL break.
					std::vector<int> chain; chain.push_back(py_siblings[iSib][iPtn]);
					bool chain_gt_1 = true;
					for(int i=0; i<connections.size(); ++i){if((connections[py_siblings[iSib][iPtn]][i]) && (current_str[i].id() == 21)){chain.push_back(i);}}
					
					for(int i=0; i<chain.size(); ++i){
						bool issib = false;
						for(int j=0;j<py_siblings[iSib].size();++j){if(chain[i] == py_siblings[iSib][j]){issib = true; break;}}
						//extra check to catch if this gluon is a sibling (case where the gluon is part of this junction, but is on a terminated leg)
							//can do this since junctions and py_siblings are 1-1 correlated. (eg, i'th junction is the i'th py_siblings)
						for(int j=0;  j<junctions[iSib].size();++j){if(chain[i] ==   junctions[iSib][j]){issib = true; break;}}
						//need yet another check - run over all the junctions and make sure that chain[i] is not in ANY junction that sib[iSib][iPtn] is also in
						for(int j=0;j<junctions.size();++j){
							bool ptn1_inc = false; bool ptn2_inc = false;
							for(int k=0;k<junctions[j].size();++k){
								if(py_siblings[iSib][iPtn] == junctions[j][k]){ptn1_inc = true;}
								else if(chain[i]           == junctions[j][k]){ptn2_inc = true;}
							}
							if(ptn1_inc && ptn2_inc){issib = true; break;}
						}
						
						//if this is not a sibling parton, then it's the next parton in the chain... right?
						if(!issib){int temp = chain[i]; chain.clear(); chain.push_back(py_siblings[iSib][iPtn]); chain.push_back(temp); break;}
						else if(i==chain.size()-1){chain.clear(); chain.push_back(py_siblings[iSib][iPtn]); chain_gt_1 = false;}
					}
					
					
					//now we move along the chain until we hit the gluon attached to the next junction
					if(chain_gt_1){
						int nconn=0; for(int i=0; i<connections.size(); ++i){if(connections[chain[1]][i]){++nconn;}}
						while(current_str[chain[chain.size()-1]].id() == 21 && nconn <= 2){
							for(int i=0; i<connections.size(); ++i){
								if(connections[chain[chain.size()-1]][i]){
									bool skip = false;
									for(int j=0; j<chain.size(); ++j){if(i == chain[j]){skip = true;}}
									if(!skip){chain.push_back(i); break;}
								}
							}
							//now need to find nconn
							nconn=0; for(int i=0; i<connections.size(); ++i){if(connections[chain[chain.size()-1]][i]){++nconn;}}
						}
					}
					
					//now adding this chain to the vector of all chains for this junction
					gluchains.push_back(chain);
					
					//need to create a fake (anti)quark, set it as the mother of the gluon chain, then set the mother of it to the junction mother
					//previous attempt to use settag failed...
					//instead, we'll first look at the fakequarktype vector; if it's nonzero, then this fake parton is predetermined
						//since the other end of the gluon chain was chosen...
					//if it's zero, then we'll want to look at the other partons in this junction
						//but first, check the other end of the chain, just in case... if that was already set, then set this one.
						//if there are two other (anti)quarks, then choose an (anti)quark
						//if there are two other quarks, and they're different (one antiquark, one quark), then choice is arbitrary.
						//if there is one other quark, then just choose the same type as it.
						//if there are no quarks, then follow one of the chains to the next junction, follow from that junction to the next quark, choose the opposite
							//if we hit yet another junction, **** it, just arbitrarily choose
						//make sure to set the other end of the gluon chain 'fakequarktype' to be the opposite of this one
					int fakeid = fakequarktype[chain[0]]; if(fakeid == 0){fakeid = -fakequarktype[chain[chain.size()-1]];}
					
					if(fakeid == 0){
						int signJ = 0;
						for(int i=0; i<py_siblings[iSib].size(); ++i){
							if(     (current_str[py_siblings[iSib][i]].id() == 21) || (current_str[py_siblings[iSib][i]].id() == 0)){continue;}
							else if((current_str[py_siblings[iSib][i]].id() > 0 && current_str[py_siblings[iSib][i]].id() <=  6) || (current_str[py_siblings[iSib][i]].id() < -6)){++signJ;}
							else if((current_str[py_siblings[iSib][i]].id() < 0 && current_str[py_siblings[iSib][i]].id() >= -6) || (current_str[py_siblings[iSib][i]].id() >  6)){--signJ;}
						}
						
						if(     signJ >= 0){fakeid =  1; fakequarktype[chain[chain.size()-1]] = -1;}
						else if(signJ <  0){fakeid = -1; fakequarktype[chain[chain.size()-1]] =  1;}
					}
					
					if(fakeid == 1){++nq;}else{++nqbar;}
					int col_id = (fakeid == 1) ? current_str[chain[0]].acol() : current_str[chain[0]].col();
					fakeptns.push_back(fakeid); fakeptns_colid.push_back(col_id);

				}
				else{
					//still need to increment nq, nqbar counters...
					//going to just treat diquarks as antiquarks (and antidiquarks as quarks)
					if((current_str[py_siblings[iSib][iPtn]].id() > 0 && current_str[py_siblings[iSib][iPtn]].id() <= 6) || (current_str[py_siblings[iSib][iPtn]].id() < -6)){++nq;}
					else{++nqbar;}
				}
			}
			
			//setting the mother to the correct id (if baryon or antibaryon)(also, if the junction is NOT a color singlet, then the mother shouldn't be one as well)
			int m_id = 0; int ifake = 0; int m_qrks[3]={0,0,0}; bool choosechg = false; int chgjun = 0;
			
			//will fail if this junction has !=3 legs (which it shouldn't be!)
			if(py_siblings[iSib].size() > 3){choosechg = true;}
			for(int iPtn=0; iPtn<py_siblings[iSib].size(); ++iPtn){
				if(current_str[py_siblings[iSib][iPtn]].id() == 21){
					m_qrks[iPtn] = fakeptns[ifake]; ++ifake;
					int chg = (m_qrks[iPtn] % 2 == 0) ? 2 : -1;
					chgjun += chg*(2*std::signbit(-m_qrks[iPtn])-1);
				}
				else if((std::abs(current_str[py_siblings[iSib][iPtn]].id()) > 0) && (std::abs(current_str[py_siblings[iSib][iPtn]].id()) <= 6)){
					m_qrks[iPtn] = current_str[py_siblings[iSib][iPtn]].id();
					int chg = (current_str[py_siblings[iSib][iPtn]].id() % 2 == 0) ? 2 : -1;
					chgjun += chg*(2*std::signbit(-current_str[py_siblings[iSib][iPtn]].id())-1);
				}
				else if(current_str[py_siblings[iSib][iPtn]].id() < -6){
					choosechg = true;
					int id1 = (current_str[py_siblings[iSib][iPtn]].id()/100 )%10; chgjun += (id1%2 == 0) ? -2 : 1;
					int id2 = (current_str[py_siblings[iSib][iPtn]].id()/1000)%10; chgjun += (id2%2 == 0) ? -2 : 1;
				}
				else{
					choosechg = true;
					int id1 = (current_str[py_siblings[iSib][iPtn]].id()/100 )%10; chgjun += (id1 % 2 == 0) ? 2 : -1;
					int id2 = (current_str[py_siblings[iSib][iPtn]].id()/1000)%10; chgjun += (id2 % 2 == 0) ? 2 : -1;
				}
			}
			
			if(!choosechg){
				if(std::abs(m_qrks[0]) < std::abs(m_qrks[1])){std::swap(m_qrks[0],m_qrks[1]);}
				if(std::abs(m_qrks[1]) < std::abs(m_qrks[2])){std::swap(m_qrks[1],m_qrks[2]);}
				if(std::abs(m_qrks[0]) < std::abs(m_qrks[1])){std::swap(m_qrks[0],m_qrks[1]);}
				
				if((m_qrks[0]*m_qrks[1]<=0) || (m_qrks[0]*m_qrks[2]<=0)/* || (m_qrks[1]*m_qrks[2]<0)*/ ){
					int agree[2]={0,0};
					if(m_qrks[0]*m_qrks[1]>0){agree[0]=m_qrks[0]; agree[1]=m_qrks[1];}
					else if(m_qrks[0]*m_qrks[2]>0){agree[0]=m_qrks[0]; agree[1]=m_qrks[2];}
					else{agree[0]=m_qrks[1]; agree[1]=m_qrks[2];}
					
					m_id = -1*(1000*agree[0] + 100*agree[1] + 3);
				}
				else{m_id = 1000*m_qrks[0] + 100*m_qrks[1] + 10*m_qrks[2] + 4; /*if(m_qrks[0]==m_qrks[2]){m_id += 2;}*/}
			}
			else{
				if(     chgjun ==  6){m_id =  2224;} // chgjun/3 ==  2
				else if(chgjun ==  3){m_id =  2214;} // chgjun/3 ==  1
				else if(chgjun ==  0){m_id =  2114;} // chgjun/3 ==  0
				else if(chgjun == -3){m_id = -2214;} // chgjun/3 == -1
				else if(chgjun == -6){m_id = -2224;} // chgjun/3 == -2
				else{JSWARN << "Junction with noninteger or >2 elecric charge present (" << chgjun << "/3). Unsure of how to properly treat ... " << chgjun; m_id = 5554;}
			}
			
			//now we add this junction's mother to the str_completed collection
			//will set its daughter entries to the first and last entries in the junction sibling list AFTER completing this junction
			HHparton mother; mother.id(m_id); mother.PY_stat(-11);/*mother.PY_stat=-12;*/ mother.string_id(prev_str); mother.orig(-1); mother.is_remnant(true);
			mother.mass(1.); mother.e(1.);
			if(     m_id ==  1103){mother.acol( ++lastused_tag ); mother.mass(2.*xmq); mother.e(2.*xmq);}
			else if(m_id == -1103){mother.col(  ++lastused_tag ); mother.mass(2.*xmq); mother.e(2.*xmq);}
			else if(m_id ==  1   ){mother.col(  ++lastused_tag ); mother.mass(   xmq); mother.e(   xmq);}
			else if(m_id == -1   ){mother.acol( ++lastused_tag ); mother.mass(   xmq); mother.e(   xmq);}
			str_completed.add(mother);
			int mother_pos = str_completed.num()-1;
			
			//now we iterate over the junction sibling list
				//for a quark listed, we set its PY_par1 entry to the junction mother
				//for a gluon listed:
					// 1) we add the corresponding fake parton to the str_completed collection
					// 2) the fake parton gets par1 entry set to junction mother, dau1 entry set to start of the corresponding gluon chain, and dau2 set to end of the gluon chain
					// 3) the corresponding gluon chain has one of its mother entries set to this fake parton (written to par1 if it hasn't been set yet, otherwise to par2)
			ifake = 0;
			for(int iPtn=0; iPtn<py_siblings[iSib].size(); ++iPtn){
				if(current_str[py_siblings[iSib][iPtn]].id() != 21){
					//need to set the mother entries not just for this quark, but also for any gluons that trail back to the junction (if any)
					//finding the partons(gluons) in the leg, trailing back to the junction
					std::vector<int> jun_leg; jun_leg.push_back(py_siblings[iSib][iPtn]);
					int new_conn=0; for(int i=0; i<connections.size(); ++i){if(connections[jun_leg.back()][i]){++new_conn;}}
					while(new_conn < 2){
						std::vector<int> next;
						for(int i=0; i<connections.size(); ++i){
							if(connections[jun_leg.back()][i]){
								bool skip = false;
								for(int j=0; j<jun_leg.size(); ++j){if(i == jun_leg[j]){skip = true;}}
								if(!skip){next.push_back(i);}
							}
						}
						if(next.size() == 1){jun_leg.push_back(next[0]); new_conn=1;}
						else{new_conn = next.size(); break;}
					}
					
					//setting the mother entry for all the partons in this leg.
					for(int i=0; i<jun_leg.size(); ++i){current_str[jun_leg[i]].PY_par1(mother_pos);}
					//current_str[py_siblings[iSib][iPtn]].PY_par1(mother_pos);
				}
				else{
					HHparton fake; fake.id(fakeptns[ifake]); fake.PY_stat(-21); fake.string_id(prev_str); fake.orig(-1); fake.is_remnant(true);
					fake.PY_par1(mother_pos); fake.PY_dau1(-2-gluchains[ifake][0]); fake.PY_dau2(-2-gluchains[ifake][gluchains[ifake].size()-1]);
					(fake.id() == 1) ? fake.col(fakeptns_colid[ifake]) : fake.acol(fakeptns_colid[ifake]);
					fake.mass(xmq); fake.e(xmq);
					str_completed.add(fake);
					for(int ichain=0;ichain<gluchains[ifake].size();++ichain){
						if(current_str[gluchains[ifake][ichain]].PY_par1() == -1){current_str[gluchains[ifake][ichain]].PY_par1( str_completed.num()-1 );}
						else{                                                     current_str[gluchains[ifake][ichain]].PY_par2( str_completed.num()-1 );}
					}
					++ifake;
				}
			}
			
			//need to do this down here; the last added parton in fake partons will correspond to the last added fake parton in str_completed
			if(current_str[py_siblings[iSib][0]].id() != 21){
				 str_completed[mother_pos].PY_dau1(-2-py_siblings[iSib][0]);
			}
			else{str_completed[mother_pos].PY_dau1(mother_pos+1);}
			if(current_str[py_siblings[iSib][py_siblings[iSib].size()-1]].id() != 21){
				 str_completed[mother_pos].PY_dau2(-2-py_siblings[iSib][py_siblings[iSib].size()-1]);
			}
			else{str_completed[mother_pos].PY_dau2(str_completed.num()-1);}
		}
		
		//reindex appropriate str_completed daughter entries, then add the current_str partons to str_completed
		for(int i=0; i<str_completed.num(); ++i){
			if(str_completed[i].PY_dau1() < -1){str_completed[i].PY_dau1( -(2+str_completed[i].PY_dau1())+str_completed.num() );}
			if(str_completed[i].PY_dau2() < -1){str_completed[i].PY_dau2( -(2+str_completed[i].PY_dau2())+str_completed.num() );}
		}
		str_completed.add(current_str);
		
		//going to remove any fake gluons that had to be added for multijunction strings - PYTHIA doesn't handle them well
		int istrfix = 0;
		while(istrfix < str_completed.num()){
			if(str_completed[istrfix].PY_stat() == -99){ //this is a fake gluon that had to be added - need to remove it and repair any appropriate py_*** entries for other partons
				//since this parton connects two fake partons, we need to repair the color entries for those partons
				str_completed[str_completed[istrfix].PY_par2()].col(  str_completed[str_completed[istrfix].PY_par1()].acol() );
				str_completed[str_completed[istrfix].PY_par2()].acol( str_completed[str_completed[istrfix].PY_par1()].col()  );
				
				//repairing the py_*** entries for all partons that have values >= istrfix
				for(int i=0; i<str_completed.num(); ++i){
					if(     str_completed[i].PY_par1() >  istrfix){str_completed[i].PY_par1(str_completed[i].PY_par1()-1);}
					else if(str_completed[i].PY_par1() == istrfix){str_completed[i].PY_par1(-1);}
					if(     str_completed[i].PY_par2() >  istrfix){str_completed[i].PY_par2(str_completed[i].PY_par2()-1);}
					else if(str_completed[i].PY_par2() == istrfix){str_completed[i].PY_par2(-1);}
					if(     str_completed[i].PY_dau1() >  istrfix){str_completed[i].PY_dau1(str_completed[i].PY_dau1()-1);}
					else if(str_completed[i].PY_dau1() == istrfix){str_completed[i].PY_dau1(-1);}
					if(     str_completed[i].PY_dau2() >  istrfix){str_completed[i].PY_dau2(str_completed[i].PY_dau2()-1);}
					else if(str_completed[i].PY_dau2() == istrfix){str_completed[i].PY_dau2(-1);}
				}
				
				//removing this fake gluon, and setting up to check the next parton (which is now in this position)
				str_completed.remove(istrfix);
				--istrfix;
			}
			++istrfix;
		}
		
		std::vector<std::vector<int>> py_rearrange;
		for(int i=0; i<str_completed.num(); ++i){
			std::vector<int> py_tmp; py_tmp.push_back(i);
			py_tmp.push_back(str_completed[i].PY_par1());
			py_tmp.push_back(str_completed[i].PY_dau1());
			py_tmp.push_back(-1);
			py_tmp.push_back(-1);
			py_tmp.push_back(-1);
			py_tmp.push_back(0);
			py_rearrange.push_back(py_tmp);
		}
		//use > to sort in descending order, < to sort in ascending
		std::stable_sort(py_rearrange.begin(), py_rearrange.end(), [](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[1]<vec2[1];});
		
		//now there are 'clusters' of daughters of mothers
			//but, may need to rearrange clusters s.t. clusters that contain one of a 'dual mother' for another cluster are next to each other
			//also, may need to rearrange the partons in those particular clusters so that the 'dual mothers' are next to each other
		//handle this by sorting the clusters as a whole, then we'll sort inside of clusters
		
		//first, we're going to set the 4th through 7th fields for each of the py_re entries
		//this will allow us to denote which daughter(s) that partons in the cluster have
		//we're going to assume a maximum of 2 daughters - if there are more, then we'll need to create another 'fake' generation so that PYTHIA can handle it
		int srt_pos = -1; int end_pos = -1; int prev_id = -1; int ncluster = -1;
		for(int i=0; i<py_rearrange.size(); ++i){
			if(py_rearrange[i][1] != prev_id){
				if(srt_pos == -1){srt_pos = i; prev_id = py_rearrange[i][1];}
				else{
					++ncluster;
					end_pos = i-1;
					int dau1 = -1; int dau2 = -1; int dau3 = 99999;
					for(int j=srt_pos; j<=end_pos; ++j){if(py_rearrange[j][2] >= 0){
						if(     dau1 == -1){dau1 = py_rearrange[j][2];}
						else if(dau2 == -1){dau2 = py_rearrange[j][2];}
						else{if(py_rearrange[j][2]<dau2){dau3=dau2; dau2=py_rearrange[j][2];} else if(py_rearrange[j][2]<dau3){dau3=py_rearrange[j][2];}}
					}}
					if(dau3==99999){dau3=-1;}
					for(int j=srt_pos; j<=end_pos; ++j){py_rearrange[j][3]=dau1; py_rearrange[j][4]=dau2; py_rearrange[j][5]=dau3; py_rearrange[j][6]=ncluster;}
					end_pos = -1; srt_pos = -1; prev_id = -1;
					if(py_rearrange[i][1] >= 0){srt_pos = i; prev_id = py_rearrange[i][1];}
				}
			}
		}
		if(srt_pos >= 0){
			++ncluster;
			end_pos = py_rearrange.size()-1;
			int dau1 = -1; int dau2 = -1; int dau3 = -1;
			for(int j=srt_pos; j<=end_pos; ++j){if(py_rearrange[j][2] >= 0){
				if(     dau1 == -1){dau1 = py_rearrange[j][2];}
				else if(dau2 == -1){dau2 = py_rearrange[j][2];}
				else{if(py_rearrange[j][2]<dau2){dau3=dau2; dau2=py_rearrange[j][2];} else if(py_rearrange[j][2]<dau3){dau3=py_rearrange[j][2];}}
			}}
			if(dau3==99999){dau3=-1;}
			for(int j=srt_pos; j<=end_pos; ++j){py_rearrange[j][3]=dau1; py_rearrange[j][4]=dau2;  py_rearrange[j][5]=dau3; py_rearrange[j][6]=ncluster;}
		}
		
		//now for each element in py_re, we have:  0-num, 1-par, 2-dau, 3-dau1, 4-dau2, 5-dau3, 6-clusternum
			//where dau1, dau2, and dau3 are >= 0 for all partons in a cluster where 1(2)(or3) partons have daughters
			//if dau3 is set, it means that the partons with dau1 and dau2 have to be at the ends of the cluster
			//and any other daughters (3+) have to instead have a fake generation added between them and their actual daughters
		//now, want to sort the clusters so that if two clusters have dau1, and one has a dau2, the single daughter one goes first
		//but, want to sort them so that if two clusters have dau2, the double daughter one goes first
		//[dau1][dau1/dau2][dau2][dau3][dau3][dau4][dau4/dau5][dau5][dau6][dau6]
		std::stable_sort(py_rearrange.begin(), py_rearrange.end(),
			[](const std::vector<int>& vec1, const std::vector<int>& vec2){
				if(     (vec1[3] == -1) && (vec2[3] == -1)){return false;}
				else if((vec1[3] == -1) && (vec2[3] != -1)){return true;}
				else if((vec1[3] != -1) && (vec2[3] == -1)){return false;}
				else if((vec1[3] != -1) && (vec2[3] != -1)){
					if(     (vec1[4] == -1) && (vec2[4] == -1)){return vec1[3]<vec2[3];}
					else if((vec1[4] == -1) && (vec2[4] != -1)){
						if(     vec1[3] == vec2[3]){return true;}
						else if(vec1[3] == vec2[4]){return false;}
						else{return vec1[3]<vec2[3];}
					}
					else if((vec1[4] != -1) && (vec2[4] == -1)){
						if(     vec1[3] == vec2[3]){return false;}
						else if(vec1[4] == vec2[3]){return true;}
						else{return vec1[3]<vec2[3];}
					}
					else if((vec1[4] != -1) && (vec2[4] != -1)){
						if((vec1[3] == vec2[3]) || (vec1[3] == vec2[4]) || (vec1[4] == vec2[3]) || (vec1[4] == vec2[4])){return false;}
						else{return vec1[3]<vec2[3];}
					}
				}
				return false; //should never make it here, but just in case
			}
		);
		
		//now we need to sort within the clusters to make sure that the 'dual mothers' sit next to one another.
		//we also need to catch any dau3 entries afterwards - these will correspond to a 'dual mother' that *can't* be paired
			//these will be fixed by making each of those unpairable dual mothers into a mother of another fake mother
			//each of those fake mothers will replace the corresponding mother of the gluon chain (or other daughter)
			//will need to update the fake mothers of each of those daughters to point to the new fake mothers
		srt_pos = -1; end_pos = -1; prev_id = -1;
		bool first = true;
		std::vector<int> clusterstofix;
		for(int i=0; i<py_rearrange.size(); ++i){
			if(py_rearrange[i][1] != prev_id){
				if(srt_pos == -1){srt_pos = i; prev_id = py_rearrange[i][6];}
				else{
					end_pos = i-1;
					if(     (py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] == -1) &&  first){
						std::stable_sort(py_rearrange.begin()+srt_pos, py_rearrange.begin()+end_pos+1,
						[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]<vec2[2];});
						first = false;
					}
					else if((py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] == -1) && !first){
						std::stable_sort(py_rearrange.begin()+srt_pos, py_rearrange.begin()+end_pos+1,
						[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]>vec2[2];});
						first = true;
					}
					else if((py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] != -1)){
						std::stable_sort(py_rearrange.begin()+srt_pos,   py_rearrange.begin()+end_pos+1,
						[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]<vec2[2];});
						std::stable_sort(py_rearrange.begin()+srt_pos+1, py_rearrange.begin()+end_pos+1,
						[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]>vec2[2];});
						if(py_rearrange[srt_pos][5] != -1){clusterstofix.push_back(py_rearrange[srt_pos][6]);}
					}
					end_pos = -1; srt_pos = -1; prev_id = -1;
					if(py_rearrange[i][6] >= 0){srt_pos = i; prev_id = py_rearrange[i][6];}
				}
			}
		}
		if(srt_pos >= 0){
			end_pos = py_rearrange.size()-1;
			if(     (py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] == -1) &&  first){
				std::stable_sort(py_rearrange.begin()+srt_pos,py_rearrange.begin()+end_pos+1,[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]<vec2[2];});
				first = false;
			}
			else if((py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] == -1) && !first){
				std::stable_sort(py_rearrange.begin()+srt_pos,py_rearrange.begin()+end_pos+1,[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]>vec2[2];});
				first = true;
			}
			else if((py_rearrange[srt_pos][3] != -1) && (py_rearrange[srt_pos][4] != -1)){
				std::stable_sort(py_rearrange.begin()+srt_pos,   py_rearrange.begin()+end_pos+1,
				[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]<vec2[2];});
				std::stable_sort(py_rearrange.begin()+srt_pos+1, py_rearrange.begin()+end_pos+1,
				[](const std::vector<int>& vec1, const std::vector<int>& vec2){return vec1[2]>vec2[2];});
				if(py_rearrange[srt_pos][5] != -1){clusterstofix.push_back(py_rearrange[srt_pos][6]);}
			}
		}
		
		//clusters with entries to fix are put into 'clusterstofix' (these have dau3 or > that need fake daughters that will act as dual mothers for the 'previous' daughters)
		//will also add in a fake daughter for the corresponding 'dual mother' which may or may not be in a cluster that needs to be fixed.
			//will fix both in one go, then make sure to check any subsequent 'clusterstofix' to make sure we didn't already fix it.
		std::vector<bool> fixed; for(int i=0;i<py_rearrange.size();++i){fixed.push_back(false);}
		while(clusterstofix.size() > 0){
			ncluster = clusterstofix.back(); clusterstofix.pop_back();
			//pulling ALL partons in this cluster that need to have a fake generation added
			std::vector<int> ptns_addgen; std::vector<int> dau_tofix;
			for(int i=0;i<py_rearrange.size();++i){
				if(!fixed[i] && (py_rearrange[i][6]==ncluster) && (py_rearrange[i][2]!=-1) && (py_rearrange[i][2]!=py_rearrange[i][3]) && (py_rearrange[i][2]!=py_rearrange[i][4])){
					ptns_addgen.push_back(i); dau_tofix.push_back(py_rearrange[i][2]);
				}
			}
			//going through the rest of the string and pulling any parton that shares a daughter with ANY of the partons in the cluster that needed to be fixed
			for(int i=0;i<py_rearrange.size();++i){
				for(int j=0;j<dau_tofix.size();++j){
					//shouldn't need to check if these have been fixed, but we'll do it anyways?
					if(!fixed[i] && (py_rearrange[i][6] != ncluster) && (py_rearrange[i][2] == dau_tofix[j])){ptns_addgen.push_back(i);}
				}
			}
			//now, each of the partons in ptns_addgen need a single fake daughter added
			//these need to be added in PAIRS(could *technically* handle more, but shouldn't be more, or less, than in pairs) that both point to the shared daughters
			for(int i=0;i<dau_tofix.size();++i){
				//grabbing the partons for this daughter
				std::vector<int> ptns_fixnow; for(int j=0;j<py_rearrange.size();++j){if(py_rearrange[j][2] == dau_tofix[i]){ptns_fixnow.push_back(j);}}
				std::vector<int> added_ptns;
				
				//the partons in ptns_fixnow need to each have a fake daughter, all added next to each other at the end of str_completed AND py_rearrange
				//the daughter entries of the partons in *fixnow need to be fixed (dau1 set to its corresponding *new* fake daughter, dau2 set to -1)
				//the mother entries of the ALL previous daughters need to be set to the *new* fake daughters
					//(all the daughters have par1 set to the first mother, and par2 set to the last)
				//the added fake daughters need to have a single mother; it's corresponding *fixnow parton.
				//the added fake daughters ALL need to be set to have the first daughter and last daughter set as the *fixnow partons did
				//the added fake daughters will NOT have the py_rearrange array entries > 2 set correctly (0,1 are set correctly)
				for(int j=0;j<ptns_fixnow.size();++j){
					//adding the fake daughter
					HHparton fakedau  = str_completed[py_rearrange[ptns_fixnow[j]][0]]; fakedau.PY_stat(-21); /*fakedau.PY_stat = -11;*/
					fakedau.PY_par1(               py_rearrange[ptns_fixnow[j]][0] ); fakedau.PY_par2(-1);
					fakedau.PY_dau1( str_completed[py_rearrange[ptns_fixnow[j]][0]].PY_dau1() ); fakedau.PY_dau2( str_completed[py_rearrange[ptns_fixnow[j]][0]].PY_dau2() );
					str_completed.add(fakedau);
					std::vector<int> py_tmp; py_tmp.push_back(str_completed.num()-1); py_tmp.push_back(fakedau.PY_par1()); py_tmp.push_back(fakedau.PY_dau1());
					py_tmp.push_back(-1); py_tmp.push_back(-1); py_tmp.push_back(-1); py_tmp.push_back(0); py_rearrange.push_back(py_tmp);
					added_ptns.push_back(str_completed.num()-1);
					
					//better to just fix the 'mother' here and now
					str_completed[py_rearrange[ptns_fixnow[j]][0]].PY_dau1( str_completed.num()-1 ); str_completed[py_rearrange[ptns_fixnow[j]][0]].PY_dau2(-1);
				}
				
				//we've added the *new* daughter, fixed the 'mother', set the new fake daughters correctly.
				//finally, repair the previous daughter(s)
				if(str_completed[added_ptns[0]].PY_dau2() > -1){
					for(int j=0;j<str_completed.num()-added_ptns.size();++j){for(int k=0;k<added_ptns.size();++k){if(str_completed[j].PY_par1() == added_ptns[k]){
						str_completed[j].PY_par1( added_ptns[0] );
						str_completed[j].PY_par2( added_ptns[added_ptns.size()-1] );
						break;
					}}}
				}
				else{
					str_completed[str_completed[added_ptns[0]].PY_dau1()].PY_par1( added_ptns[0] );
					str_completed[str_completed[added_ptns[0]].PY_dau1()].PY_par2( added_ptns[added_ptns.size()-1] );
				}
				
			}
			
		}
		
		//now, we use the py_rearrange vector to:
			//rearrange str_completed - since py_re... was sorted wrt parents, this will cluster all daughters together properly
			//relabel all !PY_parX! entries for all partons in str_completed (daughters may have been disturbed - we'll fix this after)
		//we'll do the rearrangement by making a temporary parton collection, writing str_completed partons to it in order, then overwriting str_completed with the new collection
		parton_collection temp_collection;
		for(int i=0; i<py_rearrange.size(); ++i){temp_collection.add(str_completed[py_rearrange[i][0]]);}
		str_completed = temp_collection;
		
		//py_re returns old value from new pos 'i', lookup returns new value from old pos 'i'
		//use this to relabel PY_par entries
		std::vector<int> lookup;
		for(int i=0; i<py_rearrange.size(); ++i){for(int j=0;j<py_rearrange.size();++j){if(py_rearrange[j][0] == i){lookup.push_back(j);break;}}}
		for(int i=0; i<str_completed.num(); ++i){
			if(str_completed[i].PY_par1() >= 0){str_completed[i].PY_par1( lookup[str_completed[i].PY_par1()] );}
			if(str_completed[i].PY_par2() >= 0){str_completed[i].PY_par2( lookup[str_completed[i].PY_par2()] );}
		}
		
		//lastly, need to fix PY_dau entries; do this by going down the list and fix daughter entries for the mother of each cluster of daughters (or both mothers)
		int m_srt = -1; int m_end = -1; int prev_moth = -1;
		for(int i=0; i<str_completed.num(); ++i){
			if(str_completed[i].PY_par1() != prev_moth){
				if(m_srt == -1){m_srt = i; prev_moth = str_completed[i].PY_par1();}
				else{
					m_end = i-1;
					str_completed[prev_moth].PY_dau1(m_srt);
					str_completed[prev_moth].PY_dau2(m_end);
					if(m_srt == m_end){str_completed[prev_moth].PY_dau2(-1);}
					if(str_completed[m_end].PY_par2() >= 0){
						str_completed[str_completed[m_end].PY_par2()].PY_dau1(m_srt);
						str_completed[str_completed[m_end].PY_par2()].PY_dau2(m_end);
						if(m_srt == m_end){str_completed[str_completed[m_end].PY_par2()].PY_dau2(-1);}
					}
					m_end = -1; m_srt = -1; prev_moth = -1;
					if(str_completed[i].PY_par1() >= 0){m_srt = i; prev_moth = str_completed[i].PY_par1();}
				}
			}
		}
		
		//and just in case, we'll make sure that the string didn't end on a daughter cluster (which if there's more than 0 clusters, it will!)
		if(m_srt >= 0){
			m_end = str_completed.num()-1;
			str_completed[prev_moth].PY_dau1(m_srt);
			str_completed[prev_moth].PY_dau2(m_end);
			if(m_srt == m_end){str_completed[prev_moth].PY_dau2(-1);}
			if(str_completed[m_end].PY_par2() >= 0){
				str_completed[str_completed[m_end].PY_par2()].PY_dau1(m_srt);
				str_completed[str_completed[m_end].PY_par2()].PY_dau2(m_end);
				if(m_srt == m_end){str_completed[str_completed[m_end].PY_par2()].PY_dau2(-1);}
			}
		}
		
		//enforcing E/P conservation for PYTHIA 'fake' history particles (to allow for event:check to be set true for PYTHIA)
		bool EP_conserved = false;
		while(!EP_conserved){
			EP_conserved = true;
			for(int i=0; i<str_completed.num(); ++i){
				if(str_completed[i].PY_stat() >= 0){continue;}
				FourVector P_new(0.,0.,0.,0.); FourVector pos_new(0.,0.,0.,0.);
				int jmax = (str_completed[i].PY_dau2() > str_completed[i].PY_dau1()) ? str_completed[i].PY_dau2() : str_completed[i].PY_dau1();
				for(int j=str_completed[i].PY_dau1(); j<jmax+1; ++j){
					double n = double(j-str_completed[i].PY_dau1())+1.;
					P_new.Set(P_new.x()+str_completed[j].px(),P_new.y()+str_completed[j].py(),P_new.z()+str_completed[j].pz(),P_new.t()+str_completed[j].e());
					pos_new.Set(pos_new.x()+(str_completed[j].x()-pos_new.x())/n,pos_new.y()+(str_completed[j].y()-pos_new.y())/n,pos_new.z()+(str_completed[j].z()-pos_new.z())/n,pos_new.t()+(str_completed[j].x_t()-pos_new.t())/n);
				}
				if((dif2(P_new,str_completed[i].P())+(P_new.t()-str_completed[i].e())*(P_new.t()-str_completed[i].e()) > 0.00000001/*0.0001^2*/) ||
					(dif2(pos_new,str_completed[i].pos())+(pos_new.t()-str_completed[i].x_t())*(pos_new.t()-str_completed[i].x_t()) > 0.00000001)){
					str_completed[i].P(P_new); str_completed[i].pos(pos_new);
					str_completed[i].mass( str_completed[i].e()*str_completed[i].e()
					  - str_completed[i].px()*str_completed[i].px() - str_completed[i].py()*str_completed[i].py() - str_completed[i].pz()*str_completed[i].pz() );
					str_completed[i].mass( (str_completed[i].mass() >= 0.) ? sqrt(str_completed[i].mass()) : -sqrt(-str_completed[i].mass()) ); //don't need mass >0 for reco now
					EP_conserved = false;
				}
			}
		}
		
		//need to go back though and fix the 'dual mother' quarks to the gluon chain (need to set values /2)
		//this assumes that there are no fake history particles with more than 2 mothers (will fail if this 'bad' scenario is true)
		for(int i=0; i<str_completed.num(); ++i){if((str_completed[i].PY_stat() == -21) && (str_completed[i].PY_dau2() > str_completed[i].PY_dau1())){
			str_completed[i].px(str_completed[i].px()/2.);
			str_completed[i].py(str_completed[i].py()/2.);
			str_completed[i].pz(str_completed[i].pz()/2.);
			str_completed[i].e( str_completed[i].e() /2.);
			str_completed[i].mass(str_completed[i].mass()/2.);
		}}
		
		//now need to reindex str_completed to prepare to add to final output
		for(int i=0; i<str_completed.num(); ++i){
			if(str_completed[i].PY_par1() >= 0){str_completed[i].PY_par1( str_completed[i].PY_par1()+SP_prepremn.num() );}
			if(str_completed[i].PY_par2() >= 0){str_completed[i].PY_par2( str_completed[i].PY_par2()+SP_prepremn.num() );}
			if(str_completed[i].PY_dau1() >= 0){str_completed[i].PY_dau1( str_completed[i].PY_dau1()+SP_prepremn.num() );}
			if(str_completed[i].PY_dau2() >= 0){str_completed[i].PY_dau2( str_completed[i].PY_dau2()+SP_prepremn.num() );}
		}
		
		//adding the current string to the output parton collection.
		SP_prepremn.add(str_completed);
		
		//now that this string is good, fill the next string
		start = par_i; prev_str = upd_str;
		str_completed.clear(); current_str.clear();
		if(par_i<SP_remnants.num()){current_str.add(SP_remnants[par_i]);}
	}
}

//function to hand partons/strings and hadron resonances (and various other color neutral and colored objects) to Pythia8
bool HybridHadronization::invoke_py(){
	
	Event& event = pythia.event;
	
	//should have been checked before call, but if there are no partons/hadrons to deal with, just exit without invoking pythia
	if(HH_pyremn.num() + HH_hadrons.num() == 0){return true;}
	
	//first things first, need to reindex py_remn; just need to increment ALL by one
	//also restoring original id for color octet particles (or other 'odd' colored particles)
	for(int i=0; i<HH_pyremn.num(); ++i){
		HH_pyremn[i].PY_par1(HH_pyremn[i].PY_par1()+1); HH_pyremn[i].PY_par2(HH_pyremn[i].PY_par2()+1);
		HH_pyremn[i].PY_dau1(HH_pyremn[i].PY_dau1()+1); HH_pyremn[i].PY_dau2(HH_pyremn[i].PY_dau2()+1);
		if(HH_pyremn[i].PY_origid() != 0){HH_pyremn[i].id( HH_pyremn[i].PY_origid() );}
	}
	
	bool need_hadronization = true; bool success = true;
	int attempt = 0;
	while(need_hadronization){
		//incrementing attempt number
		++attempt;
		
		//resetting PYTHIA event record, so that this event can be filled
		event.reset();
		
		//number of partons/hadrons/particles handed to pythia
		int size_input = HH_pyremn.num();
		
		//keeping track of partons from py_remn and hadrons into the event
		std::vector<int> eve_to_had; eve_to_had.push_back(0);
		
		//filling PYTHIA event record with the partons from this event
		const double mm_to_fm = 100000000000.0; const double fm_to_mm = 1./mm_to_fm;
		for(int i=0; i<HH_pyremn.num(); ++i){
			//to make sure mass is set appropriately (if E^2-P^2<0, then mass is too)
			double massnow = HH_pyremn[i].e()*HH_pyremn[i].e() - HH_pyremn[i].px()*HH_pyremn[i].px() - HH_pyremn[i].py()*HH_pyremn[i].py() - HH_pyremn[i].pz()*HH_pyremn[i].pz();
			massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
			
			//just in case... (needed for 'low' mass)
			bool recalcE = false; double Enow = HH_pyremn[i].e();
			if(     (std::abs(HH_pyremn[i].id()) == 1) && (std::abs(massnow) < 0.340)){massnow = (massnow >= 0.) ? 0.340 : -0.340; recalcE = true;}
			else if((std::abs(HH_pyremn[i].id()) == 2) && (std::abs(massnow) < 0.336)){massnow = (massnow >= 0.) ? 0.336 : -0.336; recalcE = true;}
			else if((std::abs(HH_pyremn[i].id()) == 3) && (std::abs(massnow) < 0.486)){massnow = (massnow >= 0.) ? 0.486 : -0.486; recalcE = true;}
			else if((std::abs(HH_pyremn[i].id()) == 4) && (std::abs(massnow) < 1.60 )){massnow = (massnow >= 0.) ? 1.60  : -1.60 ; recalcE = true;}
			else if((std::abs(HH_pyremn[i].id()) == 5) && (std::abs(massnow) < 4.80 )){massnow = (massnow >= 0.) ? 4.80  : -4.80 ; recalcE = true;}
			//else if((std::abs(HH_pyremn[i].id()) == 4) && (std::abs(massnow) < 1.55 )){massnow = (massnow >= 0.) ? 1.55  : -1.55 ; recalcE = true;}
			//else if((std::abs(HH_pyremn[i].id()) == 5) && (std::abs(massnow) < 4.73 )){massnow = (massnow >= 0.) ? 4.73  : -4.73 ; recalcE = true;}
			if(recalcE){
				Enow = massnow*massnow + HH_pyremn[i].px()*HH_pyremn[i].px() + HH_pyremn[i].py()*HH_pyremn[i].py() + HH_pyremn[i].pz()*HH_pyremn[i].pz();
				Enow = (Enow >= 0.) ? sqrt(Enow) : sqrt(-Enow);
			}
			
			//append( id, status, mother1, mother2, daughter1, daughter2, col, acol, px, py, pz, e, m)
			event.append(HH_pyremn[i].id(),HH_pyremn[i].PY_stat(),HH_pyremn[i].PY_par1(),HH_pyremn[i].PY_par2(),HH_pyremn[i].PY_dau1(),HH_pyremn[i].PY_dau2(),
			HH_pyremn[i].col(),HH_pyremn[i].acol(),HH_pyremn[i].px(),HH_pyremn[i].py(),HH_pyremn[i].pz(),Enow,massnow);
			event[event.size()-1].vProd(HH_pyremn[i].x(), HH_pyremn[i].y(), HH_pyremn[i].z(), HH_pyremn[i].x_t());
			eve_to_had.push_back(-i-1);
		}
		
		//adding in hadrons/leptons/other colorless particles too...
		for(int i=0; i<HH_hadrons.num(); ++i){if(!HH_hadrons[i].is_final()){
			//to make sure mass is set appropriately (if E^2-P^2<0, then mass is too)
			double massnow = HH_hadrons[i].e()*HH_hadrons[i].e() - 
			(HH_hadrons[i].px()*HH_hadrons[i].px() + HH_hadrons[i].py()*HH_hadrons[i].py() + HH_hadrons[i].pz()*HH_hadrons[i].pz());
			massnow = (massnow >= 0.) ? sqrt(massnow) : -sqrt(-massnow);
			
			//append(id, status, col, acol, px, py, pz, e, m) - 81 should be the status code for a primary hadron (81-86) produced by a hadronization process
			//event.append(HH_hadrons[i].id,81,0,0,HH_hadrons[i].px(),HH_hadrons[i].py(),HH_hadrons[i].pz(),HH_hadrons[i].e(),HH_hadrons[i].mass);
			event.append(HH_hadrons[i].id(),81,0,0,HH_hadrons[i].px(),HH_hadrons[i].py(),HH_hadrons[i].pz(),HH_hadrons[i].e(),massnow);
			event[event.size()-1].vProd(HH_hadrons[i].x()*fm_to_mm, HH_hadrons[i].y()*fm_to_mm, HH_hadrons[i].z()*fm_to_mm, HH_hadrons[i].x_t()*fm_to_mm);
			++size_input; eve_to_had.push_back(i+1);
		}}
		
		//make PYTHIA hadronize this event, if it fails then retry N=10 times... (PYTHIA can and will rarely fail, without concern)
		//if this fails more than 10 times, we may retry this event starting back before recombination (some number of times)
		if(!pythia.next()){
		  if(attempt > 9){need_hadronization = false; success = false; break;}
		  continue;
		}
		
		//this event has been successfully hadronized (hopefully), put all final state particles into hadrons; set is_final=true for these
		//deal with setting parent (and position?) data later...
		//putting in a rapidity window here - remove this later, and do the rapidity check elsewhere?
		for(int i=1; i<event.size(); ++i){if(event[i].isFinal() /*&& (std::abs(event[i].y()) <= 1.)*/){
			HHhadron hadout; //a number of other tags need to be set from the parent (either hadron or string):
			hadout.is_final(true); hadout.id(event[i].id()); hadout.mass(event[i].m());
			hadout.px(event[i].px());            hadout.py(event[i].py());            hadout.pz(event[i].pz());            hadout.e(event[i].e());
			hadout.x(event[i].xProd()*mm_to_fm); hadout.y(event[i].yProd()*mm_to_fm); hadout.z(event[i].zProd()*mm_to_fm); hadout.x_t(event[i].tProd()*mm_to_fm);
			
			//since using inbuilt pythia mother/daughter functions will segfault 'occasionally', going to code it in by hand.
			//this could probably be done more efficiently, but it's good enough for now...
			std::vector<int> mothers;
			//using a stack system to fill mothers
			std::vector<int> stack;
			//filling stack with initial mothers
			if((event[i].mother1() < event[i].mother2()) && (event[i].mother1() > 0) && (std::abs(event[i].status()) >= 81) && (std::abs(event[i].status()) <= 86)){
				for(int ipar=event[i].mother1();ipar<=event[i].mother2();++ipar){stack.push_back(ipar);}
			}else if((event[i].mother2() > 0) && (event[i].mother1() != event[i].mother2())){
				stack.push_back(event[i].mother1()); stack.push_back(event[i].mother2());
			}else if(event[i].mother1() > 0){
				stack.push_back(event[i].mother1());
			}
			else{mothers.push_back(i);} //setting it as its own mother if there are no mothers (pythia didn't decay a directly input hadron...)
			
			//filling the stack with any valid mothers of the 'current' stack element
			//then we check the 'current' stack element to see if it's a valid mother (0<element<=n_input), if so, we write it to mothers
			while(stack.size() > 0){
				int current = stack.back(); stack.pop_back();
				
				if((event[current].mother1() < event[current].mother2()) && (event[current].mother1() > 0) &&
				(std::abs(event[current].status()) >= 81) && (std::abs(event[current].status()) <= 86)){
					for(int ipar=event[current].mother1();ipar<=event[current].mother2();++ipar){stack.push_back(ipar);}
				}else if((event[current].mother2() > 0) && (event[current].mother1() != event[current].mother2())){
					stack.push_back(event[current].mother1()); stack.push_back(event[current].mother2());
				}else if(event[current].mother1() > 0){
					stack.push_back(event[current].mother1());
				}
				
				if((current > 0) && (current <= size_input)){mothers.push_back(current);}
			}
			
			//just in case...
			if(mothers.size() == 0){mothers.push_back(i);}
			
			//sorting and removing duplicate entries
			std::sort(mothers.begin(), mothers.end()); mothers.erase(std::unique(mothers.begin(), mothers.end()), mothers.end());
			
			//first, using mothers to determine if this hadron was formed via recombination, or by string fragmentation
			if(mothers[0] <= HH_pyremn.num()){hadout.is_strhad(true);} else{hadout.is_recohad(true);}
			
			//lastly, using mothers (except fake) in original input to determine if this a shower-shower or shower-thermal hadron
			bool is_therm(false);
			for(int ipar=0;ipar<mothers.size();++ipar){
				if(mothers[ipar] <= HH_pyremn.num()){
					if(HH_pyremn[mothers[ipar]-1].orig() != -1){hadout.add_par(mothers[ipar]-1); if(HH_pyremn[mothers[ipar]-1].is_thermal()){is_therm=true;}}
				}
				else if(mothers[ipar] <= size_input){ //shouldn't actually need to check, but doing so just in case.
					hadout.parh( eve_to_had[mothers[ipar]]-1 ); if(HH_hadrons[hadout.parh()].is_shth()){is_therm=true;}
				}
			}
			if(is_therm){hadout.is_shth(true);}else{hadout.is_shsh(true);}
			
			//the mother procedure might skip some partons if there are junctions involved
			//this can be 'repaired' by taking a 'mother' parton, then checking over all the partons in its string! (both adding to parents / checking if thermal)
			//this is done in hadronization calling function, after this invoke_py function is finished
			HH_hadrons.add(hadout);
		}}
		
		need_hadronization = false;
	}
	
	return success;
}
