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
#include "ThermPtnSampler.h"
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
//#include <cmath>

using namespace Jetscape;
using namespace Pythia8;

// Register the module with the base class
RegisterJetScapeModule<HybridHadronization> HybridHadronization::reg("HybridHadronization");

// Initialize static helper here
Pythia8::Pythia HybridHadronization::pythia ("IntentionallyEmpty",false);

//RNG - Mersenne Twist - 64 bit
//std::mt19937_64 eng(std::random_device{}());
//std::mt19937_64 eng(1);
//returns a random number between 0 and 1, based on above engine
double HybridHadronization::ran() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(eng);}

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

  //tinyxml2::XMLElement *hadronization= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("JetHadronization" );
  tinyxml2::XMLElement *hadronization = GetXMLElement({"JetHadronization"});

  if ( !hadronization ) {
    JSWARN << "Couldn't find tag Jet Hadronization";
    throw std::runtime_error ("Couldn't find tag Jet Hadronization");
  }
  if (hadronization) {
    string s = hadronization->FirstChildElement( "name" )->GetText();
	//std::string s = GetXMLElementText({"JetHadronization", "name"});
	JSDEBUG << s << " to be initialized ...";
	JSINFO<<"Initialize Hybrid Hadronization ...";

    JSDEBUG<<"Initialize HybridHadronization";
    VERBOSE(8);
	
	maxE_level    = 3;			//maximum energy level considered for the recombination (was 3 in recent fortran code, prev. set to 8)
	gmax          = 1.25;		//maximum allowed mass of the gluon (for q-qbar split), in GeV
	xmq           = 0.338; //0.33;		//light quark mass, in GeV
	xms           = 0.486; //0.5;		//strange quark mass, in GeV
	hbarc         = 0.197327;	// GeV*fm - maybe just set this as a constant in common?
	dist2cut      = 25.;		//maximum distance [fm] squared for recombination (involving thermal partons) - in lab frame
	sh_recofactor = 1./3.;		//suppression/enhancement factor for shower-shower recombination
	th_recofactor = 1./3.;		//suppression/enhancement factor for shower-thermal recombination
	attempts_max  = 10;			//maximum number of failed attempts to hadronize a single event before we give up.
	p_fake        = 0.;			//momentum used for fake parton, if needed
	rand_seed     = 0;			//seed for RNGs used - 0 means use a randomly determined random seed (from system time or std::random_device{}())
	
	
	//xml read in to alter settings...
	double xml_doublein = -1.; int xml_intin = -1; unsigned int xml_uintin = std::numeric_limits<unsigned int>::max();
	
	//hadronization->FirstChildElement("eCMforHadronization")->QueryDoubleText(&xml_doublein);
	xml_doublein = GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
	if(xml_doublein >= 0.){p_fake = xml_doublein;} xml_doublein = -1.;
	
	//hadronization->FirstChildElement("thermreco_distmax")->QueryDoubleText(&xml_doublein);
	xml_doublein = GetXMLElementDouble({"JetHadronization", "thermreco_distmax"});
	if(xml_doublein >= 0.){dist2cut = xml_doublein*xml_doublein;} xml_doublein = -1.;
	
	//hadronization->FirstChildElement("shower_recofactor")->QueryDoubleText(&xml_doublein);
	xml_doublein = GetXMLElementDouble({"JetHadronization", "shower_recofactor"});
	if(xml_doublein >= 0.){sh_recofactor = xml_doublein;} xml_doublein = -1.;
	
	//hadronization->FirstChildElement("thermal_recofactor")->QueryDoubleText(&xml_doublein);
	xml_doublein = GetXMLElementDouble({"JetHadronization", "thermal_recofactor"});
	if(xml_doublein >= 0.){th_recofactor = xml_doublein;} xml_doublein = -1.;
	
	//hadronization->FirstChildElement("reco_Elevelmax")->QueryIntText(&xml_intin);
	xml_intin = GetXMLElementInt({"JetHadronization", "reco_Elevelmax"});
	if(xml_intin >= 0){maxE_level = xml_intin;} xml_intin = -1;
	
	// random seed
	// xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
	//tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
	tinyxml2::XMLElement *RandomXmlDescription = GetXMLElement({"Random"});
	if ( RandomXmlDescription ){
		tinyxml2::XMLElement *xmle; xmle = RandomXmlDescription->FirstChildElement( "seed" );
		xmle->QueryUnsignedText(&xml_uintin);
		if(xml_uintin < std::numeric_limits<unsigned int>::max()){rand_seed = xml_uintin;} xml_uintin = std::numeric_limits<unsigned int>::max();
	}
	else{
		JSWARN << "No <Random> element found in xml, seeding to 0";
	}
	VERBOSE(7) <<"Seeding PYTHIA(hadronization) to "<< rand_seed;
	
	//not sure if we can seed with a negative integer...
	if(rand_seed != 0){eng.seed(rand_seed);}
	//else{eng.seed(std::random_device{}());}
	else{ //seeding the mt19937_64 object 'eng' PROPERLY!
		std::random_device rd;
		std::array<int,std::mt19937_64::state_size> seedarray; std::generate_n(seedarray.data(), seedarray.size(), std::ref(rd));
		std::seed_seq seeds(std::begin(seedarray), std::end(seedarray)); eng.seed(seeds);
	}
	
	
	//Since PYTHIA has no spacetime information, we shouldn't use recombination as it is necessary to calculate recombination probabilities
	//later, this will instead update to set a flag to attempt to artificially generate this information
	//for now, we just print a warning - but still try to run recombination.  It will just be unphysically enhanced, esp. for certain configurations
	//tinyxml2::XMLElement *XmlPythiaGun=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hard" )->FirstChildElement("PythiaGun");
/*	tinyxml2::XMLElement *XmlPythiaGun = GetXMLElement({"Hard", "PythiaGun"});
	bool PYgun_FSRon = false;
	if ( XmlPythiaGun ){
		tinyxml2::XMLElement *xmle; xmle = XmlPythiaGun->FirstChildElement( "FSR_on" );
*****		xmle->QueryIntText(&xml_intin);
		if(xml_intin == 1){PYgun_FSRon = true;} xml_intin = -1;
	}
	if(PYgun_FSRon && (sh_recofactor > 0.0000000001)){JSWARN << "Recombination with a PYTHIA FSR shower is not fully implemented.";}
*/	
	xml_intin = GetXMLElementInt({"Hard", "PythiaGun", "FSR_on"});
	if(xml_intin && (sh_recofactor > 0.0000000001)){JSWARN << "Recombination with a PYTHIA FSR shower is not fully implemented.";}
	xml_intin = -1;
	
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
	SigPi2  = SigM2_calc(R2chg_Pi,  Qm_ud, Qm_ud, chg_d, chg_u);
	SigPhi2 = SigM2_calc(R2chg_Phi, Qm_s,  Qm_s,  chg_d, -chg_d)*(2./3.); //normalizing
	SigK2   = SigM2_calc(R2chg_K,   Qm_s,  Qm_ud, chg_d, chg_u);
	SigJpi2 = SigM2_calc(R2chg_Jpi, Qm_c,  Qm_c,  chg_u, -chg_u)*(4./3.); //normalizing
	SigDs2  = SigM2_calc(R2chg_Ds,  Qm_c,  Qm_s,  chg_u, chg_d);
	SigD2   = SigM2_calc(R2chg_D,   Qm_c,  Qm_ud, chg_u, chg_d);
	SigUps2 = SigM2_calc(R2chg_Ups, Qm_b,  Qm_b,  chg_d, -chg_d)*(2./3.); //normalizing
	SigBc2  = SigM2_calc(R2chg_Bc,  Qm_b,  Qm_c,  chg_d, chg_u);
	SigB2   = SigM2_calc(R2chg_B,   Qm_b,  Qm_ud, chg_d, chg_u); // (treating B_s as B)

	//baryon width calculations (r2) - recalc if r2chg is changed on command line...
	//light/strange baryons
	SigNucR2 = SigBR2_calc(R2chg_Nuc, Qm_ud, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	SigNucL2 = SigBL2_calc(SigNucR2,  Qm_ud, Qm_ud, Qm_ud);
	SigOmgR2 = SigBR2_calc(R2chg_Omg, Qm_s,  Qm_s,  Qm_s,  chg_d, chg_d, chg_d);
	SigOmgL2 = SigBL2_calc(SigOmgR2,  Qm_s,  Qm_s,  Qm_s );
	SigXiR2  = SigBR2_calc(R2chg_Xi,  Qm_s,  Qm_s,  Qm_ud, chg_d, chg_d, chg_d);
	SigXiL2  = SigBL2_calc(SigXiR2,   Qm_s,  Qm_s,  Qm_ud);
	SigSigR2 = SigBR2_calc(R2chg_Sig, Qm_s,  Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	SigSigL2 = SigBL2_calc(SigSigR2,  Qm_s,  Qm_ud, Qm_ud);

	//charm baryons
	SigOcccR2 = SigBR2_calc(R2chg_Occc, Qm_c, Qm_c,  Qm_c,  chg_d, chg_d, chg_d); // ! maybe need to normalize? (just setting all to -1/3 for now)
	SigOcccL2 = SigBL2_calc(SigOcccR2,  Qm_c, Qm_c,  Qm_c );
	SigOccR2  = SigBR2_calc(R2chg_Occ,  Qm_c, Qm_c,  Qm_s,  chg_u, chg_u, chg_d);
	SigOccL2  = SigBL2_calc(SigOccR2,   Qm_c, Qm_c,  Qm_s );
	SigXiccR2 = SigBR2_calc(R2chg_Xicc, Qm_c, Qm_c,  Qm_ud, chg_u, chg_u, chg_d);
	SigXiccL2 = SigBL2_calc(SigXiccR2,  Qm_c, Qm_c,  Qm_ud);
	SigOcR2   = SigBR2_calc(R2chg_Oc,   Qm_c, Qm_s,  Qm_s,  chg_d, chg_d, chg_d); // ! setting all quark charges to -1/3
	SigOcL2   = SigBL2_calc(SigOcR2,    Qm_c, Qm_s,  Qm_s );
	SigXicR2  = SigBR2_calc(R2chg_Xic,  Qm_c, Qm_s,  Qm_ud, chg_u, chg_d, chg_u);
	SigXicL2  = SigBL2_calc(SigXicR2,   Qm_c, Qm_s,  Qm_ud);
	SigSigcR2 = SigBR2_calc(R2chg_Sigc, Qm_c, Qm_ud, Qm_ud, chg_u, chg_d, chg_u);
	SigSigcL2 = SigBL2_calc(SigSigcR2,  Qm_c, Qm_ud, Qm_ud);

	//bottom baryons
	SigObbbR2 = SigBR2_calc(R2chg_Obbb, Qm_b, Qm_b,  Qm_b,  chg_d, chg_d, chg_d);
	SigObbbL2 = SigBL2_calc(SigObbbR2,  Qm_b, Qm_b,  Qm_b );
	SigObbcR2 = SigBR2_calc(R2chg_Obbc, Qm_b, Qm_b,  Qm_c,  chg_d, chg_d, chg_d); // ! setting all quark charges to -1/3
	SigObbcL2 = SigBL2_calc(SigObbcR2,  Qm_b, Qm_b,  Qm_c );
	SigObbR2  = SigBR2_calc(R2chg_Obb,  Qm_b, Qm_b,  Qm_s,  chg_d, chg_d, chg_d);
	SigObbL2  = SigBL2_calc(SigObbR2,   Qm_b, Qm_b,  Qm_s );
	SigXibbR2 = SigBR2_calc(R2chg_Xibb, Qm_b, Qm_b,  Qm_ud, chg_d, chg_d, chg_d);
	SigXibbL2 = SigBL2_calc(SigXibbR2,  Qm_b, Qm_b,  Qm_ud);
	SigObccR2 = SigBR2_calc(R2chg_Obcc, Qm_b, Qm_c,  Qm_c,  chg_d, chg_u, chg_u);
	SigObccL2 = SigBL2_calc(SigObccR2,  Qm_b, Qm_c,  Qm_c );
	SigObcR2  = SigBR2_calc(R2chg_Obc,  Qm_b, Qm_c,  Qm_s,  chg_d, chg_d, chg_d); // ! flipping c quark charge (all to -1/3)
	SigObcL2  = SigBL2_calc(SigObcR2,   Qm_b, Qm_c,  Qm_s );
	SigXibcR2 = SigBR2_calc(R2chg_Xibc, Qm_b, Qm_c,  Qm_ud, chg_d, chg_u, chg_u);
	SigXibcL2 = SigBL2_calc(SigXibcR2,  Qm_b, Qm_c,  Qm_ud);
	SigObR2   = SigBR2_calc(R2chg_Ob,   Qm_b, Qm_s,  Qm_s,  chg_d, chg_d, chg_d);
	SigObL2   = SigBL2_calc(SigObR2,    Qm_b, Qm_s,  Qm_s );
	SigXibR2  = SigBR2_calc(R2chg_Xib,  Qm_b, Qm_s,  Qm_ud, chg_d, chg_d, chg_d);
	SigXibL2  = SigBL2_calc(SigXibR2,   Qm_b, Qm_s,  Qm_ud);
	SigSigbR2 = SigBR2_calc(R2chg_Sigb, Qm_b, Qm_ud, Qm_ud, chg_d, chg_u, chg_u);
	SigSigbL2 = SigBL2_calc(SigSigbR2,  Qm_b, Qm_ud, Qm_ud);

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
	
	//setting seed, or using random seed
	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = " + std::to_string(rand_seed));
	
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
	//pythia.readString("PartonVertex:setVertex = on");        //this might allow PYTHIA to keep track of partonic space-time information (default was for 'rope hadronization')
	
	//setting hadron color tags so that spacetime information can be reconstructed
	pythia.readString("StringFragmentation:TraceColours = on");
	
	//using QCD based color reconnection (original PYTHIA MPI based CR can't be used at hadron level)
	pythia.readString("ColourReconnection:reconnect = off");          //allowing color reconnections (should have been default on, but doing it here for clarity)
/*	pythia.readString("ColourReconnection:mode = 1");                //sets the color reconnection scheme to 'new' QCD based scheme (TODO: make sure this is better than (2)gluon move)
	pythia.readString("ColourReconnection:forceHadronLevelCR = on"); //allowing color reconnections for these constructed strings!
	//a few params for the QCD based color reconnection scheme are set below.
	pythia.readString("MultipartonInteractions:pT0Ref = 2.15");      //not sure if this is needed for this setup, but is part of the recommended 'default'
	pythia.readString("ColourReconnection:allowDoubleJunRem = off"); //default on - allows directly connected double junction systems to split into two strings
	pythia.readString("ColourReconnection:junctionCorrection = 1.15");
	pythia.readString("ColourReconnection:timeDilationMode = 3");    //allow reconnection if single pair of dipoles are in causal contact (maybe try 5 as well?)
	pythia.readString("ColourReconnection:timeDilationPar = 0.18");  //parameter used in causal interaction between strings (mode set above)(maybe try 0.073?)
*/
    // And initialize
    pythia.init();
	
	//setting up thermal partons...
	//tinyxml2::XMLElement *XmlMat=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("Matter");
/*	tinyxml2::XMLElement *XmlMat = GetXMLElement({"Eloss", "Matter"});
	bool brickptns = false; double brick_L = -1;
	if ( XmlMat ){
		tinyxml2::XMLElement *xmle; xmle = XmlMat->FirstChildElement( "brick_med" );
		xmle->QueryIntText(&xml_intin);
		if(xml_intin == 1){brickptns = true;} xml_intin = -1;
		if(brickptns){xmle = XmlMat->FirstChildElement( "brick_length" ); xmle->QueryDoubleText(&xml_doublein); if(xml_doublein >= 0.){brick_L = xml_doublein;} xml_doublein = -1.;}
	}
*/	
	bool brickptns = false; double brick_L = -1;
	xml_intin = GetXMLElementInt({"Eloss", "Matter", "brick_med"});
	if(xml_intin == 1){brickptns = true;} xml_intin = -1;
	xml_doublein = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});
	if(brickptns && (xml_doublein >= 0.)){brick_L = xml_doublein;} xml_doublein = -1.;
	if(brick_L < 0.){brick_L = 0.; brickptns = false;}
	
	if(brickptns){
		
		thermptnsampler brick; //creating a thermal brick
		brick.brick_LWo(2.*brick_L,2.*brick_L,2.*brick_L); //(double len_bri(z), double wid_bri(x,y), double offset_in(unused)) //Just setting as a 'normal' brick 'cube'
		//brick.brick_flow(0., 0., tanh(0.6)); //(double vx_in, double vy_in, double vz_in)
		brick.brick_flow(0.6, 0., 0.); //(double vx_in, double vy_in, double vz_in)
		brick.brick_seed(rand_seed); //(unsigned int seed_in)
		brick.samplebrick();
		
		JSINFO << "A " << brick_L << " fm brick was sampled, generating " << brick.nTot() << " partons (" << brick.th_nL() << " light, " << brick.th_nS() << " strange).";
		
		for(int ith=0; ith<brick.nTot(); ++ith){
			HHparton thparton;
			
			thparton.is_thermal(true); thparton.id(brick.th_pid(ith)); thparton.orig(1); //sh_parton.string_id(str)
			thparton.px(brick.th_px(ith)); thparton.py(brick.th_py(ith)); thparton.pz(brick.th_pz(ith)); thparton.e(  brick.th_e(ith));
			thparton.x( brick.th_x(ith) ); thparton.y( brick.th_y(ith) ); thparton.z( brick.th_z(ith) ); thparton.x_t(brick.th_t(ith));
			thparton.mass( thparton.e()*thparton.e() - thparton.px()*thparton.px() - thparton.py()*thparton.py() - thparton.pz()*thparton.pz() );
			thparton.mass( (thparton.mass() >= 0.) ? sqrt(thparton.mass()) : sqrt(-thparton.mass()) );
			thparton.pos_str(1);
			
			HH_thermal.add(thparton); //adding this parton to thermal collection
		}
		
	}
	//read in thermal partons, THEN do sibling setup...
	for(int i=0;i<HH_thermal.num();++i){
		HH_thermal[i].sibling(i); HH_thermal[i].string_id(-i);
		HH_thermal[i].is_used(true); HH_thermal[i].sibling( findthermalsibling(i, HH_thermal) ); HH_thermal[i].is_used(false);
	}
	
	//temp. outputting thermal partons to ensure they're sampled...
	std::cout << "\n\n\n";
	for(int i=0;i<HH_thermal.num();++i){
		std::cout << HH_thermal[i].id() << ", " << HH_thermal[i].px() << ", " << HH_thermal[i].py() << ", " << HH_thermal[i].pz() << ", " << HH_thermal[i].e() << ", " << 
		  HH_thermal[i].x() << ", " << HH_thermal[i].y() << ", " << HH_thermal[i].z() << ", " << HH_thermal[i].x_t() << "\n";
	}
	std::cout << "\n\n\n";
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
		if(shower.at(ishower).at(ipart)->pstat()<0){continue;}
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
  
	//checking to see if in brick, and if so then propagate all partons to hypersurface!
	//this really shouldn't have to be done; the partonic showering modules should *not* give partons inside the brick... but they will.
	//tinyxml2::XMLElement *XmlMat=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" )->FirstChildElement("Matter");
	//tinyxml2::XMLElement *XmlMat = GetXMLElement({"Eloss", "Matter"});
	bool in_brick = false; double brick_L = -1;
	(GetXMLElementInt({"Eloss", "Matter", "brick_med"}) > 0) ? in_brick = true : in_brick = false;
	if(in_brick){brick_L = GetXMLElementDouble({"Eloss", "Matter", "brick_length"});}
	if(brick_L < 0.){brick_L = 0.; in_brick = false;}
	if(in_brick){
		for(int i_sh=0; i_sh<HH_shower.num(); ++i_sh){ if(HH_shower[i_sh].x_t()<brick_L){
			JSINFO<<"Parton found inside brick at: " << HH_shower[i_sh].x() << ", " << HH_shower[i_sh].y() << ", " << HH_shower[i_sh].z() << ", " << HH_shower[i_sh].x_t();
			double t_dif = brick_L - HH_shower[i_sh].x_t();
			double vel[3]; vel[0]=HH_shower[i_sh].px()/HH_shower[i_sh].e(); vel[1]=HH_shower[i_sh].py()/HH_shower[i_sh].e(); vel[2]=HH_shower[i_sh].pz()/HH_shower[i_sh].e();
			HH_shower[i_sh].x(HH_shower[i_sh].x() + vel[0]*t_dif); HH_shower[i_sh].y(HH_shower[i_sh].y() + vel[1]*t_dif); HH_shower[i_sh].z(HH_shower[i_sh].z() + vel[2]*t_dif);
			HH_shower[i_sh].x_t(brick_L);
			JSINFO<<"Parton propagated to: " << HH_shower[i_sh].x() << ", " << HH_shower[i_sh].y() << ", " << HH_shower[i_sh].z() << ", " << HH_shower[i_sh].x_t();
		}}
	}
  
  int num_strings = 0;
  
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
	num_strings = HH_showerptns[HH_showerptns.num()-1].string_id();
	
	//running recombination
	recomb();
	
	//temporary workaround to force all formed hadrons to be final state - so pythia won't overwrite spacetime info
	for(int i=0; i<HH_hadrons.num(); ++i){HH_hadrons[i].is_final(true);}
	
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
		int stat = (HH_hadrons[iHad].is_recohad()) ? 999 : -999;
		int lab = 99;
		if(HH_hadrons[iHad].is_shth()){/*stat*=100;*/lab = -99;}
		//int had_parstr = HH_showerptns[HH_hadrons[iHad].par(0)].string_id();
		int idH = HH_hadrons[iHad].id(); double mH = HH_hadrons[iHad].mass();
		FourVector p(HH_hadrons[iHad].P()); FourVector x(HH_hadrons[iHad].pos());
		//hOut.push_back(std::make_shared<Hadron> (Hadron (0,idH,1,p,x,mH)));
		hOut.push_back(std::make_shared<Hadron> (Hadron (stat,idH,lab,p,x,mH)));
	}
  }
  JSINFO<<"#Showers hadronized together: " << shower.size() << " ( " << num_strings << " initial strings ). There are " <<
    hOut.size() << " hadrons and " << pOut.size() << " partons after Hybrid Hadronization";
  
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

  //parton list for treating thermal siblings for string repair
  parton_collection Extraparton;

  std::vector<int> fakepartoninfo;

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

//std::cout <<endl;
//std::cout <<"Below is Color information of all particles in the String (Col, Acol)"<<endl;
std::vector<int*> ColInfo3;
      for(int i=0; i<HH_showerptns.num(); i++) {
         int colinfo3[2] = {HH_showerptns[i].col() , HH_showerptns[i].acol()};
         ColInfo3.push_back(colinfo3);

        //std::cout <<" ( "<<HH_showerptns[i].string_id()<<" , "<<ColInfo3.at(i)[0]<<"   "<<ColInfo3.at(i)[1]<<"  ) ";
      }
//std::cout <<endl;
//std::cout <<endl;

//while running PythiaBrickTest, error{same col tags and acol tag for the different particles, somehow the cause would be in MATTER in determining the color tag} detected,
//So, If there are problem from the initial structure, correct it based on the color flow.
int temptag = 1;
for(int i = 0; i < HH_showerptns.num(); i++){
  for(int j = 0; j < HH_showerptns.num(); j++){
    if(HH_showerptns[i].col() != 0 && HH_showerptns[i].col() == HH_showerptns[j].col()){ //same col tag for different particles detected.
      for(int k = 0; k < HH_showerptns.num(); k++){
        if(HH_showerptns[i].col() == HH_showerptns[k].acol()){
          vector<int> Info;// location of i, location of k, temp color tag will be saved in this vector.
          Info.push_back(i);
          Info.push_back(k);
          Info.push_back(temptag);

          HH_showerptns[i].col(temptag);
          HH_showerptns[k].acol(temptag); //apply the change for the col, acol pair to preserve the color flow with least impact to the final result.
          temptag++;
          Info.clear();
        }
      }
    }
  }
}

for(int i = 0; i < HH_showerptns.num(); i++){
  for(int j = 0; j < HH_showerptns.num(); j++){
    if(HH_showerptns[i].acol() != 0 && HH_showerptns[i].acol() == HH_showerptns[j].acol()){ //same col tag for different particles detected.
      for(int k = 0; k < HH_showerptns.num(); k++){
        if(HH_showerptns[i].acol() == HH_showerptns[k].col()){
          vector<int> Info;// location of i, location of k, temp color tag will be saved in this vector.
          Info.push_back(i);
          Info.push_back(k);
          Info.push_back(temptag);

          HH_showerptns[i].acol(temptag);
          HH_showerptns[k].col(temptag); //apply the change for the col, acol pair to preserve the color flow with least impact to the final result.
          temptag++;
          Info.clear();
        }
      }
    }
  }
}

//std::cout <<endl;
//std::cout <<"Below is Color information of all particles in the Revised String (Col, Acol)"<<endl;
std::vector<int*> ColInfo2;
      for(int i=0; i<HH_showerptns.num(); i++) {
         int colinfo2[2] = {HH_showerptns[i].col() , HH_showerptns[i].acol()};
         ColInfo2.push_back(colinfo2);
        //std::cout <<" ( "<<HH_showerptns[i].string_id()<<" , "<<ColInfo2.at(i)[0]<<"   "<<ColInfo2.at(i)[1]<<"  ) ";
      }
//std::cout <<endl;
//std::cout <<endl;




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
			//HH_showerptns[i_pt].mass( 2.*xmq + (gmax - 2.*xmq)*ran() );	// gluon virtuality
			if(HH_showerptns[i_pt].mass()<2.*xmq+0.001){HH_showerptns[i_pt].mass(2.*xmq+0.001);}
			
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
//TODO : Make Color Recombination Vector for Meson
//std::cout <<endl<<endl;
//std::cout <<" Let's Check What We Do!!"<<endl;
    double MesonColRecomb[4] = {1, 0, 1/8, 1/9};
//std::cout <<"Color Recombination array below"<<endl;
//Check anti-color tags from all partons
//std::cout <<"This is Array of Acol Tags"<<endl;
std::vector<int> Acol;
for(int iacol=0;iacol<HH_showerptns.num();++iacol){Acol.push_back(HH_showerptns[iacol].acol());
     //std::cout <<Acol.at(iacol)<<"    ";
   }
//std::cout <<endl;
//std::cout <<"This is Array of Col Tags"<<endl;
std::vector<int> Col;
for(int icol=0;icol<HH_showerptns.num();++icol){Col.push_back(HH_showerptns[icol].col());
     //std::cout <<Col.at(icol)<<"    ";
   }
//std::cout <<endl;
//std::cout <<"Below is Color information of all particles in the String (Col, Acol)"<<endl;
std::vector<int*> ColInfo;
      for(int i=0; i<HH_showerptns.num(); i++) {
         int colinfo[2] = {HH_showerptns[i].col() , HH_showerptns[i].acol()};
         ColInfo.push_back(colinfo);
        //std::cout <<" ( "<<HH_showerptns[i].string_id()<<" , "<<ColInfo.at(i)[0]<<"   "<<ColInfo.at(i)[1]<<"  ) ";
      }
//std::cout <<endl;
//std::cout <<endl;




//do the same task for acol tag checking.


//TODO: We need to form Recombination matrix based on color tag
// Work flow

// make all color tags to be indices for the matrix (so the indice represent the color tag, don't need to classify anti or not)
std::vector<vector<int>> IndiceForCol1; // this vector will form vector of vectors(location in the string, anti or not, color tag)
std::vector<int> IndiceForCol2; // this vector is one component of vector above
std::vector<int> IndiceForColFin; // this vector is final vector that'll be used for process

for(int iIndice=0; iIndice < HH_showerptns.num(); iIndice++) {
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
} // As a result, two vectors with 3 component keep being added to bigger vector, now exclude same tag and zero
//dignostic measure
//std::cout <<endl<<" Below is list of all color tags " <<endl;
     //std::cout <<" ( ";
for (int icheck=0; icheck < IndiceForCol1.size(); icheck++) {
     //std::cout << IndiceForCol1.at(icheck).at(2) << " , ";
}
     //std::cout <<" ) " <<endl;



vector<int> tempcol;
// now put all the values in ColInfo into IndiceForColFin
for(int icol = 0; icol < IndiceForCol1.size(); icol++){
  tempcol.push_back(IndiceForCol1.at(icol).at(2));
}
if(tempcol.at(0) != 0){
IndiceForColFin.push_back(tempcol.at(0));
}
if(tempcol.at(1) != 0){
IndiceForColFin.push_back(tempcol.at(1));
}

for(int iclear1 = 0; iclear1 < tempcol.size(); iclear1++){
  bool foundsame = false;
  for(int iclear2 = 0; iclear2 < IndiceForColFin.size(); iclear2++){
    if(tempcol.at(iclear1) == IndiceForColFin.at(iclear2) || tempcol.at(iclear1) == 0){
    foundsame = true;
     }
    }
    if(!foundsame){
      IndiceForColFin.push_back(tempcol.at(iclear1));
    }
    foundsame = false;
  }

//last checking whether there is zero in the color tag list
for(int i = 0; i < IndiceForColFin.size(); i++){
  if(IndiceForColFin[i] == 0){
    std::vector<int>::iterator i1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), 0 );
    int distance = std::distance(IndiceForColFin.begin() , i1);
    //std::cout <<endl<<" zero is located at "<<distance<<endl;
    IndiceForColFin.erase(i1);
  }
}




//dignostic syntax
//std::cout <<endl<<" Below is list of valid color tags " <<endl;
     //std::cout <<" ( ";
for (int icheck=0; icheck < IndiceForColFin.size(); icheck++) {
     //std::cout << IndiceForColFin.at(icheck) << " , ";
}
     //std::cout <<" ) " <<endl;

//when we includes the partons from LBT, color tags from them are both zero in col and acol tags,
//Therefore, based on the maximum color tag from the parton with col , acol tags, reassign the color tags for LBT partons. As a result, one fake string should be formed based on these color Tags
//First, check the maximum color tags in the vector of valid color tags
int limit;
int maxtag;
maxtag = *max_element(IndiceForColFin.begin() , IndiceForColFin.end()); //this will be so important in dealing thermal partons, beacause they get color tag incremented from this max tag
limit = maxtag; // save initial maxtag for the future usage{It should be used for color reconnection matrix}
//std::cout <<endl<<" max tag is "<<maxtag<<endl;


// before determine the matrix, need to link color tag with the matrix indices( by color tag and trace back to the particle, then determine the probability) It's done above(IndiceForCol1)
// to define the probability, parameters for distance should exist. possibly, that could be the par() or string number the last
std::vector<vector<double>> MesonrecoMatrix1; // vector of double (prob) , ,,,// in a big scheme, these two vector form a Matrix
std::vector<double> MesonrecoMatrix2; // this is 1d vector of (probabilty)

//test statement
//double test = 0.125;
////std::cout <<endl<<"test number is : "<<test<<endl;
//first, Form Meson recombination Matrix
for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
        double Mesonfactor;
        int distance; // distance between the color tags in string(not c++ function)
        int tag1 = IndiceForColFin.at(irow);
        int tag2 = IndiceForColFin.at(icol);
        std::vector<int>::iterator it1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
        std::vector<int>::iterator it2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag2);// set up for finding the location of co tag in vector that will be used to find distance
        int pos1 = std::distance(IndiceForColFin.begin(), it1); //location of col tag of irow
        int pos2 = std::distance(IndiceForColFin.begin(), it2); //location of col tag of icol
//        std::cout <<endl<<"position of "<<IndiceForColFin.at(irow)<<" is " << pos1 <<endl;
//        std::cout <<endl<<"position of "<<IndiceForColFin.at(icol)<<" is " << pos2 <<endl;

        distance = abs(pos1 - pos2);
 //       std::cout <<endl<<" At ( "<<irow<<","<<icol<<" ) , distance is "<<distance<<endl;

        if(distance == 0) {
        Mesonfactor = 1;
 //       std::cout <<endl<<"Temp Meson Factor is "<<Mesonfactor<<endl;
        MesonrecoMatrix2.push_back(Mesonfactor);
        };
        if(distance > 0) {
        Mesonfactor = 0.111;
//        std::cout <<endl<<"Temp Meson Factor is "<<Mesonfactor<<endl;
        MesonrecoMatrix2.push_back(Mesonfactor);
        };
    }
   MesonrecoMatrix1.push_back(MesonrecoMatrix2);
   MesonrecoMatrix2.clear();
}

// now correct the component in MesonMatrix1 by searching through HH_showerptns{checking gluon}
for(int icheck1 = 0; icheck1 < HH_showerptns.num(); icheck1++ ){
  if(HH_showerptns[icheck1].col() != 0 && HH_showerptns[icheck1].acol() != 0 ){//finding gluon means that finding partons that can't form color singlet, but octet
    int tag1 = HH_showerptns[icheck1].col();
    int tag2 = HH_showerptns[icheck1].acol();
    std::vector<int>::iterator it1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
    std::vector<int>::iterator it2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag2);// set up for finding the location of co tag in vector that will be used to find distance
    int pos1 = std::distance(IndiceForColFin.begin(), it1); //location of col tag of irow
    int pos2 = std::distance(IndiceForColFin.begin(), it2); //location of col tag of icol
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
//MesonrecoMatrix is Formed so far, Now make that of Baryon!!

std::vector<vector<vector<double>>> BaryonrecoMatrix1; // vector of double (prob) , ,,,// in a big scheme, these three items form a Matrix of vector
std::vector<vector<double>> BaryonrecoMatrix2; // this is vector of 2d vector below (probabilty, required color charge for color neutrality), ( , ) , ( , ) ...
std::vector<double> BaryonrecoMatrix3; // this is 2d vector of (probabilty, required color charge for color neutrality)

for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
        double Baryonfactor;
        BaryonrecoMatrix3.push_back((double)1/27);//1/27: probabiity to get singlet state
        BaryonrecoMatrix3.push_back(0); // now the probability factor and color(not specified yet, it'll be done with baryon formation) for neutrality are saved.

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
*/


//std::cout <<endl<<endl;



// TODO: GO line near about 1390 to Check Correction derived by this


	//looping over all the quarks that we have...
	//'q1' loops over all quarks in the shower
	//'q2' loops over all quarks in the event
	//'q3' loops over all quarks in the event, starting from 'q2' and ending at the last quark
	//when 'q2' is at the last quark, we will not consider quark 'q3' - can only make a meson at that point...

	parton_collection considering; int element[3];

	for(int q1=0;q1<showerquarks.num();++q1){
		//accessing first considered quark
		//set q1 variables here

		if(sh_recofactor < 0.0000000001){continue;} //turning off recombination in a computationally friendly way...

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
			else{/*is thermal quark*/ element[1] = perm2[q2] + 1;}

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
				//else if(showerquarks[element[1]].string_id() != showerquarks[element[0]].string_id()){continue;}
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
				//turning off recombination for thermal partons in a computationally friendly way...
				else if(th_recofactor < 0.0000000001){continue;}

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
				else{/*is thermal quark*/ element[2] = perm2[q3] + 1;}

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

//TODO: Here we need to establish the probability about baryon formation. first, need to find 3sets of col tag pair are included junction info at the same time. if all 3 does, prob is 1/27, if 1 does, find other component's tag with element[2](3rd one)'s in Meson Matrix!
//second, scan all junction lists
//before that, we need to evaluate temporary junction, because some junctions could be eliminated, So need to make tempjunction list. then we need to proocedure to erase the factor in the vectors!
//Now evaluating the Tempjunctions element first, before this, declare bool variables for the next proocedure
bool element1 = false;
bool element2 = false;
bool element3 = false;
int juncnum1 = 999999999;
int juncnum2 = 999999999;
int juncnum3 = 999999999;

int standard = IndiceForColFin.size(); // to apply the "if" statement, set the condition for later syntax
int tagformatrix = 999999999; // if all three tags are in same junction, we need to set large number not to confuse, this will be fitered by later syntax
int loc1 = standard; //location in the MesonrecoMatrix
std::vector<int>::iterator I1;
std::vector<int>::iterator I2;
int loc2 = standard;

if(considering[0].id()*considering[1].id()*considering[2].id() > 0){
  // baryon case, and this will check whether two particles are in same junction and the other is not.
  // and evaluate(by mesonrecomatrix) the other particle's color tag with the the color tag of 3rd particle collected to be baryon.
    for(int ijunc=0; ijunc < Tempjunctions.size(); ijunc++){
      if((considering[0].col() == Tempjunctions.at(ijunc).at(1).at(1)) ||
         (considering[0].col() == Tempjunctions.at(ijunc).at(2).at(1)) ||
         (considering[0].col() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
          element1 = true;
          juncnum1 = ijunc;
          }
     if((considering[1].col() == Tempjunctions.at(ijunc).at(1).at(1)) ||
          (considering[1].col() == Tempjunctions.at(ijunc).at(2).at(1)) ||
          (considering[1].col() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
          element2 = true;
          juncnum2 = ijunc;
          }
     if((considering[2].col() == Tempjunctions.at(ijunc).at(1).at(1)) ||
          (considering[2].col() == Tempjunctions.at(ijunc).at(2).at(1)) ||
          (considering[2].col() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
          element2 = true;
          juncnum2 = ijunc;
          }
    }
    if((juncnum1 == juncnum2) && (juncnum1 != juncnum3) && juncnum1 != 999999999 && considering[2].col() != 0){
       tagformatrix = Tempjunctions.at(juncnum1).at(3).at(1);
       I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[2].col());
       loc2 = std::distance(IndiceForColFin.begin(), I2);
       //std::cout <<endl<<"chosen color tag1 is "<<HH_showerptns[showerquarks[element[2]].par()].col();
       //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
    };
    if((juncnum2 == juncnum3) && (juncnum2 != juncnum1) && juncnum2 != 999999999 && considering[0].col() != 0){
      tagformatrix = Tempjunctions.at(juncnum2).at(1).at(1);
      I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[0].col());
      loc2 = std::distance(IndiceForColFin.begin(), I2);
      //std::cout <<endl<<"chosen color tag2 is "<<HH_showerptns[showerquarks[element[0]].par()].col();
      //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
    };
    if((juncnum1 == juncnum3) && (juncnum1 != juncnum2) && juncnum3 != 999999999 && considering[1].col() != 0){
      tagformatrix = Tempjunctions.at(juncnum3).at(2).at(1);
      I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[1].col());
      loc2 = std::distance(IndiceForColFin.begin(), I2);
      //std::cout <<endl<<"chosen color tag3 is "<<HH_showerptns[showerquarks[element[1]].par()].col();
      //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
    };

    I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tagformatrix);
    loc1 = std::distance(IndiceForColFin.begin(), I1);
    //std::cout <<endl<<"chosen orginal color tag is "<<tagformatrix;
    //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc1<<endl;// now we find the locations of the color tags.

//std::cout <<endl<<"loc2 is given as "<<loc2<<endl;
// now, search the mesonrecomatrix
    if (tagformatrix < 999999999 && loc1 < standard && loc2 < standard) { // Not all the tags are in same junction so loc1,2 has always same or smaller value than the size of INdiceForColFin vector
      if(MesonrecoMatrix1.at(loc1).at(loc2) == 1) {
      recofactor3 = 1; // one specific tag, which is located in different location from other two strings from a junction, is neighbored with the tag in the not used tag in junction
      }
      else{recofactor3 = (double) 0;} // general case
    }
    else{recofactor3 = (double) 1/27;} // general case


 }
else if(considering[0].id()*considering[1].id()*considering[2].id() < 0){

  for(int ijunc=0; ijunc < Tempjunctions.size(); ijunc++){
    if((considering[0].acol() == Tempjunctions.at(ijunc).at(1).at(1)) ||
       (considering[0].acol() == Tempjunctions.at(ijunc).at(2).at(1)) ||
       (considering[0].acol() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
        element1 = true;
        juncnum1 = ijunc;
        }
   if((considering[1].acol() == Tempjunctions.at(ijunc).at(1).at(1)) ||
        (considering[1].acol() == Tempjunctions.at(ijunc).at(2).at(1)) ||
        (considering[1].acol() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
        element2 = true;
        juncnum2 = ijunc;
        }
   if((considering[2].acol() == Tempjunctions.at(ijunc).at(1).at(1)) ||
        (considering[2].acol() == Tempjunctions.at(ijunc).at(2).at(1)) ||
        (considering[2].acol() == Tempjunctions.at(ijunc).at(3).at(1)) ) {
        element2 = true;
        juncnum2 = ijunc;
        }
  }
  if((juncnum1 == juncnum2) && (juncnum1 != juncnum3) && juncnum1 != 999999999){
     tagformatrix = Tempjunctions.at(juncnum1).at(3).at(1);
     I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[2].acol());
     loc2 = std::distance(IndiceForColFin.begin(), I2);
     //std::cout <<endl<<"chosen color tag1 is "<<HH_showerptns[showerquarks[element[2]].par()].col();
     //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
  };
  if((juncnum2 == juncnum3) && (juncnum2 != juncnum1) && juncnum2 != 999999999){
    tagformatrix = Tempjunctions.at(juncnum2).at(1).at(1);
    I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[0].acol());
    loc2 = std::distance(IndiceForColFin.begin(), I2);
    //std::cout <<endl<<"chosen color tag2 is "<<HH_showerptns[showerquarks[element[0]].par()].col();
    //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
  };
  if((juncnum1 == juncnum3) && (juncnum1 != juncnum2) && juncnum3 != 999999999){
    tagformatrix = Tempjunctions.at(juncnum3).at(2).at(1);
    I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), considering[1].acol());
    loc2 = std::distance(IndiceForColFin.begin(), I2);
    //std::cout <<endl<<"chosen color tag3 is "<<HH_showerptns[showerquarks[element[1]].par()].col();
    //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc2<<endl;
  };

  I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tagformatrix);
  loc1 = std::distance(IndiceForColFin.begin(), I1);
  //std::cout <<endl<<"chosen orginal color tag is "<<tagformatrix;
  //std::cout <<endl<<"and corresponding indice in the matrix is  "<<loc1<<endl;// now we find the locations of the color tags.

//std::cout <<endl<<"loc2 is given as "<<loc2<<endl;
// now, search the mesonrecomatrix
  if (tagformatrix < 999999999 && loc1 < standard && loc2 < standard) { // Not all the tags are in same junction so loc1,2 has always same or smaller value than the size of INdiceForColFin vector
    if(MesonrecoMatrix1.at(loc1).at(loc2) == 1) {
    recofactor3 = 1; // one specific tag, which is located in different location from other two strings from a junction, is neighbored with the tag in the not used tag in junction
    }
    else{recofactor3 = (double) 0;} // general case
  }
  else{recofactor3 = (double) 1/27;} // general case

 }


					//recofactor3 = sh_recofactor*recofactor2;
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
					//turning off recombination for thermal partons in a computationally friendly way...
					else if(th_recofactor < 0.0000000001){continue;}

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
/*below is for convenience of working-->12/30/2019, should be deleted at last.
if(considering[0].id() > 0 && considering[1].id() < 0){//case of first parton is q and second is q-bar
                    std::cout <<endl<<"chosen partons are "<< int(considering[0].id()) << " and " << int(considering[1].id()) << endl <<endl;
                    std::cout <<endl<<"and their color tag is "<< int(HH_showerptns[showerquarks[element[0]].par()].col()) << " and " << int(HH_showerptns[showerquarks[element[1]].par()].acol()) << endl <<endl;
                    std::cout <<"color correction implemented as Follows: "<< HH_showerptns[showerquarks[element[0]].par()].col() <<" = " << HH_showerptns[showerquarks[element[1]].par()].acol()<<endl;

//before changing the tags, revise the MesonrecoMatrix elements!
                    int coltag1 = HH_showerptns[showerquarks[element[0]].par()].col();
                    int coltag2 = HH_showerptns[showerquarks[element[1]].par()].acol();
                    std::vector<int>::iterator I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                    std::vector<int>::iterator I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                    int loc1 = std::distance(IndiceForColFin.begin(), I1);
                    int loc2 = std::distance(IndiceForColFin.begin(), I2); //set up for find matrix indices corresponding to the color tags

                    MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc1) = 0; //Matrix revised

                    std::cout <<endl<<"Revised Matrix is same as below"<<endl;
                    for(int irow=0; irow < IndiceForColFin.size(); irow++){
                    for(int icol=0; icol < IndiceForColFin.size(); icol++){
                    std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                    }
                    std::cout <<endl<<endl;
                    }


                    HH_showerptns[showerquarks[element[0]].par()].col(HH_showerptns[showerquarks[element[1]].par()].acol());//now color tags from both partons are same

                    }
*/
//declare the vector needed to form fake parton to meet the need of PYTHIA {color neutrality}


//TODO:Starting point of the string repair after Baryon Formation!
		 	if(WigB[1]*recofactor3 >= rndbaryon){
				//if(rndbaryon >= 0.99){
/*This Variables are in header file
                    vector<vector<vector<int>>> Tempjunctions; // vector of all tempjunctions
                    vector<vector<int>> JunctionInfo; // vector of one junction's color tag(3) and particle info
                    vector<int> IdColInfo1;
                    vector<int> IdColInfo2;
                    vector<int> IdColInfo3;
                    vector<int> IdColInfo4;
                    int candidateNum = 1;
*/

/*
          for(int ibary = 0; ibary < 3; ibary++){
              if(HH_showerptns[0].col() > 0 && HH_showerptns[0].acol() == 0){ //check initiating particle{in this case, it's quark, so endpoint fakeparton should be anti quark}
                if(considering[ibary].col() == HH_showerptns[HH_showerptns.num()-1].col()){
                fakepartoninfo.push_back(-1);
                fakepartoninfo.push_back(HH_showerptns[HH_showerptns.num()-1].col());
                //dignostic measure
                std::cout <<endl<<"fake parton "<<" ( "<<fakepartoninfo.at(0)<<" , "<<fakepartoninfo.at(1)<<" ) added!"<<endl;
                break;
                }
              }
              else if(HH_showerptns[0].col() == 0 && HH_showerptns[0].acol() > 0){ //check initiating particle{in this case, it's anti-quark, so endpoint fakeparton should be quark}
              if(considering[ibary].acol() == HH_showerptns[HH_showerptns.num()-1].acol()){
                fakepartoninfo.push_back(1);
                fakepartoninfo.push_back(HH_showerptns[HH_showerptns.num()-1].acol());
                //dignostic measure
                std::cout <<endl<<"fake parton "<<" ( "<<fakepartoninfo.at(0)<<" , "<<fakepartoninfo.at(1)<<" ) added!"<<endl;
                break;
               }
              } // this data will be used for appending fake partons {at the last part of recomb function, pid, pstat col, acol needed, }


            }
            */


			        if(considering[0].id() * considering[1].id() * considering[2].id() > 0){
                //TODO:If baryon is formed from Tempjunction, corresponding temp junction should be eliminated.
                                    //std::cout <<endl<<"Baryon is made from "<<" ( "<<considering[0].col()<<" , "<<considering[1].col()<<" , "<<considering[2].col()<<" ) "<<endl;
                                    for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                                        if( (Tempjunctions.at(ijunc).at(1).at(1) == considering[0].col()) ||
                                        (Tempjunctions.at(ijunc).at(1).at(1) == considering[1].col()) ||
                                        (Tempjunctions.at(ijunc).at(1).at(1) == considering[2].col())){
                                          Tempjunctions.at(ijunc).at(1).pop_back();
                                          Tempjunctions.at(ijunc).at(1).push_back(999);
                                          Tempjunctions.at(ijunc).at(2).pop_back();
                                          Tempjunctions.at(ijunc).at(2).push_back(999);
                                          Tempjunctions.at(ijunc).at(3).pop_back();
                                          Tempjunctions.at(ijunc).at(3).push_back(999);// replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                                        }
                                      }
                                      for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                                        if( (Tempjunctions.at(ijunc).at(2).at(1) == considering[0].col()) ||
                                        (Tempjunctions.at(ijunc).at(2).at(1) == considering[1].col()) ||
                                        (Tempjunctions.at(ijunc).at(2).at(1) == considering[2].col()) ){
                                          Tempjunctions.at(ijunc).at(1).pop_back();
                                          Tempjunctions.at(ijunc).at(1).push_back(888);
                                          Tempjunctions.at(ijunc).at(2).pop_back();
                                          Tempjunctions.at(ijunc).at(2).push_back(888);
                                          Tempjunctions.at(ijunc).at(3).pop_back();
                                          Tempjunctions.at(ijunc).at(3).push_back(888);// replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                                        }
                                      }
                                      for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                                        if( (Tempjunctions.at(ijunc).at(3).at(1) == considering[0].col()) ||
                                        (Tempjunctions.at(ijunc).at(3).at(1) == considering[1].col() ||
                                        (Tempjunctions.at(ijunc).at(3).at(1) == considering[2].col())) ){
                                          Tempjunctions.at(ijunc).at(1).pop_back();
                                          Tempjunctions.at(ijunc).at(1).push_back(777);
                                          Tempjunctions.at(ijunc).at(2).pop_back();
                                          Tempjunctions.at(ijunc).at(2).push_back(777);
                                          Tempjunctions.at(ijunc).at(3).pop_back();
                                          Tempjunctions.at(ijunc).at(3).push_back(777);// replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                                      }
                                  } // now corresponding Tempjunction vector is cleared!
                // baryon is to be formed, so temporary anti junction is defined
                    if(considering[0].col() > 0 && considering[1].col() > 0 && considering[2].col() > 0) {
                    IdColInfo1.push_back(-1); // antijunction tag(-1)
                    IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st
                    IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                    IdColInfo2.push_back(considering[0].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                    IdColInfo3.push_back(-1);
                    IdColInfo3.push_back(considering[1].col());
                    IdColInfo4.push_back(-1);
                    IdColInfo4.push_back(considering[2].col());

                    JunctionInfo.push_back(IdColInfo1);
                    JunctionInfo.push_back(IdColInfo2);
                    JunctionInfo.push_back(IdColInfo3);
                    JunctionInfo.push_back(IdColInfo4);

                    Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                    IdColInfo1.clear();
                    IdColInfo2.clear();
                    IdColInfo3.clear();
                    IdColInfo4.clear();
                    JunctionInfo.clear();



//Since Baryon is formed, color tag for neutrality should be added, casting int into double is needed!!(ex (double) intvalue
                    int coltag1 = considering[0].col();
                    int coltag2 = considering[1].col();
                    int coltag3 = considering[2].col();
                    if(coltag1 > 0 && coltag2 > 0 && coltag3 && coltag1 < limit && coltag2 < limit && coltag3 < limit ){
                    double tag1 = (double)coltag1;  // they are casted to be inserted into the matrix(since it's vector of double
                    double tag2 = (double)coltag2;
                    double tag3 = (double)coltag3;
                    std::vector<int>::iterator I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                    std::vector<int>::iterator I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                    std::vector<int>::iterator I3 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag3);
                    int loc1 = std::distance(IndiceForColFin.begin(), I1);
                    int loc2 = std::distance(IndiceForColFin.begin(), I2);
                    int loc3 = std::distance(IndiceForColFin.begin(), I3); //set up for find matrix indices corresponding to the color tags (we just found the corresponding indice in BaryonrecoMatrix1 with col tags )

                    BaryonrecoMatrix1.at(loc1).at(loc2).at(1) = tag3;
                    BaryonrecoMatrix1.at(loc2).at(loc1).at(1) = tag3;
                    BaryonrecoMatrix1.at(loc2).at(loc3).at(1) = tag1;
                    BaryonrecoMatrix1.at(loc3).at(loc2).at(1) = tag1;
                    BaryonrecoMatrix1.at(loc1).at(loc3).at(1) = tag2;
                    BaryonrecoMatrix1.at(loc3).at(loc1).at(1) = tag2; // now the color tag info for color neutrality is saved in Matrix, also we need to consider the impact from this to Meson Formation

//MesonrecoMatrix1 is modified by below
                    MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc1) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc3) = 0;
                    MesonrecoMatrix1.at(loc3).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc1).at(loc3) = 0;
                    MesonrecoMatrix1.at(loc3).at(loc1) = 0; // since three color tags of baryon are different from each other, so that these tags can't form meson with each other.
                  }
/*
//dignostic measure
                   std::cout <<endl<<"Meson reco Matrix revised by Baryon formation is same as below"<<endl;
                   for(int irow=0; irow < IndiceForColFin.size(); irow++){
                       for(int icol=0; icol < IndiceForColFin.size(); icol++){
                       std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                       }
                   std::cout <<endl<<endl;
                   }
                   */
//dignostic measure
/*
                   std::cout <<endl<<"Revised Baryon reco Matrix is same as below"<<endl;
                   for(int irow=0; irow < IndiceForColFin.size(); irow++){
                       for(int icol=0; icol < IndiceForColFin.size(); icol++){
                       std::cout <<"  ( "<<BaryonrecoMatrix1.at(irow).at(icol).at(1)<<" )  ";
                       }
                   std::cout <<endl<<endl;
                   }
*/

                       }
                    if(considering[0].col() > 0 && considering[1].col() == 0 && considering[2].col() > 0){ // MAT + Therm/LBT + MAT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      maxtag++;
                      int loc = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                      if(loc == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc > 0){showerquarks[loc - 1].acol(maxtag); }
                      else if(loc < 0){HH_thermal[-loc - 1].acol(maxtag); }

                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                      IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo2.push_back(considering[0].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(maxtag);
                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(considering[2].col());

                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    }

                    if(considering[0].col() == 0 && considering[1].col() > 0 && considering[2].col() > 0){ // LBT + MAT + MAT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      maxtag++;
                      int loc = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                      if(loc == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc > 0){showerquarks[loc - 1].acol(maxtag); }
                      else if(loc < 0){HH_thermal[-loc - 1].acol(maxtag); }
                      //TODO: purpose of statement above is to give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                      IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo2.push_back(maxtag);// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(considering[1].col());
                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(considering[2].col());

                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    }

                    if(considering[0].col() > 0 && considering[1].col() > 0 && considering[2].col() == 0){ // MAT + MAT + Therm/LBT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      maxtag++;
                      int loc = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                      if(loc == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc > 0){showerquarks[loc - 1].acol(maxtag); }
                      else if(loc < 0){HH_thermal[-loc - 1].acol(maxtag); }

                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                      IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo2.push_back(considering[0].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(considering[1].col());
                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(maxtag);

                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    } // so far, we are finished with the list with one zero component.


                    if(considering[0].col() > 0 && considering[1].col() == 0 && considering[2].col() == 0){ // MAT + Therm/LBT + Therm/LBT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st
                      IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo2.push_back(considering[0].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th

                      maxtag++;
                      int loc1 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                      if(loc1 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc1 > 0){showerquarks[loc1 - 1].acol(maxtag); }
                      else if(loc1 < 0){HH_thermal[-loc1 - 1].acol(maxtag); }
                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(maxtag);

                      maxtag++;
                      int loc2 = findcloserepl(considering[2], perm2[q3], true, true, showerquarks, HH_thermal);
                      if(loc2 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc2 > 0){showerquarks[loc2 - 1].acol(maxtag); }
                      else if(loc2 < 0){HH_thermal[-loc2 - 1].acol(maxtag); }

                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(maxtag);
                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    }

                    if(considering[0].col() == 0 && considering[1].col() == 0 && considering[2].col() > 0){ //  Therm/LBT + Therm/LBT + MAT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                      if(loc1 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc1 > 0){showerquarks[loc1 - 1].acol(++maxtag); }
                      else if(loc1 < 0){HH_thermal[-loc1 - 1].acol(++maxtag); }
                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(maxtag);



                      int loc2 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                      if(loc2 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc2 > 0){showerquarks[loc2 - 1].acol(++maxtag); }

                      else if(loc2 < 0){HH_thermal[-loc2 - 1].acol(++maxtag); }

                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(maxtag);
                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.



                      IdColInfo2.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo2.push_back(considering[2].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th


                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    }

                    if(considering[0].col() == 0 && considering[1].col() > 0 && considering[2].col() == 0){ //  Therm/LBT + Therm/LBT + MAT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      maxtag++;
                      int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                      if(loc1 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc1 > 0){showerquarks[loc1 - 1].acol(maxtag); }
                      else if(loc1 < 0){HH_thermal[-loc1 - 1].acol(maxtag); }
                      IdColInfo2.push_back(-1);
                      IdColInfo2.push_back(maxtag);


                      IdColInfo3.push_back(-1); // means anticolor tag(negative color charge)
                      IdColInfo3.push_back(considering[1].col());// {-1, anticolor tag} will be at 2nd, 3rd, 4th



                      int loc2 = findcloserepl(considering[2], perm2[q3], true, true, showerquarks, HH_thermal);
                      if(loc2 == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[2]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc2 > 0){showerquarks[loc2 - 1].acol(++maxtag); }
                      else if(loc2 < 0){HH_thermal[-loc2 - 1].acol(++maxtag); }
                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(maxtag);
                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.





                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    } // so far, we are finished with the list with two zeros

/*
                    if(considering[0].col() == 0 && considering[1].col() == 0 && considering[2].col() == 0){ //  Therm/LBT + Therm/LBT + MAT case
                      IdColInfo1.push_back(-1); // antijunction tag(-1)
                      IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                      maxtag++;
                      int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                      if (loc1 > 0){showerquarks[loc1 - 1].acol(maxtag); }
                      else if(loc1 < 0){HH_thermal[-loc1 - 1].acol(maxtag); }
                      IdColInfo2.push_back(-1);
                      IdColInfo2.push_back(maxtag);


                      maxtag++;
                      int loc2 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                      if(loc2 = 999999999){
                         std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(maxtag);
                         //Extraparton collection.add{fakep};
                         //somewhere, we need to make for loop to toss all partons to remnants list.


                      }

                      if (loc2 > 0){showerquarks[loc2 - 1].acol(maxtag); }
                      else if(loc2 < 0){HH_thermal[-loc2 - 1].acol(maxtag); }


                      IdColInfo3.push_back(-1);
                      IdColInfo3.push_back(maxtag);


                      int loc3 = findcloserepl(considering[2], perm2[q3], true, true, showerquarks, HH_thermal);
                      if (loc3 > 0){showerquarks[loc3 - 1].acol(++maxtag); }
                      else if(loc3 < 0){HH_thermal[-loc3 - 1].acol(++maxtag); }
                      IdColInfo4.push_back(-1);
                      IdColInfo4.push_back(maxtag);
                      //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.





                      JunctionInfo.push_back(IdColInfo1);
                      JunctionInfo.push_back(IdColInfo2);
                      JunctionInfo.push_back(IdColInfo3);
                      JunctionInfo.push_back(IdColInfo4);

                      Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                      IdColInfo1.clear();
                      IdColInfo2.clear();
                      IdColInfo3.clear();
                      IdColInfo4.clear();
                      JunctionInfo.clear();

                    } // so far, we are finished with the list with three zeros
                    */




                    }

                    else if(considering[0].id() * considering[1].id() * considering[2].id() < 0){ //anti baryon would be made so temporary junction is defined
                      int col1 = HH_showerptns[showerquarks[element[0]].par()].acol();
                      int col2 = 0;
                      int col3 = 0;
                    if(perm2[q2] > 0){
                      col2 = HH_showerptns[showerquarks[element[1]].par()].acol();
                    }
                    if(perm2[q3] > 0){
                      col3 = HH_showerptns[showerquarks[element[2]].par()].acol();
                    }

                    //std::cout <<endl<<"Anti Baryon is made from "<<" ( "<<considering[0].acol()<<" , "<<considering[1].acol()<<" , "<<considering[2].acol()<<" ) "<<endl;
                  //  std::cout <<endl<<"Anti Baryon is made from "<<" ( "<<HH_showerptns[showerquarks[element[0]].par()].acol()<<" , "<<HH_showerptns[showerquarks[element[1]].par()].acol()<<" , "<<HH_showerptns[showerquarks[element[2]].par()].acol()<<" ) "<<endl;


                    for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                        if( (Tempjunctions.at(ijunc).at(1).at(1) == considering[0].acol()) ||
                        (Tempjunctions.at(ijunc).at(1).at(1) == considering[1].acol()) ||
                        (Tempjunctions.at(ijunc).at(1).at(1) == considering[2].acol())){
                         Tempjunctions.at(ijunc).at(1).pop_back();
                         Tempjunctions.at(ijunc).at(1).push_back(999);
                         Tempjunctions.at(ijunc).at(2).pop_back();
                         Tempjunctions.at(ijunc).at(2).push_back(999);
                         Tempjunctions.at(ijunc).at(3).pop_back();
                         Tempjunctions.at(ijunc).at(3).push_back(999); // replace all the col tag in junction candidate to preserve the size of the vector, or the code would be broken
                        }
                      }
                      for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                        if( (Tempjunctions.at(ijunc).at(2).at(1) == considering[0].acol()) ||
                        (Tempjunctions.at(ijunc).at(2).at(1) == considering[1].acol()) ||
                        (Tempjunctions.at(ijunc).at(2).at(1) == considering[2].acol()) ){
                          Tempjunctions.at(ijunc).at(1).pop_back();
                          Tempjunctions.at(ijunc).at(1).push_back(999);
                          Tempjunctions.at(ijunc).at(2).pop_back();
                          Tempjunctions.at(ijunc).at(2).push_back(999);
                          Tempjunctions.at(ijunc).at(3).pop_back();
                          Tempjunctions.at(ijunc).at(3).push_back(999);// replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                        }
                      }
                      for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                        if( (Tempjunctions.at(ijunc).at(3).at(1) == considering[0].acol()) ||
                        (Tempjunctions.at(ijunc).at(3).at(1) == considering[1].acol() ||
                        (Tempjunctions.at(ijunc).at(3).at(1) == considering[2].acol())) ){
                          Tempjunctions.at(ijunc).at(1).pop_back();
                          Tempjunctions.at(ijunc).at(1).push_back(999);
                          Tempjunctions.at(ijunc).at(2).pop_back();
                          Tempjunctions.at(ijunc).at(2).push_back(999);
                          Tempjunctions.at(ijunc).at(3).pop_back();
                          Tempjunctions.at(ijunc).at(3).push_back(999);// replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                      }
                  } // now corresponding Tempjunction vector is cleared!
                  if(considering[0].acol() > 0 && considering[1].acol() > 0 && considering[2].acol() > 0){
                  IdColInfo1.push_back(1); // junction tag(-1)
                  IdColInfo1.push_back(0); // zero(just room for the other usage)   : {1, 0} at 1st
                  IdColInfo2.push_back(1); // means color tag(positive color charge)   : { 1, color tag } at 2nd, 3rd, 4th in the vector
                  IdColInfo2.push_back(considering[0].acol());
                  IdColInfo3.push_back(1);
                  IdColInfo3.push_back(considering[1].acol());
                  IdColInfo4.push_back(1);
                  IdColInfo4.push_back(considering[2].acol());

                  JunctionInfo.push_back(IdColInfo1);
                  JunctionInfo.push_back(IdColInfo2);
                  JunctionInfo.push_back(IdColInfo3);
                  JunctionInfo.push_back(IdColInfo4);

                  Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                  IdColInfo1.clear();
                  IdColInfo2.clear();
                  IdColInfo3.clear();
                  IdColInfo4.clear();
                  JunctionInfo.clear();

//Since Baryon is formed, color tag for neutrality should be added, casting int into double is needed!!(ex (double) intvalue
                    int coltag1 = considering[0].acol();
                    int coltag2 = considering[1].acol();
                    int coltag3 = considering[2].acol();
                    if(coltag1 > 0 && coltag2 > 0 && coltag3 > 0 && coltag1 < limit  && coltag2 < limit  && coltag3 < limit ){
                    double tag1 = (double)coltag1;  // they are casted to be inserted into the matrix(since it's vector of double
                    double tag2 = (double)coltag2;
                    double tag3 = (double)coltag3;
                    std::vector<int>::iterator I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                    std::vector<int>::iterator I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                    std::vector<int>::iterator I3 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag3);
                    int loc1 = std::distance(IndiceForColFin.begin(), I1);
                    int loc2 = std::distance(IndiceForColFin.begin(), I2);
                    int loc3 = std::distance(IndiceForColFin.begin(), I3); //set up for find matrix indices corresponding to the color tags (we just found the corresponding indice in BaryonrecoMatrix1 with col tags )

                    BaryonrecoMatrix1.at(loc1).at(loc2).at(1) = tag3;
                    BaryonrecoMatrix1.at(loc2).at(loc1).at(1) = tag3;
                    BaryonrecoMatrix1.at(loc2).at(loc3).at(1) = tag1;
                    BaryonrecoMatrix1.at(loc3).at(loc2).at(1) = tag1;
                    BaryonrecoMatrix1.at(loc1).at(loc3).at(1) = tag2;
                    BaryonrecoMatrix1.at(loc3).at(loc1).at(1) = tag2; // now the color tag info for color neutrality is saved in Matrix, also we need to consider the impact from this to Meson Formation

//MesonrecoMatrix1 is modified by below
                    MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc1) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc3) = 0;
                    MesonrecoMatrix1.at(loc3).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc1).at(loc3) = 0;
                    MesonrecoMatrix1.at(loc3).at(loc1) = 0; // since three color tags of baryon are different from each other, so that these tags can't form meson with each other.

                   }
                 }
                 if(considering[0].acol() > 0 && considering[1].acol() == 0 && considering[2].acol() > 0){ // MAT + Therm/LBT + MAT case
                   IdColInfo1.push_back(1); // junction tag(-1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                   maxtag++;
                   int loc = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                   if(loc == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc > 0){showerquarks[loc - 1].col(maxtag); }
                   else if(loc < 0){HH_thermal[-loc - 1].col(maxtag); }
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                   IdColInfo2.push_back(1); // color tag(negative color charge)
                   IdColInfo2.push_back(considering[0].acol());// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(maxtag);
                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(considering[2].acol());

                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 }

                 if(considering[0].acol() == 0 && considering[1].acol() > 0 && considering[2].acol() > 0){ // LBT + MAT + MAT case
                   IdColInfo1.push_back(1); // junction tag(1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {1, 0} will be in 1st

                   maxtag++;
                   int loc = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                   if(loc == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc > 0){showerquarks[loc - 1].col(maxtag); }
                   else if(loc < 0){HH_thermal[-loc - 1].col(maxtag); }
                   //TODO: purpose of statement above is to give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                   IdColInfo2.push_back(1); // means anticolor tag(negative color charge)
                   IdColInfo2.push_back(maxtag);// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(considering[1].acol());
                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(considering[2].acol());

                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 }

                 if(considering[0].col() > 0 && considering[1].col() > 0 && considering[2].col() == 0){ // MAT + MAT + Therm/LBT case
                   IdColInfo1.push_back(1); // junction tag(1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                   maxtag++;
                   int loc = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                   if(loc == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc > 0){showerquarks[loc - 1].col(maxtag); }
                   else if(loc < 0){HH_thermal[-loc - 1].col(maxtag); }
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                   IdColInfo2.push_back(1); // means anticolor tag(negative color charge)
                   IdColInfo2.push_back(considering[0].acol());// {-1, anticolor tag} will be at 2nd, 3rd, 4th
                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(considering[1].acol());
                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(maxtag);

                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 } // so far, we are finished with the list with one zero component.


                 if(considering[0].acol() > 0 && considering[1].acol() == 0 && considering[2].acol() == 0){ // MAT + Therm/LBT + Therm/LBT case
                   IdColInfo1.push_back(1); // antijunction tag(-1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st
                   IdColInfo2.push_back(1); // means anticolor tag(negative color charge)
                   IdColInfo2.push_back(considering[0].acol());// {-1, anticolor tag} will be at 2nd, 3rd, 4th

                   maxtag++;
                   int loc1 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                   if(loc1 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc1 > 0){showerquarks[loc1 - 1].col(maxtag); }
                   else if(loc1 < 0){HH_thermal[-loc1 - 1].col(maxtag); }
                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(maxtag);

                   maxtag++;
                   int loc2 = findcloserepl(considering[2], perm2[q3], true, true, showerquarks, HH_thermal);
                   if(loc2 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[2]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc2 > 0){showerquarks[loc2 - 1].col(maxtag); }
                   else if(loc2 < 0){HH_thermal[-loc2 - 1].col(maxtag); }

                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(maxtag);
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.


                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 }

                 if(considering[0].acol() == 0 && considering[1].acol() == 0 && considering[2].acol() > 0){ //  Therm/LBT + Therm/LBT + MAT case
                   IdColInfo1.push_back(1); // junction tag(1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                   int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                   if(loc1 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc1 > 0){showerquarks[loc1 - 1].col(++maxtag); }
                   else if(loc1 < 0){HH_thermal[-loc1 - 1].col(++maxtag); }
                   IdColInfo2.push_back(1);
                   IdColInfo2.push_back(maxtag);



                   int loc2 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                   if(loc2 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc2 > 0){showerquarks[loc2 - 1].col(++maxtag); }
                   else if(loc2 < 0){HH_thermal[-loc2 - 1].col(++maxtag); }

                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(maxtag);
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.



                   IdColInfo4.push_back(-1); // means anticolor tag(negative color charge)
                   IdColInfo4.push_back(considering[2].acol());// {-1, anticolor tag} will be at 2nd, 3rd, 4th


                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 }

                 if(considering[0].acol() == 0 && considering[1].acol() > 0 && considering[2].acol() == 0){ //  Therm/LBT + Therm/LBT + MAT case
                   IdColInfo1.push_back(-1); // junction tag(1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {-1, 0} will be in 1st

                   maxtag++;
                   int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                   if(loc1 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc1 > 0){showerquarks[loc1 - 1].col(maxtag); }
                   else if(loc1 < 0){HH_thermal[-loc1 - 1].col(maxtag); }
                   IdColInfo2.push_back(1);
                   IdColInfo2.push_back(maxtag);


                   IdColInfo3.push_back(1); // means anticolor tag(negative color charge)
                   IdColInfo3.push_back(considering[1].acol());// {-1, anticolor tag} will be at 2nd, 3rd, 4th



                   int loc2 = findcloserepl(considering[2], perm2[q2], true, true, showerquarks, HH_thermal);
                   if(loc2 == 999999999){
                      //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                      HHparton fakep = considering[2]; fakep.id(-fakep.id());  fakep.col(maxtag);
                      Extraparton.add(fakep);
                      //somewhere, we need to make for loop to toss all partons to remnants list.
                   }
                   else if (loc2 > 0){showerquarks[loc2 - 1].col(++maxtag); }
                   else if(loc2 < 0){HH_thermal[-loc2 - 1].col(++maxtag); }
                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(maxtag);
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.





                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 } // so far, we are finished with the list with two zeros



                 /*
                 if(considering[0].acol() == 0 && considering[1].acol() == 0 && considering[2].acol() == 0){ //  Therm/LBT + Therm/LBT + MAT case
                   IdColInfo1.push_back(1); // junction tag(1)
                   IdColInfo1.push_back(0); // zero(just room for the other usage) {1, 0} will be in 1st

                   maxtag++;
                   int loc1 = findcloserepl(considering[0], element[0] + 1, true, true, showerquarks, HH_thermal);
                   if (loc1 > 0){showerquarks[loc1 - 1].col(maxtag); }
                   else if(loc1 < 0){HH_thermal[-loc1 - 1].col(maxtag); }
                   IdColInfo2.push_back(1);
                   IdColInfo2.push_back(maxtag);


                   maxtag++;
                   int loc2 = findcloserepl(considering[1], perm2[q2], true, true, showerquarks, HH_thermal);
                   if (loc2 > 0){showerquarks[loc2 - 1].col(maxtag); }
                   else if(loc2 < 0){HH_thermal[-loc2 - 1].col(maxtag); }
                   IdColInfo3.push_back(1);
                   IdColInfo3.push_back(maxtag);


                   int loc3 = findcloserepl(considering[2], perm2[q3], true, true, showerquarks, HH_thermal);
                   if (loc3 > 0){showerquarks[loc3 - 1].col(++maxtag); }
                   else if(loc3 < 0){HH_thermal[-loc3 - 1].col(++maxtag); }
                   IdColInfo4.push_back(1);
                   IdColInfo4.push_back(maxtag);
                   //TODO: give thermal partons "ANTI COLOR TAGS" to form anti junction and conserve baryon number.





                   JunctionInfo.push_back(IdColInfo1);
                   JunctionInfo.push_back(IdColInfo2);
                   JunctionInfo.push_back(IdColInfo3);
                   JunctionInfo.push_back(IdColInfo4);

                   Tempjunctions.push_back(JunctionInfo); // information of tempjunction is saved so et clear subordinate vector for next entry

                   IdColInfo1.clear();
                   IdColInfo2.clear();
                   IdColInfo3.clear();
                   IdColInfo4.clear();
                   JunctionInfo.clear();

                 } // so far, we are finished with the list with three zeros

                 */
//dignostic measure
/*
std::cout <<endl<<"Meson reco Matrix revised by Baryon formation is same as below"<<endl;
for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
    std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
    }
std::cout <<endl<<endl;
}
*/
//dignostic measure
/*
std::cout <<endl<<"Revised Baryon reco Matrix is same as below"<<endl;
for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
    std::cout <<"  ( "<<BaryonrecoMatrix1.at(irow).at(icol).at(1)<<" )  ";
    }
std::cout <<endl<<endl;
}
*/

                    }


//std::cout <<endl<< " I'm Working!! from line 1207 "<<endl;
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



/*below is for reference for the process below(12/23/2019)--shoudl be deleted at last

//dignostic measure
std::cout <<endl<<"Matrix is same as below"<<endl;
for(int irow=0; irow < IndiceForColFin.size(); irow++){
    for(int icol=0; icol < IndiceForColFin.size(); icol++){
    std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
    }
std::cout <<endl<<endl;
}

*/

			else if(considering[0].id()*considering[1].id() < 0){
//TODO: Key point is determing recofactor2, which is same as "squared" Value in MesonColRecomb Vector! Keep it in mind first of all
               if(considering[0].id() > 0 && considering[1].id() < 0){
                       int tag0 = considering[0].col();
                       int tag1 = considering[1].acol();
                       if(tag0 > 0 && tag1 > 0 && tag0 < limit  && tag1 < limit ){
                       std::vector<int>::iterator L1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag0);
                       std::vector<int>::iterator L2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
                       int indexMatrix1 = std::distance(IndiceForColFin.begin(), L1);
                       int indexMatrix2 = std::distance(IndiceForColFin.begin(), L2);

                       recofactor2 = (MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2));
//                       std::cout <<endl<<" the recofactor for meson is "<<(MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2))<<endl;
                       }
                      else{recofactor2 = (double) 1/9;}
                  }


               else if(considering[0].id() > 0 && considering[1].id() < 0){

                       int tag0 = considering[0].col();
                       int tag1 = considering[1].acol();
                       if(tag0 > 0 && tag1 > 0 && tag0 < limit && tag1 < limit){
                       std::vector<int>::iterator L1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag0);
                       std::vector<int>::iterator L2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), tag1);
                       int indexMatrix1 = std::distance(IndiceForColFin.begin(), L1);
                       int indexMatrix2 = std::distance(IndiceForColFin.begin(), L2);

                       recofactor2 = (MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2));
                       }
                       else{recofactor2 = (double) 1/9;}
//                       std::cout <<endl<<" the recofactor for meson is "<<(MesonrecoMatrix1.at(indexMatrix1).at(indexMatrix2))<<endl;
                }
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
//				if(WigM[1]*recofactor2 >= 1){
				if(WigM[1]*recofactor2 >= rndmeson){
//TODO:confirm if the partons not used into Meson get back to the original permutation list
//std::cout <<endl<<"Let's find what was recombined!"<<endl<<endl;


                    if(considering[0].id() > 0 && considering[1].id() < 0){//case of first parton is q and second is q-bar
                    //std::cout <<endl<<"chosen partons are "<< int(considering[0].id()) << " and " << int(considering[1].id()) << endl <<endl;
                    //std::cout <<endl<<"and their color tag is "<< int(considering[0].col()) << " and " << int(considering[1].acol()) << endl <<endl;
                    //std::cout <<"color correction implemented as Follows: "<< considering[0].col() <<" = " << considering[1].acol()<<endl;

//TODO:If Meson is formed from Tempjunction, corresponding temp junction should be eliminated.
                    for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(1).at(1) == considering[0].col()) ||
                           (Tempjunctions.at(ijunc).at(1).at(1) == considering[1].acol()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(999);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(999);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(999);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                        }
                        for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(2).at(1) == considering[0].col()) ||
                           (Tempjunctions.at(ijunc).at(2).at(1) == considering[1].acol()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(888);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(888);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(888);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                        }
                        for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(3).at(1) == considering[0].col()) ||
                           (Tempjunctions.at(ijunc).at(3).at(1) == considering[1].acol()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(777);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(777);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(777);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                    } // now corresponding Tempjunction vector is cleared!

//before changing the tags, revise the MesonrecoMatrix elements!
                    int coltag1 = considering[0].col();
                    int coltag2 = considering[1].acol();
                    if(coltag1 > 0 && coltag2 && coltag1 < limit && coltag2 < limit){
                    std::vector<int>::iterator I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                    std::vector<int>::iterator I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                    int loc1 = std::distance(IndiceForColFin.begin(), I1);
                    int loc2 = std::distance(IndiceForColFin.begin(), I2); //set up for find matrix indices corresponding to the color tags

                    MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc1) = 0; //Matrix revised
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

            //TODO: Need to be changed to reflect thermal parton(04232020)
            //possible case  MAT+LBT, MAT+THERM,
            //treatment 1 : based on distance.
                   if(considering[0].col() > 0 && considering[1].acol() > 0){ //MAT/lbt or therm with color tags + MAT/lbt or therm with color tags
                     if(perm2[q2] > 0){
                      HH_showerptns[showerquarks[element[1]].par()].acol(considering[0].col());//now color tags from both partons are same
                     }
                     else{
                     HH_thermal[-element[1]].acol(considering[0].col());
                     }
                     }


                    else if(considering[0].col() > 0){ // MAT + LBT/THERM
                    int thermsib = findcloserepl(considering[0] , perm1[q1]+1, true, true, HH_showerptns, HH_thermal); // functon to find
                    if(thermsib == 999999999){
                       //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                       HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
                       Extraparton.add(fakep);
                       //somewhere, we need to make for loop to toss all partons to remnants list.
                    }
                    else if(thermsib < 0){
                    HH_thermal[-thermsib-1].acol(considering[0].col());
                    }
                    else{HH_showerptns[thermsib-1].acol(considering[0].col()); }
                    }
                    else if(considering[1].acol() > 0){ //LBT + MAT
                    int thermsib = findcloserepl(considering[1] , perm2[q2], true, true, HH_showerptns, HH_thermal);
                    if(thermsib == 999999999){
                       //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                       HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
                       Extraparton.add(fakep);
                       //somewhere, we need to make for loop to toss all partons to remnants list.
                    }
                    else if(thermsib < 0){
                    HH_thermal[-thermsib-1].col(considering[1].acol());
                    }
                    else{HH_showerptns[thermsib-1].col(considering[1].acol()); }
                    }
                  }

                    else if(considering[0].id() < 0 && considering[1].id() > 0){//case of first parton is q-bar and second is q
                    //std::cout <<endl<<"chosen partons are "<< int(considering[0].id()) << " and " << int(considering[1].id()) << endl <<endl;
                    //std::cout <<endl<<"and their color tag is "<< int(considering[0].acol()) << " and " << int(considering[1].col()) << endl <<endl;
                    //std::cout <<"color correction implemented as Follows: "<< considering[0].acol()<<" = "<<considering[1].col()<<endl;

                    for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(1).at(1) == considering[0].acol()) ||
                           (Tempjunctions.at(ijunc).at(1).at(1) == considering[1].col()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(999);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(999);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(999);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                        }
                        for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(2).at(1) == considering[0].acol()) ||
                           (Tempjunctions.at(ijunc).at(2).at(1) == considering[1].col()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(999);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(999);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(999);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                        }
                        for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
                       if( (Tempjunctions.at(ijunc).at(3).at(1) == considering[0].acol()) ||
                           (Tempjunctions.at(ijunc).at(3).at(1) == considering[1].col()) ){
                             Tempjunctions.at(ijunc).at(1).pop_back();
                             Tempjunctions.at(ijunc).at(1).push_back(999);
                             Tempjunctions.at(ijunc).at(2).pop_back();
                             Tempjunctions.at(ijunc).at(2).push_back(999);
                             Tempjunctions.at(ijunc).at(3).pop_back();
                             Tempjunctions.at(ijunc).at(3).push_back(999);
                             // replace all the col tag in junction candidate with zero to preserve the size of the vector, or the code would be broken
                          }
                    } // now corresponding Tempjunction vector is cleared!

//before changing the tags, revise the MesonrecoMatrix elements!
                    int coltag1 = considering[0].acol();
                    int coltag2 = considering[1].col();
                    if(coltag1 > 0 && coltag2 && coltag1 < limit && coltag2 < limit){
                    std::vector<int>::iterator I1 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag1);
                    std::vector<int>::iterator I2 = std::find(IndiceForColFin.begin(), IndiceForColFin.end(), coltag2);
                    int loc1 = std::distance(IndiceForColFin.begin(), I1);
                    int loc2 = std::distance(IndiceForColFin.begin(), I2); //set up for find matrix indices corresponding to the color tags

                    MesonrecoMatrix1.at(loc1).at(loc2) = 0;
                    MesonrecoMatrix1.at(loc2).at(loc1) = 0; //Matrix revised
                    }
/*
                    std::cout <<endl<<"Revised Matrix is same as below"<<endl;
                    for(int irow=0; irow < IndiceForColFin.size(); irow++){
                    for(int icol=0; icol < IndiceForColFin.size(); icol++){
                    std::cout <<"  "<<MesonrecoMatrix1.at(irow).at(icol)<<"  ";
                    }
                    std::cout <<endl<<endl;
                    }
                    int loc = findcloserepl(HH_showerptns[i], i+1, true, true, HH_showerptns, HH_thermal );
                    if(loc == 999999999){
                       std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                       HHparton fakep = HH_showerptns[i]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
                       Extraparton.add(fakep);
                       //somewhere, we need to make for loop to toss all partons to remnants list.
                    }
                    else if (loc > 0){HH_showerptns[loc - 1].acol(++maxtag); }
                    else if(loc < 0){HH_thermal[-loc - 1].acol(++maxtag); }

                    HH_showerptns[i].col(maxtag);
                  }



*/ //TODO:NEED TO BE DEALT IN THERMAL PARTON PROCESSING{04232020}
                    if(considering[0].acol() > 0 && considering[1].col() > 0){ //MAT + MAT
                    if(perm2[q2] > 0){
                    HH_showerptns[showerquarks[element[1]].par()].col(considering[0].acol());//now color tags from both partons are same
                    }
                    else{
                    HH_thermal[-element[1]].col(considering[0].acol());
                    }

                    }
                     else if(considering[0].acol() > 0){ // MAT + LBT/THERM

                     int loc = findcloserepl(considering[0], perm1[q1]+1, true, true, HH_showerptns, HH_thermal );
                     if(loc == 999999999){
                        //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                        HHparton fakep = considering[0]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
                        Extraparton.add(fakep);
                        //somewhere, we need to make for loop to toss all partons to remnants list.
                     }
                     else if (loc > 0){HH_showerptns[loc - 1].col(considering[0].acol()); }
                     else if(loc < 0){HH_thermal[-loc - 1].col(considering[0].acol()); }
                     /*
                     int thermsib = findthermalsibling(-element[1] , HH_thermal); // functon to find
                     if(thermsib < 0){
                     HH_thermal[thermsib-1].col(considering[0].acol());
                     }
                     else{HH_showerptns[thermsib+1].col(considering[0].acol()); }
                     */
                     }
                    else if(considering[1].col() > 0){ //LBT + MAT

                      int loc = findcloserepl(considering[1], perm2[q2], true, true, HH_showerptns, HH_thermal );
                      if(loc == 999999999){
                         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
                         HHparton fakep = considering[1]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
                         Extraparton.add(fakep);
                         //somewhere, we need to make for loop to toss all partons to remnants list.
                      }
                      else if (loc > 0){HH_showerptns[loc - 1].col(considering[0].acol()); }
                      else if(loc < 0){HH_thermal[-loc - 1].col(considering[0].acol()); }
                     /*
                     int thermsib = findthermalsibling(-element[0] , HH_thermal);
                     if(thermsib < 0){
                     HH_thermal[thermsib-1].acol(considering[1].col());
                     }
                     else{HH_showerptns[thermsib+1].acol(considering[1].col()); }
                     }
                     */
                    }
                    //now color tags from both partons are same
                  }


// DO the test without part below
					//*******************************************string repair functionality below*******************************************
					//determining which partons, of the 2 being considered, are endpoints
/*
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


*/
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

  //TODO: here we need new function to find thermal sibling for gluon
  // sol 1 : gluon loops.
  //sol 2 : decay gluon find pairs.{concern : energy conservation violated}
  // sol 3 : find two theraml siblings for gluon. declate fincloserepl twice. first pick quark and antiquarks to the 2nd.--> this is chosen
  for(int i = 0; i < HH_showerptns.num(); i++){
    if(HH_showerptns[i].col() == 0 && HH_showerptns[i].acol() == 0 && !HH_showerptns[i].is_used() ){
    if(HH_showerptns[i].id() == 21){
       int sel_out[2] = { 0 , 0 };

       findcloserepl_glu(HH_showerptns[i], i+1, true, true, HH_showerptns, HH_thermal, sel_out);

       if(sel_out[0] == 999999999 || sel_out[1] == 999999999){
         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
         HHparton fakeg = HH_showerptns[i]; fakeg.id(21);  fakeg.acol(++maxtag); fakeg.col(++maxtag);
         Extraparton.add(fakeg);
       }
       else if(sel_out[0] > 0){HH_showerptns[sel_out[0]-1].col(++maxtag); HH_showerptns[i].acol(maxtag);}
       else if(sel_out[0] < 0){HH_thermal[-sel_out[0]-1].col(++maxtag); HH_showerptns[i].acol(maxtag);}

       if(sel_out[1] > 0 && sel_out[1] != 999999999){HH_showerptns[sel_out[1]-1].acol(++maxtag); HH_showerptns[i].col(maxtag);}
       else if(sel_out[1] < 0){HH_thermal[-sel_out[1]-1].acol(++maxtag); HH_showerptns[i].col(maxtag);}

      }

    else if(HH_showerptns[i].id() > 0){
      int loc = findcloserepl(HH_showerptns[i], i+1, true, true, HH_showerptns, HH_thermal );
      if(loc == 999999999){
         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
         HHparton fakep = HH_showerptns[i]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
         Extraparton.add(fakep);
         //somewhere, we need to make for loop to toss all partons to remnants list.
      }
      else if (loc > 0){HH_showerptns[loc - 1].acol(++maxtag); }
      else if(loc < 0){HH_thermal[-loc - 1].acol(++maxtag); }

      HH_showerptns[i].col(maxtag);
    }
    else if(HH_showerptns[i].id() < 0){
      int loc = findcloserepl(HH_showerptns[i], i+1, true, true, HH_showerptns, HH_thermal );
      if(loc == 999999999){
         //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
         HHparton fakep = HH_showerptns[i]; fakep.id(-fakep.id());  fakep.col(++maxtag);
         Extraparton.add(fakep);
         //somewhere, we need to make for loop to toss all partons to remnants list.
      }
      else if (loc > 0){HH_showerptns[loc - 1].col(++maxtag); }
      else if(loc < 0){HH_thermal[-loc - 1].col(++maxtag); }

      HH_showerptns[i].acol(maxtag);
    }
    }

  }
	//sticking all unused partons into remnants; keeping order intact (a partially used gluon is replaced with it's unused quark)
	for(int i=0; i<HH_showerptns.num(); ++i){
		//if unused parton, write into remnants
		if(HH_showerptns[i].status() == 0){HH_remnants.add(HH_showerptns[i]); HH_remnants[HH_remnants.num() - 1].par(i); HH_showerptns[i].is_remnant(true);}
		//if 'partially' used gluon, write unused daughter quark into remnants
		else if(HH_showerptns[i].status() == -1){
			//finding the unused quark for this gluon and adding it to remnants (have to loop over as we only keep track of parents, not daughters)
			for(int j=0; j<showerquarks.num(); ++j){if(showerquarks[j].par() == i && !showerquarks[j].is_used()){
        if(showerquarks[j].col() == 0 && showerquarks[j].id() > 0){ //quark with zero colg tag is left,
          int loc = findcloserepl(showerquarks[j], j+1, true, true, HH_showerptns, HH_thermal );
          if(loc == 999999999){
             //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
             HHparton fakep = showerquarks[j]; fakep.id(-fakep.id());  fakep.acol(++maxtag);
             Extraparton.add(fakep);
             //somewhere, we need to make for loop to toss all partons to remnants list.
          }
          else if (loc > 0){HH_showerptns[loc - 1].acol(++maxtag); }
          else if(loc < 0){HH_thermal[-loc - 1].acol(++maxtag); }

          showerquarks[j].col(maxtag);
        }

        else if(showerquarks[j].acol() == 0 && showerquarks[j].id() < 0){ //antiquark with zero colg tag is left,
          int loc = findcloserepl(showerquarks[j], j+1, true, true, HH_showerptns, HH_thermal );
          if(loc == 999999999){
             //std::cout <<endl<<"Warning : extra parton used for string repair!!"<<endl;
             HHparton fakep = showerquarks[j]; fakep.id(-fakep.id());  fakep.col(++maxtag);
             Extraparton.add(fakep);
             //somewhere, we need to make for loop to toss all partons to remnants list.
          }
          else if (loc > 0){HH_showerptns[loc - 1].col(++maxtag); }
          else if(loc < 0){HH_thermal[-loc - 1].col(++maxtag); }

          showerquarks[j].acol(maxtag);
        }

        HH_remnants.add(showerquarks[j]); break;}}
			HH_showerptns[i].is_remnant(true);
		}
	}

	//appending the thermal partons used in the string repair functionality into remnants - order is NOT preserved...
	//is later sorted based not only on the string, but also on the position of the parton IN the string
	//can use the fact that the thermal partons needed either will be endpoints, or will be in a string with id > 0
	for(int i=0; i<HH_thermal.num(); ++i){
		//if this thermal parton is ... then add it to the remnants collection
		if(HH_thermal[i].is_used()){HH_thermal[i].status(1); HH_thermal[i].used_reco(true); continue;}
		if(HH_thermal[i].col() > 0 || HH_thermal[i].acol() > 0){HH_remnants.add(HH_thermal[i]); HH_remnants[HH_remnants.num() - 1].par(-i-1); HH_thermal[i].is_remnant(true);}
	}

  for(int i = 0; i < Extraparton.num(); ++i ){
    HH_remnants.add(Extraparton[i]);
  }

  // add fakeparton for color neutrality
  HHparton partonforneutrality;
  parton_collection fpartonforneutrality;
  fpartonforneutrality.add(partonforneutrality);
  if(fakepartoninfo.size() > 1){
    if(fakepartoninfo.at(0) > 0){
    fpartonforneutrality[0].col(fakepartoninfo.at(1));
    fpartonforneutrality[0].is_remnant(true);
    fpartonforneutrality[0].id(1);// generate fake up-quark
    //std::cout <<endl<<"fake parton is added with following tags : ("<<fpartonforneutrality[0].col()<<" , "<<fpartonforneutrality[0].acol()<<" ) "<<endl;
    HH_remnants.add(fpartonforneutrality[0]);
  }
  else if(fakepartoninfo.at(0) < 0){
    fpartonforneutrality[0].acol(fakepartoninfo.at(1));
    fpartonforneutrality[0].is_remnant(true);
    fpartonforneutrality[0].id(-1); //generate anti-up quark
    //std::cout <<endl<<"fake parton is added with following tags : ("<<fpartonforneutrality[0].col()<<" , "<<fpartonforneutrality[0].acol()<<" ) "<<endl;
    HH_remnants.add(fpartonforneutrality[0]);
  }
}
  fakepartoninfo.clear();
  fpartonforneutrality.clear();





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
  q1.col(glu.col()); q2.acol(glu.acol());

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

//!!!!! place into HybridHadronization.cc (after HybridHadronization::findthermalsibling(...){...} function block)
int HybridHadronization::findcloserepl(HHparton ptn, int iptn, bool lbt, bool thm, parton_collection& sh_lbt, parton_collection& therm){

	//should not happen
	if(iptn == 0 || ptn.id() == 21){throw std::runtime_error ("Parton index is incorrect (should not be 0 or gluon)");}

	//if the parton is thermal, and we're only looking at thermal partons, and the sibling was already found & not used, then return that.
	if((iptn<0) && thm && !lbt && !therm[ptn.sibling()].is_used() &&
	  (((therm[ptn.sibling()].id() > 0) && (therm[ptn.sibling()].col() != 0)) || ((therm[ptn.sibling()].id() < 0) && (therm[ptn.sibling()].acol() != 0))) &&
	  (iptn+1 != ptn.sibling())){return ptn.sibling();}

	//initializing vars
	int qrk_close = 999999999; double dist2min = 999999999999.;
	//checking thermal partons for closest parton
	if(thm){for(int i=0;i<therm.num();++i){
		//if the parton's color tag is not 0, then skip (is already repairing a string)
		if(((therm[i].id() > 0) && (therm[i].col() != 0)) || ((therm[i].id() < 0) && (therm[i].acol() != 0))){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if((ptn.id() * therm[i].id() > 0) || therm[i].is_used() || ((iptn<0) && (-i == iptn+1))){continue;}
		double distnow = ptn.posDif2(therm[i]) + (ptn.x_t()-therm[i].x_t())*(ptn.x_t()-therm[i].x_t());
		if(distnow < dist2min){qrk_close = -i-1; dist2min = distnow;}
	}}
	//checking lbt partons for closest parton
	if(lbt){for(int i=0;i<sh_lbt.num();++i){
		//if gluon, skip
		if(sh_lbt[i].id() == 21){continue;}
		//if the parton's color tag is not 0, then skip (not an lbt/martini? parton)
		if(((sh_lbt[i].id() > 0) && (sh_lbt[i].col() != 0)) || ((sh_lbt[i].id() < 0) && (sh_lbt[i].acol() != 0))){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if((ptn.id() * sh_lbt[i].id() > 0) || sh_lbt[i].is_used() || ((iptn>0) && (i == iptn-1))){continue;}
		double distnow = ptn.posDif2(sh_lbt[i]) + (ptn.x_t()-sh_lbt[i].x_t())*(ptn.x_t()-sh_lbt[i].x_t());
		if(distnow < dist2min){qrk_close = i+1; dist2min = distnow;}
	}}

	//positive return values indicate an lbt parton in the shower was found (the (i-1)th parton)
	//negative return values indicate a thermal parton was found (the -(i+1)th parton)
	return qrk_close;
}


//version to handle finding a q-qbar pair for lbt gluons
void HybridHadronization::findcloserepl_glu(HHparton ptn, int iptn, bool lbt, bool thm, parton_collection& sh_lbt, parton_collection& therm, int sel_out[]){

	//should not happen
	if(iptn == 0 || ptn.id() != 21){throw std::runtime_error ("Parton index is incorrect (should not be 0, or anything other than gluon)");}
	if(ptn.is_thermal()){throw std::runtime_error ("Parton is thermal (should not be so)");}

	//initializing vars
	int qrk_close = 999999999; double dist2min = 999999999999.;
	//checking thermal partons for closest quark
	if(thm){for(int i=0;i<therm.num();++i){
		//if the parton's color tag is not 0, then skip (is already repairing a string), or skip if antiquark
		if(((therm[i].id() > 0) && (therm[i].col() != 0)) || (therm[i].id() < 0)){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if(therm[i].is_used() || ((iptn<0) && (-i == iptn+1))){continue;}
		double distnow = ptn.posDif2(therm[i]) + (ptn.x_t()-therm[i].x_t())*(ptn.x_t()-therm[i].x_t());
		if(distnow < dist2min){qrk_close = -i-1; dist2min = distnow;}
	}}
	//checking lbt partons for closest parton
	if(lbt){for(int i=0;i<sh_lbt.num();++i){
		//if gluon, skip
		if(sh_lbt[i].id() == 21){continue;}
		//if the parton's color tag is not 0, then skip (not an lbt/martini? parton), or skip if antiquark
		if(((sh_lbt[i].id() > 0) && (sh_lbt[i].col() != 0)) || (sh_lbt[i].id() < 0)){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if(sh_lbt[i].is_used() || ((iptn>0) && (i == iptn-1))){continue;}
		double distnow = ptn.posDif2(sh_lbt[i]) + (ptn.x_t()-sh_lbt[i].x_t())*(ptn.x_t()-sh_lbt[i].x_t());
		if(distnow < dist2min){qrk_close = i+1; dist2min = distnow;}
	}}

	//positive return values indicate an lbt parton in the shower was found (the (i-1)th parton)
	//negative return values indicate a thermal parton was found (the -(i+1)th parton)
	sel_out[0]=qrk_close;

	qrk_close = 999999999; dist2min = 999999999999.;
	//checking thermal partons for closest antiquark
	if(thm){for(int i=0;i<therm.num();++i){
		//if the parton's color tag is not 0, then skip (is already repairing a string), or skip if quark
		if((therm[i].id() > 0) || ((therm[i].id() < 0) && (therm[i].acol() != 0))){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if(therm[i].is_used() || ((iptn<0) && (-i == iptn+1))){continue;}
		double distnow = ptn.posDif2(therm[i]) + (ptn.x_t()-therm[i].x_t())*(ptn.x_t()-therm[i].x_t());
		if(distnow < dist2min){qrk_close = -i-1; dist2min = distnow;}
	}}
	//checking lbt partons for closest parton
	if(lbt){for(int i=0;i<sh_lbt.num();++i){
		//if gluon, skip
		if(sh_lbt[i].id() == 21){continue;}
		//if the parton's color tag is not 0, then skip (not an lbt/martini? parton), or skip if quark
		if((sh_lbt[i].id() > 0) || ((sh_lbt[i].id() < 0) && (sh_lbt[i].acol() != 0))){continue;}
		//if the parton is not a partner (eg. both quarks/antiquarks), or is used, or is the same as the input parton, skip it.
		if(sh_lbt[i].is_used() || ((iptn>0) && (i == iptn-1))){continue;}
		double distnow = ptn.posDif2(sh_lbt[i]) + (ptn.x_t()-sh_lbt[i].x_t())*(ptn.x_t()-sh_lbt[i].x_t());
		if(distnow < dist2min){qrk_close = i+1; dist2min = distnow;}
	}}

	//positive return values indicate an lbt parton in the shower was found (the (i-1)th parton)
	//negative return values indicate a thermal parton was found (the -(i+1)th parton)
	sel_out[1]=qrk_close;
}

//prepares remnant partons/strings for PYTHIA string hadronization
//sorts strings, ensures strings are in 'valid' configurations, assigns color/anticolor tags
//TODO: this might be where to use thermal partons to enforce color neutrality
void HybridHadronization::stringprep(parton_collection& SP_remnants, parton_collection& SP_prepremn, bool cutstr){

  //dignostic measure
  //std::cout <<endl<<"below is all the colors in the Tempjunctin lists"<<endl;
  for(int ijunc=0; ijunc<Tempjunctions.size(); ijunc++){
    //std::cout <<" { "<<Tempjunctions.at(ijunc).at(0).at(0)<<" , "<<Tempjunctions.at(ijunc).at(1).at(1)<<" , "<<Tempjunctions.at(ijunc).at(2).at(1)<<" , "<<Tempjunctions.at(ijunc).at(3).at(1)<<" } "<<endl;
}
//Declare essential vectors for string repair
    vector<vector<vector<HHparton>>> JuncStructure;
    vector<vector<HHparton>> JuncLegs; // vector of all junction ( Junction Num, Leg1, Leg2, Leg3)
    vector<HHparton> Leg1; // partons in the legs for juncion formation
    vector<HHparton> Leg2;
    vector<HHparton> Leg3;
    vector<int> Legconsidering;//address of partons in SP_remnants
    vector<vector<vector<int>>> IMStructure1; //Intermediate structure for saving the indices in Tempjunction and Corresponding Leg Number, this will be used for cutting string for appending info to PYTHIA
    vector<vector<int>> IMStructure2; // this will contain Tempjunction Indice and -+1 and leg numbers
    //ex. 2nd lef in 3rd temp anti junction is same as 1st leg 2nd tempjunction. 3,-1,0,0   3,1,2,0 would be saved in the vectorIMS1
    vector<int> IMStructure3; // 4 element vector of Tempjunction order, kind and leg address in JuncStructure vector
    vector<vector<HHparton>> Recombearly1; //these partons are in the junction with three or two shared legs with others , and they will be recombined into baryon to cut? or arrange the string for PYTHIA to understand input come from this code
    vector<HHparton> Recombearly2;
    vector<vector<vector<HHparton>>> Dijunction1; // these partons are in the junction with two legs shared legs!
    vector<vector<HHparton>> Dijunction2;
    //Vector of Indices for Dijunction1 vector to be ordered for invoking Pythia
    //vector<vector<vector<int>>> DinjunctionInfo1;
    vector<vector<int>> DijunctionInfo1;
    vector<int> DijunctionInfo2; //  { -1: antiJ, +1 : J , 0 : shared leg}
    //these tags should be assigned to corresponding legs in Dijunction1 vector! so DijunctionInfo2 has five elements{since, dijunction structure has five legs!}


    vector<vector<vector<HHparton>>> Singlejunction1; // easiest case! just single string
    vector<vector<HHparton>> Singlejunction2;
    vector<vector<HHparton>> Tailoredstring1; // when there are the junction with two shared legs, it also should be recombined into baryon and one string will remain after, which would be saved in this vector
    vector<HHparton> Tailoredstring2;
    vector<int> realjuncindice; // vector to save the indice of Tempjunction when real junction is formed by three initiating particle

    vector<HHparton> finalstring; // final space for all of the remnant particles being corrected by fake parton addition


//Since MATTER shower is not color neutral, we need add condition for the color neutrality.
//First of all, check if the Tempjunction contains color tags of the gluon located at endpoint of string.


/* just reference below
                    HH_showerptns[showerquarks[element[0]].par()].col(HH_showerptns[showerquarks[element[1]].par()].acol());
                    vector<vector<vector<int>>> Tempjunctions; // vector of all tempjunctions
                    vector<vector<int>> JunctionInfo; // vector of one junction's color tag and particle info
                    vector<int> IdColInfo1;
                    vector<int> IdColInfo2;
                    vector<int> IdColInfo3;
                    vector<int> IdColInfo4;
*/

   for(int ijunc=0; ijunc < Tempjunctions.size(); ++ijunc) { //Let's make Legs for Final string by verifying Tempjunctions! starting from the kind of junction(-1 or +1)
         //std::cout <<endl<<endl<<"  Let's see candidates for (anti)junction  "<<endl;
         //std::cout <<" Candidate " << ijunc+1 <<endl;
         //std::cout <<" [ (  " << Tempjunctions.at(ijunc).at(0).at(0) <<" , " << Tempjunctions.at(ijunc).at(0).at(1) <<" ), " ;
         //std::cout <<"  (  " << Tempjunctions.at(ijunc).at(1).at(0) <<" , " << Tempjunctions.at(ijunc).at(1).at(1) <<" ), " ;
         //std::cout <<"  (  " << Tempjunctions.at(ijunc).at(2).at(0) <<" , " << Tempjunctions.at(ijunc).at(2).at(1) <<" ), " ;
         //std::cout <<"  (  " << Tempjunctions.at(ijunc).at(3).at(0) <<" , " << Tempjunctions.at(ijunc).at(3).at(1) <<" ), " ;
         //std::cout <<" ] "<<endl;

 //std::cout <<endl<<" I'm Working !! from line 2027" <<endl;
           vector<int> correction; // for the case of 1 or 2 initiating particles, we need to remember the partons and get the used_junction tag to the original
           if(Tempjunctions.at(ijunc).at(0).at(0) == -1) { // check anti-color tags in SP_remnants partons!! to form anti-junction
               for(int irem=0; irem < SP_remnants.num(); ++irem) { // searching through all remnant particles
//std::cout <<endl<<" Let's see " << irem+1 <<"th remnant particle with "<< "pid : " << SP_remnants[irem].pid() <<" color : "<< SP_remnants[irem].col()<< " anti-color : " << SP_remnants[irem].acol() <<endl;
                  if(SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(1).at(1)){
                  Leg1.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  //std::cout <<endl<<" Leg 1 has initiating particle " <<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].acol()<<endl<<endl;
                  Legconsidering.push_back(irem);
                  }
                  if(SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(2).at(1)){
                  Leg2.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  //std::cout <<endl<<" Leg 2 has initiating particle " <<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].acol()<<endl<<endl;
                  Legconsidering.push_back(irem);
                  }
                  if(SP_remnants[irem].acol() == Tempjunctions.at(ijunc).at(3).at(1)){
                  Leg3.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  //std::cout <<endl<<" Leg 3 has initiating particle " <<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].acol()<<endl<<endl;
                  Legconsidering.push_back(irem);
                  }
                  if(Legconsidering.size() !=3 ){// if only one or two particles found in remnants, we need to set the tag to the original
                    for(int icor = 0; icor < correction.size(); icor++){
                    SP_remnants[correction.at(icor)].used_junction(false);
                    //dignostic measure
                    //std::cout <<endl<<"the "<<icor<<" th particle with color tag of "<<SP_remnants[icor].col()<<" , "<<SP_remnants[icor].acol()<<" are turned into unused particle"<<endl;
                    }
                    correction.clear();
                  }

               }
           }
           else if(Tempjunctions.at(ijunc).at(0).at(0) == 1) { // check color tags in SP_remnants partons!! to form junction
               for(int irem=0; irem < SP_remnants.num(); ++irem) { // searching through all remnant particles
//std::cout <<endl<<" Let's see " << irem+1 <<"th remnant particle with "<< "pid : " << SP_remnants[irem].pid() <<" color : "<< SP_remnants[irem].col()<< " anti-color : " << SP_remnants[irem].acol() <<endl;
                  if(SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(1).at(1)){
                  Leg1.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  Legconsidering.push_back(irem);
                  //std::cout <<endl<<" Leg 1 has initiating particle "<<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].col()<<endl<<endl;
                  }

                  if(SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(2).at(1)){
                  Leg2.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  //std::cout <<endl<<" Leg 2 has initiating particle " <<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].col()<<endl<<endl;
                  Legconsidering.push_back(irem);
                  }

                  if(SP_remnants[irem].col() == Tempjunctions.at(ijunc).at(3).at(1)){
                  Leg3.push_back(SP_remnants[irem]);
                  correction.push_back(irem);
                  SP_remnants[irem].used_junction(true);
                  //std::cout <<endl<<" Leg 3 has initiating particle " <<" and the pID is "<< SP_remnants[irem].pid()<<" and the color tag : "<<SP_remnants[irem].col()<<endl<<endl;
                  Legconsidering.push_back(irem);
                  }
                 if(Legconsidering.size() !=3 ){
                    for(int icor = 0; icor < correction.size(); icor++){
                    SP_remnants[correction.at(icor)].used_junction(false);
//dignostic measure
                    //std::cout <<endl<<"the "<<icor<<" th particle with color tag of "<<SP_remnants[icor].col()<<" , "<<SP_remnants[icor].acol()<<" are turned into unused particle"<<endl;
                    }
                    correction.clear();
                  }

               }
           }
// So far, we put all the starting particles in three legs for (anti)junction
 //std::cout <<endl<<" I'm Working !! from line 2061" <<endl;
 //std::cout <<endl<<" number of quarks relevant to real-junction is "<< Legconsidering.size()<<" " <<endl;

         if(Legconsidering.size() == 3) {
           realjuncindice.push_back(ijunc);
           //std::cout <<endl<<"the indice added for string repair : "<<ijunc<<endl;
           if(Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg1.at(0).col() != 0 && Leg1.at(0).acol() != 0){ // starting to form the Leg1 for anti_junction by Color Tag Tracing. If first particle is (anti)quark, no need to trace
               for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                      if( (Leg1.back().col() != 0) && (Leg1.back().col() == SP_remnants[icf].acol()) ){
                        //std::cout <<endl<<"particle added to Leg 1 with color tag  ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg1.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }
                      if( Leg1.back().col() == 0 || Leg1.back().acol() == 0 ){ break; }
                  }
                }

           }
           else if(Tempjunctions.at(ijunc).at(0).at(0) == 1 && Leg1.at(0).col() != 0 && Leg1.at(0).acol() != 0){ // starting to form the Leg1 junction by Color Tag Tracing. If first particle is quark, no need to trace
               for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                      if( (Leg1.back().acol() != 0) && (Leg1.back().acol() == SP_remnants[icf].col()) ){
                        //std::cout <<endl<<"particle added to Leg 1 with color tag ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg1.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }
                      if(Leg1.back().col() == 0 || Leg1.back().acol() == 0){break;}
                  }
                }

           }

           if(Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg2.at(0).col() != 0 && Leg2.at(0).acol() != 0){ // starting to form the Leg2 for anti_junction by Color Tag Tracing. If first particle is (anti)quark, no need to trace
               for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                      if(  ( Leg2.back().col() != 0 )  && (Leg2.back().col() == SP_remnants[icf].acol()) ){
                        //std::cout <<endl<<"particle added to Leg 2 with color tag ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg2.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }
                      if( Leg2.back().col() == 0 || Leg2.back().acol() == 0){break;}

                  }
                }
           }
           else if(Tempjunctions.at(ijunc).at(0).at(0) == 1 && Leg2.at(0).col() != 0 && Leg2.at(0).acol() != 0){ // starting to form the Leg2 junction by Color Tag Tracing. If first particle is quark, no need to trace
                for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                      if(   (Leg2.back().acol() != 0) && (Leg2.back().acol() == SP_remnants[icf].col())  ){
                        //std::cout <<endl<<"particle added to Leg 2 with color tag ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg2.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }
                      if( Leg2.back().col() == 0 || Leg2.back().acol() == 0){break;}
                  }
                }

           }

           if(Tempjunctions.at(ijunc).at(0).at(0) == -1 && Leg3.at(0).col() != 0 && Leg3.at(0).acol() != 0){ // starting to form the Leg3 for anti_junction by Color Tag Tracing. If first particle is (anti)quark, no need to trace
                for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                      if(  (Leg3.back().col() != 0) && (Leg3.back().col() == SP_remnants[icf].acol()) ){
                        //std::cout <<endl<<"particle added to Leg 3 with color tag ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg3.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }

                     if( Leg3.back().col() == 0 || Leg3.back().acol() == 0){break;}

              }
            }

           }
           else if(Tempjunctions.at(ijunc).at(0).at(0) == 1 && Leg3.at(0).col() != 0 && Leg3.at(0).acol() != 0){ // starting to form the Leg3 junction by Color Tag Tracing. If first particle is quark, no need to trace
                for(int iloop=0; iloop < SP_remnants.num(); iloop++){
                  for(int icf = 0; icf < SP_remnants.num(); ++icf) {
                    //std::cout <<"candidate particle is "<<SP_remnants[icf].pid()<<" , "<<SP_remnants[icf].col()<<" , "<<SP_remnants[icf].acol()<<" , "<<SP_remnants[icf].used_junction()<<endl;
                    //std::cout <<"and standard particle is "<<Leg3.back().pid()<<" , "<<Leg3.back().col()<<" , "<<Leg3.back().acol()<<" , "<<Leg3.back().used_junction()<<endl;
                      if( (Leg3.back().acol() != 0) && (Leg3.back().acol() == SP_remnants[icf].col()) ){
                        //std::cout <<endl<<"particle added to Leg 3 with color tag ( "<<SP_remnants[icf].col() <<" and "<<SP_remnants[icf].acol()<<" ) "<<endl;
                      Leg3.push_back(SP_remnants[icf]);
                      SP_remnants[icf].used_junction(true);
                      }

                      if( Leg3.back().col() == 0 || Leg3.back().acol() == 0 ){break;}
                  }
                }

           }

         JuncLegs.push_back(Leg1);
         JuncLegs.push_back(Leg2);
         JuncLegs.push_back(Leg3);

         JuncStructure.push_back(JuncLegs); // now junction structure with three complete legs is saved, so clear previous infos for considering next junction candidate.
        }
        /*
        else { // if not all three legs have initiating particles, make the tags to the original
            for(int icon=0; icon < Legconsidering.size(); ++icon){
            SP_remnants[icon].used_junction(false);
            }
        }
        */



         JuncLegs.clear();
         Leg1.clear();
         Leg2.clear();
         Leg3.clear();
         Legconsidering.clear();
      }
//TODO: just testing new procedure below
//So far, all the information for (Anti)Junction and Legs are formed by color tag tracing. Let's Check this out


//dignostic measure
    for(int ijuncstr=0; ijuncstr < JuncStructure.size(); ++ijuncstr) {
         //std::cout <<endl<<endl<<"  Let's see Temporary junction structure   "<<endl;
         //std::cout <<" Temp Junction :" << ijuncstr+1 <<endl;
         //std::cout <<"Leg 1 : [ ";
         for(int ileg1=0; ileg1 < JuncStructure.at(ijuncstr).at(0).size(); ++ileg1){
         //std::cout <<" (  " << JuncStructure.at(ijuncstr).at(0)[ileg1].col() <<" , " << JuncStructure.at(ijuncstr).at(0)[ileg1].acol() <<" ), " ;
         }
         //std::cout <<" ] " <<endl;
         //std::cout <<"Leg 2 : [ ";
         for(int ileg2=0; ileg2 < JuncStructure.at(ijuncstr).at(1).size(); ++ileg2){
         //std::cout <<" (  " << JuncStructure.at(ijuncstr).at(1)[ileg2].col() <<" , " << JuncStructure.at(ijuncstr).at(1)[ileg2].acol() <<" ), " ;
         }
         //std::cout <<" ] " <<endl;
         //std::cout <<"Leg 3 : [ ";
         for(int ileg3=0; ileg3 < JuncStructure.at(ijuncstr).at(2).size(); ++ileg3){
         //std::cout <<" (  " << JuncStructure.at(ijuncstr).at(2)[ileg3].col() <<" , " << JuncStructure.at(ijuncstr).at(2)[ileg3].acol() <<" ), " ;
         }
         //std::cout <<" ] " <<endl;

    }

    // Now, we will find the shared legs and sorting them into the corresponding vectors, the information from JuncStructure vector and realjuncindice vector work together here.
    // if the tags of first particle in a string and end particle in another string are same, they are shared Leg!!
    for(int irep1 = 0; irep1 < JuncStructure.size(); irep1++){
      IMStructure3.push_back(realjuncindice.at(irep1));
      IMStructure3.push_back(Tempjunctions.at(realjuncindice.at(irep1)).at(0).at(0));
      IMStructure3.push_back(irep1);// room for other usage
      IMStructure3.push_back(0); // tag 1: dijinction, tag0: other cases
      IMStructure2.push_back(IMStructure3);
      IMStructure3.clear();
      for(int irep2 = 0; irep2 < 3; irep2++){
      for(int irep3 = 0; irep3 < JuncStructure.size(); irep3++){
        for(int irep4 = 0;  irep4 < 3; irep4++){
          if( (irep1 != irep3) && (JuncStructure.at(irep1).at(irep2).at(0).col() == JuncStructure.at(irep3).at(irep4).back().col()) &&
          (JuncStructure.at(irep1).at(irep2).at(0).acol() == JuncStructure.at(irep3).at(irep4).back().acol())){
            IMStructure3.push_back(irep1);
            IMStructure3.push_back(irep2);
            IMStructure3.push_back(irep3);
            IMStructure3.push_back(irep4); // the first pair means irep2_th Leg in irep1_th junction is linked with irep4_th Leg in irep3_th Tempjunction
            IMStructure2.push_back(IMStructure3);
            IMStructure3.clear();
          }

        }
      }
    }
    IMStructure1.push_back(IMStructure2);
    IMStructure2.clear();
  }
//dignostic measure!
//std::cout <<endl;
for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
  for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
    //std::cout <<" ( ";
    for(int icheck3 = 0; icheck3 < IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){
      //std::cout <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
    }
  //std::cout <<" ) ";

  }
    //std::cout <<endl;
}

//TODO:since the relation between junctions is defined, Let's sorting Junctions based on shared legs

for(int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++){
  //structure for 3-shared particles find shared legs and tossing them into the baryon formation
  if(IMStructure1.at(iloop1).size() == 4){
    for(int iloop2 = 1; iloop2 < IMStructure1.at(iloop1).size(); iloop2++){
      for(int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++){
        for(int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size(); iloop4++){
          if( (IMStructure1.at(iloop1).at(iloop2).at(2) == IMStructure1.at(iloop3).at(iloop4).at(0)) &&
               (IMStructure1.at(iloop1).at(iloop2).at(3) == IMStructure1.at(iloop3).at(iloop4).at(1))  ){
                 vector<int> testing;
                 testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
                 testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
                 testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
                 testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(1)); //now temporary vector is created by this procedure.
                 std::vector<vector<int>>::iterator it = std::find(IMStructure1.at(iloop3).begin(), IMStructure1.at(iloop3).end(), testing);
                 int loc = std::distance(IMStructure1.at(iloop3).begin(), it);//dignositc measure for the Process
                 IMStructure1.at(iloop3).erase(it);
                 testing.clear(); //now the one of the leg is eliminated in the list, so there is no overlapped leg to be tossed into baryon formation list.
               }
               if( (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0).col() != 0 ) &&
               (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0).acol() != 0 ) ) {
                 parton_collection tempqpair;
                 gluon_decay(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0), tempqpair ); // the first particle has col tag and the second had acol tag
                 //std::vector<HHparton>::iterator tempit1 = JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).begin();
                 //JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).erase(tempit1);
                 if(IMStructure1.at(iloop1).at(0).at(1) == -1){
                   Recombearly2.push_back(tempqpair[1]); //add anti-particle to form anti baryon
                   //std::cout <<endl<<" the particle with tags { "<< JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().col()<<" , "
                   //<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).back().acol()<<" } is changed into ";

                   JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).pop_back();

                   JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).push_back(tempqpair[0]); // add remain parton to the original leg
                   //std::cout <<" { " << tempqpair[0].col()<<" , "<< tempqpair[0].acol()<<" }"<<endl;
                   //std::cout <<" and "<<" { " << tempqpair[1].col()<<" , "<< tempqpair[1].acol()<<" } "<<" is gone to be recombined"<<endl;

                 }
                 if(IMStructure1.at(iloop1).at(0).at(1) == 1){
                   Recombearly2.push_back(tempqpair[0]); //add particle to form baryon
                   JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).pop_back();
                   JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).push_back(tempqpair[1]); // add remain parton to the original leg
                 }
                 tempqpair.clear();


               }

        }
      }
    }
    Recombearly1.push_back(Recombearly2);
    Recombearly2.clear();
  }
}

//One thing to remind : JuncStructure has information of parton class, and IMStructure has informations about these partons or junctions with same indice with JuncStructure
for(int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++){

 if(IMStructure1.at(iloop1).size() == 3){ // 2-shared leg case {basic info, shared leg1's info, shared leg2's info}
 // loop for pop out the unshared leg in the junction structure.
     int leftleg; // indice of unshared leg which will remain after the Early recombinaton
     int SI = 0; // Sum of the Indice of Leg{if 1 = 0+1, Leg3 is not shared, if 3 = 1+2, Leg 0 is not shared}
     for(int itail = 1; itail < IMStructure1.at(iloop1).size(); itail++){
       SI = SI + IMStructure1.at(iloop1).at(itail).at(1);
     }
     //std::cout <<endl<<"SI = "<<SI<<endl;
     vector<int> indicevec; // vector for indice, we could use if statement, but when if and for statement are repeated, possibly forced execution occurs regardless of if statement, so take most defensive measure.
     indicevec.push_back(2);
     indicevec.push_back(1);
     indicevec.push_back(0);
     leftleg = indicevec[SI-1];
     //below was tested, but showed forced execution{that is all of the Three statements executed so that left leg was always 0}
     //if(SI = 1 ){ leftleg = 2; SI = 0;}
     //if(SI = 2 ){ leftleg = 1; SI = 0;}
     //if(SI = 3 ){ leftleg = 0; SI = 0;}
     //std::cout <<endl<<"leftleg number is "<<leftleg<<endl;

     int jidentity = IMStructure1.at(iloop1).at(0).at(1); //junction identity 1 or -1 {junction or anti juncttion}
     //std::cout <<endl<<"junction identity is "<<jidentity<<endl;
     int juncnum = IMStructure1.at(iloop1).at(1).at(0); //junction number to be dealt with
     vector<HHparton> tempstring = JuncStructure.at(juncnum).at(leftleg); //declare the unshared leg.
     if(jidentity = 1 && tempstring.at(0).col() != 0 && tempstring.at(0).acol() != 0){
     parton_collection tempqpair;
     gluon_decay( tempstring.at(0), tempqpair);
     tempstring.erase(tempstring.begin());
     tempstring.insert(tempstring.begin() , tempqpair[1] );
     Tailoredstring1.push_back(tempstring);
    }
    if(jidentity = -1 && tempstring.at(0).col() != 0 && tempstring.at(0).acol() != 0){
     parton_collection tempqpair;
     gluon_decay( tempstring.at(0), tempqpair);
     tempstring.erase(tempstring.begin());
     tempstring.insert(tempstring.begin() , tempqpair[0] );
     Tailoredstring1.push_back(tempstring);
    }
    indicevec.clear();


     for(int iloop2 = 1; iloop2 < IMStructure1.at(iloop1).size(); iloop2++){
       for(int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++){
         for(int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size(); iloop4++){
           if( (IMStructure1.at(iloop1).at(iloop2).at(2) == IMStructure1.at(iloop3).at(iloop4).at(0)) &&
                (IMStructure1.at(iloop1).at(iloop2).at(3) == IMStructure1.at(iloop3).at(iloop4).at(1))  ){
                  vector<int> testing;
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(1)); //now temporary vector is created by this procedure.
                  std::vector<vector<int>>::iterator it = std::find(IMStructure1.at(iloop3).begin(), IMStructure1.at(iloop3).end(), testing);
                  int loc = std::distance(IMStructure1.at(iloop3).begin(), it);
                  IMStructure1.at(iloop3).erase(it);
                  ////std::cout <<endl<<" the IMS vector ( "<<testing[1]<<" , "<<testing[2]<<" , "<<testing[3]<<" , "<<testing[4]<<" is to be deleted!"<<endl;
                  testing.clear(); //now the one of the leg is eliminated in the list, so there is no overlapped leg to be tossed into baryon formation list.
                }
                if( (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0).col() != 0 ) &&
                (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0).acol() != 0 ) ) {
                  parton_collection tempqpair;
                  gluon_decay(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).at(0), tempqpair ); // the first particle has col tag and the second had acol tag
                  //std::vector<HHparton>::iterator tempit1 = JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).begin();
                  //JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(IMStructure1.at(iloop1).at(iloop2).at(1)).erase(tempit1);
                  if(IMStructure1.at(iloop1).at(0).at(1) == -1){
                    Recombearly2.push_back(tempqpair[1]); //add anti-particle to form anti baryon
                    JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).pop_back();
                    JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).push_back(tempqpair[0]); // add remain parton to the original leg
                  }
                  if(IMStructure1.at(iloop1).at(0).at(1) == 1){
                    Recombearly2.push_back(tempqpair[0]); //add particle to form baryon
                    JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).pop_back();
                    JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(IMStructure1.at(iloop1).at(iloop2).at(3)).push_back(tempqpair[1]); // add remain parton to the original leg
                  }
                }
                /*
                //dignostic measure!
                std::cout <<endl;
                for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
                  for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
                    std::cout <<" ( ";
                    for(int icheck3 = 0; icheck3 < IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){
                      std::cout <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
                    }
                  std::cout <<" ) ";

                  }
                    std::cout <<endl;
                }
                */

                for(int icheck1 = 0; icheck1 < 3; icheck1++){ // procedure to find non-shared leg
                  if( (IMStructure1.at(iloop1).size() == 3) && (icheck1 != IMStructure1.at(iloop1).at(1).at(1)) && (icheck1 != IMStructure1.at(iloop1).at(2).at(1)) ){
                    //std::cout <<endl<<"Let's check iloop1 through iloop4! which are ( "<<iloop1<<" , "<<iloop2<<" , "<<iloop3<<" , "<<iloop4<<" , "<<icheck1<<" ) "<<endl;

                    if ( (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).col() != 0) &&
                  (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).acol() != 0) ){//first particle is gluon, in this case, we need to decay it and rearrange the particles in leg
                    parton_collection tempqpair2;
                    gluon_decay( JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0) , tempqpair2 );
                    if(IMStructure1.at(iloop1).at(0).at(1) == -1){
                      Recombearly2.push_back(tempqpair2[1]); //add anti-particle to form anti baryon
                      JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).pop_back();
                      JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).push_back(tempqpair2[0]); // add remain parton to the original leg
                    }
                    if(IMStructure1.at(iloop1).at(0).at(1) == 1){
                      Recombearly2.push_back(tempqpair2[0]); //add particle to form baryon
                      JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).pop_back();
                      JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).push_back(tempqpair2[1]); // add remain parton to the original leg
                    }
                      tempqpair2.clear();
                  }

                    if ( (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).col() == 0) ||
                  (JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).acol() == 0) ){ // the remained leg only has quark or antiquark, so directly toss this for baryon formation
                    //std::cout <<endl<<" the particle at "<<IMStructure1.at(iloop1).at(iloop2).at(0)<< " , "<<icheck1<<" is attached to be baryon"<<endl;
                    //std::cout <<"and it's col tag :"<<JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).col()<<" , "<<
                    //JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0).acol()<<" ) "<<endl;
                  Recombearly2.push_back(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1).at(0));
                  }

                    //Tailoredstring1.push_back(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(icheck1));



                  }
                }
         }
       }
     }
     Recombearly1.push_back(Recombearly2);
     Recombearly2.clear();


   }
 }

 for(int iloop1 = 0; iloop1 < IMStructure1.size(); iloop1++){

 if(IMStructure1.at(iloop1).size() == 2 && IMStructure1.at(iloop1).at(0).at(3) == 0){ // Dijunction Structure


     for(int iloop2 = 1; iloop2< IMStructure1.at(iloop1).size(); iloop2++){
       for(int iloop3 = 0; iloop3 < IMStructure1.size(); iloop3++){
         for(int iloop4 = 1; iloop4 < IMStructure1.at(iloop3).size(); iloop4++){
           if( (IMStructure1.at(iloop1).at(iloop2).at(2) == IMStructure1.at(iloop3).at(iloop4).at(0)) &&
                (IMStructure1.at(iloop1).at(iloop2).at(3) == IMStructure1.at(iloop3).at(iloop4).at(1))  ){ // the case where we found the shared legs pair.
              if(IMStructure1.at(iloop1).at(0).at(3) == 0) {
                  for(int idijunc1 = 0; idijunc1 < 3; idijunc1++){ // first, put three legs of first junction
                      Dijunction2.push_back(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(0)).at(idijunc1));
                  }
                  for(int idijunc2 = 0; idijunc2 < 3; idijunc2++){ // second, exclude shared leg in second junction to be merged
                    if( idijunc2 != IMStructure1.at(iloop1).at(iloop2).at(3)){
                      Dijunction2.push_back(JuncStructure.at(IMStructure1.at(iloop1).at(iloop2).at(2)).at(idijunc2));
                    }
                  }
                }
              if( IMStructure1.at(iloop3).size() == 2 && IMStructure1.at(IMStructure1.at(iloop1).at(1).at(2)).size() == 2){ // this means that shared junction is also shares only one legs with othe junction
                  //std::cout <<endl<<"the shared leg of "<<iloop1+1<<"th junction with "<< IMStructure1.at(iloop1).at(1).at(2) + 1 <<" is two"<<endl;
                  vector<int> testing;
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(2));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(3));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(0));
                  testing.push_back(IMStructure1.at(iloop1).at(iloop2).at(1)); //now temporary vector is created by this procedure.
                  std::vector<vector<int>>::iterator it = std::find(IMStructure1.at(iloop3).begin(), IMStructure1.at(iloop3).end(), testing);
                  int loc = std::distance(IMStructure1.at(iloop3).begin(), it);//dignositc measure for the Process
                  IMStructure1.at(iloop3).erase(it);
                  testing.clear();
                 //now the one of the leg is eliminated in the list, so there is no overlapped leg to be tossed into baryon formation list.

                  Dijunction1.push_back(Dijunction2); //finally, five legs are added in the vector to form dijunction structure.
                  Dijunction2.clear(); // clear temp vector for next procedure.
                  IMStructure1.at(iloop1).at(0).pop_back();
                  IMStructure1.at(iloop1).at(0).push_back(1); // tag representing dijunction structure.
                  IMStructure1.at(iloop3).at(0).pop_back();
                  IMStructure1.at(iloop3).at(0).push_back(1); // tag representing dijunction structure.
                }
              if(IMStructure1.at(IMStructure1.at(iloop1).at(1).at(2)).size() == 3 ||
                 IMStructure1.at(IMStructure1.at(iloop1).at(1).at(2)).size() == 4 ){
                   //std::cout <<endl<<"the shared leg of "<<iloop1+1<<"th junction with " <<IMStructure1.at(iloop1).at(1).at(2) + 1 <<" is not two"<<endl;
                   IMStructure1.at(iloop1).at(0).pop_back();
                   IMStructure1.at(iloop1).at(0).push_back(0); // tag representing dijunction structure.
                   IMStructure1.at(iloop3).at(0).pop_back();
                   IMStructure1.at(iloop3).at(0).push_back(0); // tag representing dijunction structure.
                    Dijunction2.clear(); }
            }
         }
       }
     }


 }
}


//}// end parentheses IMStructure1 list loop

for(int irepair = 0; irepair < IMStructure1.size(); irepair++) {
if(IMStructure1.at(irepair).size() == 1 && IMStructure1.at(irepair).at(0).at(3) == 0  ){ // size one means there are no shared legs to the junction, so put them directly to the single junction vector
  Singlejunction2.push_back(JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(0));
  Singlejunction2.push_back(JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(1));
  Singlejunction2.push_back(JuncStructure.at(IMStructure1.at(irepair).at(0).at(2)).at(2));
  Singlejunction1.push_back(Singlejunction2);
  Singlejunction2.clear();
}//After all, all these vectors of junction legs will be tossed into PYTHIA to be hadronized. before that we need to set Mother Daughter tag for PYTHIA

} // collecting all single junctions



//dignostic measure{IMStructure1}
//std::cout <<endl;
for(int icheck1 = 0; icheck1 < IMStructure1.size(); icheck1++){
  for(int icheck2 = 0; icheck2 < IMStructure1.at(icheck1).size(); icheck2++){
    //std::cout <<" ( ";
    for(int icheck3 = 0; icheck3 < IMStructure1.at(icheck1).at(icheck2).size(); icheck3++){
      //std::cout <<IMStructure1.at(icheck1).at(icheck2).at(icheck3)<<" , ";
    }
  //std::cout <<" ) ";

  }
    //std::cout <<endl;
}


//now we need to form fake partons and baryons for appending particles to PYTHIA{for Dijunction Structure}
//declaration of required Variables
bool J = false;
bool antiJ = false;
bool bothJ = false;


for(int idj1 = 0; idj1 < Dijunction1.size(); idj1++){ //searching dijunction
  for(int idj2 = 0; idj2 < Dijunction1.at(idj1).size(); idj2++){ //searching dijunction's leg
  //std::cout <<endl<<"let's check these particles "<<endl;
  //std::cout <<" { "<<Dijunction1.at(idj1).at(idj2).at(0).col()<<" , "<<Dijunction1.at(idj1).at(idj2).at(0).acol()<<" } "<<endl;
  //std::cout <<" { "<<Dijunction1.at(idj1).at(idj2).back().col()<<" , "<<Dijunction1.at(idj1).at(idj2).back().acol()<<" } "<<endl;

      for(int itpj1 = 0; itpj1 < IMStructure1.size(); itpj1++){ //searching Tempjunction
        for(int itpj2 = 1; itpj2 < Tempjunctions.at(itpj1).size(); itpj2++){ //searching color tags in tempjunction{ {_+1, 0 }, {_+1, tag1}, {_+1, tag2}, {_+1, tag3} }
          if(Dijunction1.at(idj1).at(idj2).at(0).col() == Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0)).at(itpj2).at(1) || //check whether first or last particle has correponding color tags for J
             Dijunction1.at(idj1).at(idj2).back().col() == Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0)).at(itpj2).at(1) ){
               J = true;
             }
          if(Dijunction1.at(idj1).at(idj2).at(0).acol() == Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0)).at(itpj2).at(1) || //check whether first or last particle has correponding color tags for J
             Dijunction1.at(idj1).at(idj2).back().acol() == Tempjunctions.at(IMStructure1.at(itpj1).at(0).at(0)).at(itpj2).at(1) ){
               antiJ = true;
            }

          }
        }
        if(J == true && antiJ == false){ // this leg is a part of junction
          DijunctionInfo2.push_back(1);
        }
        if(J == false && antiJ == true){ // this leg is a part of junction
          DijunctionInfo2.push_back(-1);
        }
        if ( J == true && antiJ == true){
          bothJ = true;
        }
        if( bothJ ){ //this leg is shared leg, so that fake parton shoud be added after all
          DijunctionInfo2.push_back(0);
        }
        J = false; //reseting variables to check other legs in Dijunction structure
        antiJ = false;
        bothJ = false;
      }
      DijunctionInfo1.push_back(DijunctionInfo2);
      DijunctionInfo2.clear();
  }//so far, vector about the information of dijunction structure is formed,{ so, we could put identity 1-leg first, and put identity-0 leg and -1 leg}
  //based on the information above, we need to add fake partons, which are used for telling pythia about the color flow{fake particle id is ignored, since we'll set
// negative status flag}, with this color flow information, PYTHIA can detect J and anti-J

//additional proceudre to change the configuration of shared leg by searching through dijunction1 and dijunction1info vectors
for(int ileg1 = 0 ; ileg1 < Dijunction1.size(); ileg1++){
 for(int ileg2 = 0; ileg2 < Dijunction1.at(ileg1).size(); ileg2++){
   bool needflip = false;
   if((DijunctionInfo1.at(ileg1).at(ileg2) == 0 && Dijunction1.at(ileg1).at(ileg2).size() != 1) &&
      (Dijunction1.at(ileg1).at(ileg2).at(0).col() == Dijunction1.at(ileg1).at(ileg2).at(1).acol()) ){ // In this case, we need to flip the order of the partons in this shared leg
      needflip = true;
     }
   if(needflip){
     std::reverse(Dijunction1.at(ileg1).at(ileg2).begin() , Dijunction1.at(ileg1).at(ileg2).end() );
   }
     needflip = false;
 }
}


  //dignostic measure{DijunctionInfo1}
  //std::cout <<endl;
  //std::cout <<" DijunctionInfo Checking"<<endl;
  for(int icheck1 = 0; icheck1 < DijunctionInfo1.size(); icheck1++){
    //std::cout <<" Let's Check all the elements! : ";
    for(int icheck2 = 0; icheck2 < DijunctionInfo1.at(icheck1).size(); icheck2++){
    //std::cout <<DijunctionInfo1.at(icheck1).at(icheck2)<<" , ";
    }
    //std::cout <<endl;
  }

  //dignostic measure{Dijunction1}
  //std::cout <<endl;
  //std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Dijunction1.size(); icheck1++){
    //std::cout <<" Dijunction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Dijunction1.at(icheck1).size(); icheck2++){
      //std::cout <<" Leg "<<icheck2<<" : "<<"{ identity : "<< DijunctionInfo1.at(icheck1).at(icheck2)<<" }  ";
      for(int icheck3 = 0; icheck3 < Dijunction1.at(icheck1).at(icheck2).size(); icheck3++){
        //std::cout <<" ( "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" , "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      //std::cout <<endl;
    }
  }
  //dignostic measure{Singlejunction1}
  //std::cout <<endl;
  //std::cout <<" List of Single Junction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Singlejunction1.size(); icheck1++){
    //std::cout <<" Single junction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Singlejunction1.at(icheck1).size(); icheck2++){
      //std::cout <<" Leg "<<icheck2<<" : ";
      for(int icheck3 = 0; icheck3 < Singlejunction1.at(icheck1).at(icheck2).size(); icheck3++){
        //std::cout <<" ( "<<Singlejunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" , "<<Singlejunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      //std::cout <<endl;
    }
  }

//For the validity with PYTHIA running, each Leg in Single or Dijunction Structure should have q or q-bar as endpoint particle, so checking whether the gluon is located at the endpoint of the Leg
//Starting from Leg's in single junction.

for(int is1 = 0; is1 < Singlejunction1.size(); is1++){
  HHparton p1 = Singlejunction1.at(is1).at(0).back(); //endpoint particle of first leg in is1_st Singlejunction1
  HHparton p2 = Singlejunction1.at(is1).at(1).back(); //endpoint particle of second leg in is1_st Singlejunction1
  HHparton p3 = Singlejunction1.at(is1).at(2).back();//endpoint particle of third leg in is1_st Singlejunction1

  // basically, assume the junction. but check the id's of them and change identity of the condition for antijunction is meeted
  bool ajunction = false;

  if( p1.id() < 0 || p2.id() < 0 || p3.id() < 0){
    ajunction = true;
  }

  for(int is2 = 0; is2 < Singlejunction1.at(is1).size(); is2++){
   HHparton endpoint = Singlejunction1.at(is1).at(is2).back();

   if(endpoint.col() != 0 && endpoint.acol() != 0 && ajunction){
     HHparton fakeaq;
     fakeaq.id(-1);
     fakeaq.col(0); fakeaq.acol(endpoint.col());
     Singlejunction1.at(is1).at(is2).push_back(fakeaq);
   }
  }

  for(int is3 = 0; is3 < Singlejunction1.at(is1).size(); is3++){
   HHparton endpoint = Singlejunction1.at(is1).at(is3).back();

   if(endpoint.col() != 0 && endpoint.acol() != 0 && !ajunction){
     HHparton fakeq;
     fakeq.id(1);
     fakeq.acol(0); fakeq.col(endpoint.acol());
     Singlejunction1.at(is1).at(is3).push_back(fakeq);
   }
  }
  ajunction = false; //resetting the variable.
}


//Working for Leg's in dijunction structure.
for(int idj1 = 0; idj1 < Dijunction1.size(); idj1++){
  for(int idj2 = 0; idj2 < Dijunction1.at(idj1).size(); idj2++){
    HHparton endpoint = Dijunction1.at(idj1).at(idj2).back();

    if(DijunctionInfo1.at(idj1).at(idj2) == -1 && endpoint.col() != 0 && endpoint.acol() != 0){
      HHparton fakeaq;
      fakeaq.id(-1);
      fakeaq.col(0); fakeaq.acol(endpoint.col());
      Dijunction1.at(idj1).at(idj2).push_back(fakeaq);
    }
  }

  for(int idj3 = 0; idj3 < Dijunction1.at(idj1).size(); idj3++){
    HHparton endpoint = Dijunction1.at(idj1).at(idj3).back();

    if(DijunctionInfo1.at(idj1).at(idj3) == 1 && endpoint.col() != 0 && endpoint.acol() != 0){
      HHparton fakeq;
      fakeq.id(1);
      fakeq.col(endpoint.acol()); fakeq.col(0);
      Dijunction1.at(idj1).at(idj3).push_back(fakeq);
    }
  }
}




//since the repeating is needed for completing the junction structure, It is inevitable that there are repeated particle in Recombeary1, It's compicated to explain, and useless to focus on.
// so that just delete repeated particles by check color Tags

 for(int irb1 = 0; irb1 < Recombearly1.size(); irb1++){ // pick one baryon
  for(int irb2 = 0; irb2 < Recombearly1.at(irb1).size(); irb2++){ // pick one particle in baryon
    for(int irb3 = 0; irb3 < Recombearly1.at(irb1).size(); irb3++){ // checking through all particles in baryon
      if( (Recombearly1.at(irb1).at(irb2).col() == Recombearly1.at(irb1).at(irb3).col()) &&
      Recombearly1.at(irb1).at(irb2).acol() == Recombearly1.at(irb1).at(irb3).acol()){
        Recombearly1.at(irb1).erase(Recombearly1.at(irb1).begin()+irb3);
      }
    }
  }
 }




  //dignositic measure{ealry recombined baryon}
  //std::cout <<endl;
  //std::cout <<" number of early recombined baryons : "<<Recombearly1.size()<<endl;
  for(int ibary1 = 0; ibary1 < Recombearly1.size(); ibary1++){
    //std::cout <<"the "<< ibary1 <<" st baryon is made of ";
    for(int ibary2 = 0; ibary2 < Recombearly1.at(ibary1).size(); ibary2++){
      //std::cout <<" { "<<Recombearly1.at(ibary1).at(ibary2).col()<<" , "<< Recombearly1.at(ibary1).at(ibary2).acol()  <<" } , ";
    }
    //std::cout <<endl;
  }




//now we need to form fake partons and baryons for appending particles to PYTHIA
//first, set loop for pre-recombined hadrons from Recombearly1
for(int irb = 0; irb < Recombearly1.size(); irb++){
  double hbac2 = hbarc*hbarc;

  FourVector Pbaryon;
  Pbaryon.Set(Recombearly1.at(irb)[0].px()+Recombearly1.at(irb)[1].px()+Recombearly1.at(irb)[2].px(),Recombearly1.at(irb)[0].py()+Recombearly1.at(irb)[1].py()+Recombearly1.at(irb)[2].py(),Recombearly1.at(irb)[0].pz()+Recombearly1.at(irb)[1].pz()+Recombearly1.at(irb)[2].pz(),0.);

  //baryon(CM) velocity
  FourVector betaB; //really p[i]/e below
  betaB.Set(Pbaryon.x()/(Recombearly1.at(irb)[0].e()+Recombearly1.at(irb)[1].e()+Recombearly1.at(irb)[2].e()),Pbaryon.y()/(Recombearly1.at(irb)[0].e()+Recombearly1.at(irb)[1].e()+Recombearly1.at(irb)[2].e()),Pbaryon.z()/(Recombearly1.at(irb)[0].e()+Recombearly1.at(irb)[1].e()+Recombearly1.at(irb)[2].e()),0.);
  betaB.Set(betaB.x(),betaB.y(),betaB.z(),1./(sqrt(1. - (betaB.x()*betaB.x() + betaB.y()*betaB.y() + betaB.z()*betaB.z()))));

  //boosting into CM frame
  FourVector pos_BCM[3], p_BCM[3];
  pos_BCM[0] = Recombearly1.at(irb)[0].boost_pos(betaB); pos_BCM[1] = Recombearly1.at(irb)[1].boost_pos(betaB); pos_BCM[2] = Recombearly1.at(irb)[2].boost_pos(betaB);
    p_BCM[0] = Recombearly1.at(irb)[0].boost_P(betaB);     p_BCM[1] = Recombearly1.at(irb)[1].boost_P(betaB);     p_BCM[2] = Recombearly1.at(irb)[2].boost_P(betaB);

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
  (cur_pos[0].x()*Recombearly1.at(irb)[0].mass()+cur_pos[1].x()*Recombearly1.at(irb)[1].mass()+cur_pos[2].x()*Recombearly1.at(irb)[2].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  (cur_pos[0].y()*Recombearly1.at(irb)[0].mass()+cur_pos[1].y()*Recombearly1.at(irb)[1].mass()+cur_pos[2].y()*Recombearly1.at(irb)[2].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  (cur_pos[0].z()*Recombearly1.at(irb)[0].mass()+cur_pos[1].z()*Recombearly1.at(irb)[1].mass()+cur_pos[2].z()*Recombearly1.at(irb)[2].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  curtime
  );

  //finding position of baryon in lab frame
  betaB.Set(-betaB.x(),-betaB.y(),-betaB.z(),betaB.t());
  FourVector pos_lab = HHboost(betaB, pos_CM);


  //finding relative momenta of partons in CM frame
  FourVector k_rel[2];
  k_rel[0].Set(
  (Recombearly1.at(irb)[1].mass()*p_BCM[0].x()-Recombearly1.at(irb)[0].mass()*p_BCM[1].x())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()),
  (Recombearly1.at(irb)[1].mass()*p_BCM[0].y()-Recombearly1.at(irb)[0].mass()*p_BCM[1].y())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()),
  (Recombearly1.at(irb)[1].mass()*p_BCM[0].z()-Recombearly1.at(irb)[0].mass()*p_BCM[1].z())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()),
  0.);
  k_rel[1].Set(
  (Recombearly1.at(irb)[2].mass()*(p_BCM[0].x()+p_BCM[1].x())-(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())*p_BCM[2].x())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  (Recombearly1.at(irb)[2].mass()*(p_BCM[0].y()+p_BCM[1].y())-(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())*p_BCM[2].y())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  (Recombearly1.at(irb)[2].mass()*(p_BCM[0].z()+p_BCM[1].z())-(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())*p_BCM[2].z())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass()+Recombearly1.at(irb)[2].mass()),
  0.);

  //finding relative positions of partons in CM frame
  FourVector pos_rel[2];
  pos_rel[0].Set((cur_pos[0].x()-cur_pos[1].x())/sqrt(2.),(cur_pos[0].y()-cur_pos[1].y())/sqrt(2.),(cur_pos[0].z()-cur_pos[1].z())/sqrt(2.),0.);
  pos_rel[1].Set(
  ((cur_pos[0].x()*Recombearly1.at(irb)[0].mass()+cur_pos[1].x()*Recombearly1.at(irb)[1].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())-cur_pos[2].x())*sqrt(2./3.),
  ((cur_pos[0].y()*Recombearly1.at(irb)[0].mass()+cur_pos[1].y()*Recombearly1.at(irb)[1].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())-cur_pos[2].y())*sqrt(2./3.),
  ((cur_pos[0].z()*Recombearly1.at(irb)[0].mass()+cur_pos[1].z()*Recombearly1.at(irb)[1].mass())/(Recombearly1.at(irb)[0].mass()+Recombearly1.at(irb)[1].mass())-cur_pos[2].z())*sqrt(2./3.),
  0.);

  //so far, values to form a baryon are set, now make hadron objects based on this information
  //now we're forming the hadron
  HHhadron formedhadron;
  //setting the hadron values: is a recombined hadron, mass, and parents
  formedhadron.is_recohad(true); formedhadron.mass( p_BCM[0].t() + p_BCM[1].t() + p_BCM[2].t() );

  //setting hadron position and momentum vectors
  Pbaryon.Set(Pbaryon.x(),Pbaryon.y(),Pbaryon.z(),sqrt(Pbaryon.x()*Pbaryon.x() + Pbaryon.y()*Pbaryon.y() + Pbaryon.z()*Pbaryon.z() + formedhadron.mass()*formedhadron.mass()));
  formedhadron.pos(pos_lab); formedhadron.P(Pbaryon);

  //need to choose *what* hadron we've formed... base this on the parton id's, mass, & if excited
  //might want to do this differently? void f'n(partoncollection, formedhadron)?
  //since set_baryon_id function works with object of parton_colection class, make temporary parton_collection
  parton_collection tempobject; //since the generator resets all components, don't need to reset the vactors in it
  tempobject.add(Recombearly1.at(irb)[0]);
  tempobject.add(Recombearly1.at(irb)[1]);
  tempobject.add(Recombearly1.at(irb)[2]);

  set_baryon_id(tempobject, formedhadron);

  //need to add the hadron to the collection
  HH_hadrons.add(formedhadron);
  //std::cout <<endl<<"finally Early recombined baryon with id : "<< formedhadron.id()  <<" is added to the system!"<<endl;

}


//from now we need re indexing a of the particles in SP_remnants by comparing color tag in Juncstructure component.
//Maybe there would not be a problem in the used_junction value, but take most defensive measure for the possible error.
std::vector<int> checking;
bool indicecheck = false;
for (int irem = 0; irem < SP_remnants.num(); irem++){
  for(int ijs1 = 0; ijs1 < JuncStructure.size(); ijs1++){
    for(int ijs2 = 0; ijs2 < JuncStructure.at(ijs1).size(); ijs2++){
      for(int ijs3 = 0; ijs3 < JuncStructure.at(ijs1).at(ijs2).size(); ijs3++){
        if( (SP_remnants[irem].col() == JuncStructure.at(ijs1).at(ijs2).at(ijs3).col()) &&
        (SP_remnants[irem].acol() == JuncStructure.at(ijs1).at(ijs2).at(ijs3).acol()) ){
          SP_remnants[irem].used_junction(true);
          checking.push_back(irem);
        }
      }
    }
  }
}
for (int irem = 0; irem < SP_remnants.num(); irem++){
  for(int icheck = 0; icheck < checking.size(); icheck++){
    if(irem == checking.at(icheck)){ indicecheck = true; }
  }
  if(indicecheck = false){
    SP_remnants[irem].used_junction(false);
  }
}




// TODO: many things to be cleared for the next running, remember that!!

// Now set up for junction is done, let's check the left partons to establish string structure


//std::cout <<endl<<"let's check left particles!"<<endl;
 for(int irp = 0 ; irp < SP_remnants.num(); irp++) {
   if(SP_remnants[irp].used_junction() == false){
     //std::cout << "( "<<SP_remnants[irp].pid()<<" , "<<SP_remnants[irp].col()<<" , "<<SP_remnants[irp].acol()<<" ) "<<endl;
   }
 }

 // for the string repair among leftquarks and tailoredstring, check whether all color tags exist in pair

for(int ileft = 0; ileft < SP_remnants.num(); ileft++){
  if(SP_remnants[ileft].used_junction() == false){
    finalstring.push_back(SP_remnants[ileft]);
  }
}

for(int itail1 = 0; itail1 < Tailoredstring1.size(); itail1++ ){
  for(int itail2 = 0; itail2 < Tailoredstring1.at(itail1).size(); itail2++){
    finalstring.push_back(Tailoredstring1.at(itail1).at(itail2));
  }
}

//std::cout <<endl<<"check final particles!"<<endl;
  for(int ifin = 0; ifin < finalstring.size(); ifin++){
    //std::cout <<" { "<< finalstring.at(ifin).col() << " , " << finalstring.at(ifin).acol() << " } "<<endl;
  }

//to check the pairty of color tags, declare some vectors for set up

  vector<int> cols;
  vector<int> acols;

// first, set up col tags vector.
for(int icol = 0; icol < finalstring.size(); icol++){
  if(finalstring.at(icol).col() != 0){
    cols.push_back(finalstring.at(icol).col());
  }
}

for(int iacol = 0; iacol < finalstring.size(); iacol++){
  if(finalstring.at(iacol).acol() != 0){
    acols.push_back(finalstring.at(iacol).acol());
  }
}

//based on this information. check cols and generate required fakeparton
for(int ifq1 = 0; ifq1 < cols.size(); ifq1++){
  bool paired = false;
  for(int ifq2 = 0; ifq2 < acols.size(); ifq2++){
    if( cols.at(ifq1) == acols.at(ifq2)){
      paired = true;
    }
  }
  if( !paired ){
    HHparton fakeq;
    fakeq.id(-1);//set fake anti up quark
    fakeq.set_color(0);
    fakeq.set_anti_color(cols.at(ifq1));
    finalstring.push_back(fakeq);
  }
  paired = false;
}

for(int ifaq1 = 0; ifaq1 < acols.size(); ifaq1++){
  bool paired = false;
  for(int ifaq2 = 0; ifaq2 < cols.size(); ifaq2++){
    if( acols.at(ifaq1) == cols.at(ifaq2)){
      paired = true;
    }
  }
  if( !paired ){
    HHparton fakeaq;
    fakeaq.id(1);//add fake up quark id
    fakeaq.set_color(acols.at(ifaq1));
    fakeaq.set_anti_color(0);
    finalstring.push_back(fakeaq);
  }
  paired = false;
}
//reseting cols and acols vector for next running
cols.clear();
acols.clear();


//std::cout <<endl<<"check final particles!"<<endl;
  for(int ifin = 0; ifin < finalstring.size(); ifin++){
    //std::cout <<" { "<< finalstring.at(ifin).col() << " , " << finalstring.at(ifin).acol() << " } "<<endl;
  }

//so far, particles in finalstring are ready for being transfered into PYTHIA, since we don't need to set mother, daughter tag.
//but, still need to care about single junction and dijunction system, especially for the mother daughter tag!
// mother, daughter tags are working with the execution order{= order of being declared in pythia} so before tossing the datum to PYTHIA, rearrange order in the vector
// declare the transit space for these particles, while moving onto vector to vector, partons are assigned mother, daughter tag tor PYTHIA running{definition of that tags are given well in main21.cc file in   /installed pythia folder/bin/example/ }

  vector<vector<HHparton>> Transitdijunction1;
  vector<HHparton> Transitdijunction2; // two fake mother B, B-bar and fake q,qbar added and give mother daughter tags!
  vector<vector<vector<HHparton>>> Tempsorting1; // vectors for rearranging legs in dijunction{Desired arrangement : identity 1, 1, 0 , -1, -1}
  vector<vector<HHparton>> Tempsorting2;
  vector<vector<HHparton>> Transitsinglejunction1;
  vector<HHparton> Transitsinglejunction2;

  vector<HHparton> WaitingLineforPY; // final list of partons with complete M,D tag and col tags




  //set the value for incrementation in mother daughter tag.


  //and start from dijunction, first, check the identity{1: junction leg, -1: anti_junction leg, 0:shared leg}, since we'll read col tag first, consider identity=1 leg first and then 0, -1
  for(int dijuncfin1 = 0; dijuncfin1 < Dijunction1.size(); dijuncfin1++){ // Going to find shared leg and attach fake partons to the begin and the end of the leg
    for(int dijuncfin2= 0; dijuncfin2 < 5; dijuncfin2++){ //inspect through five legs
      if(DijunctionInfo1.at(dijuncfin1).at(dijuncfin2) == 0){ //check whether it's shared leg{identity =0} and add fakepartons to the first and second posittion in the vector
        HHparton fakeq;
        HHparton fakeqbar;
        fakeq.id(1); fakeq.set_color(Dijunction1.at(dijuncfin1).at(dijuncfin2).at(0).col());
        fakeq.PY_stat(-21); fakeq.mass(xmq); fakeq.e(xmq); fakeq.orig(-1);
        int endpoint = Dijunction1.at(dijuncfin1).at(dijuncfin2).size() - 1;
        fakeqbar.id(-1); fakeqbar.set_anti_color(Dijunction1.at(dijuncfin1).at(dijuncfin2).at(endpoint).acol());
        fakeqbar.PY_stat(-21); fakeqbar.mass(xmq); fakeqbar.e(xmq); fakeqbar.orig(-1);
        std::vector<HHparton>::iterator it2 = Dijunction1.at(dijuncfin1).at(dijuncfin2).begin();
        Dijunction1.at(dijuncfin1).at(dijuncfin2).insert(it2, fakeqbar);
        std::vector<HHparton>::iterator it1 = Dijunction1.at(dijuncfin1).at(dijuncfin2).begin();
        Dijunction1.at(dijuncfin1).at(dijuncfin2).insert(it1, fakeq);
      }
    }
    for(int isort1 = 0; isort1 < 5; isort1++){
      if(DijunctionInfo1.at(dijuncfin1).at(isort1) == 1){
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    for(int isort1 = 0; isort1 < 5; isort1++){
      if(DijunctionInfo1.at(dijuncfin1).at(isort1) == 0){
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    for(int isort1 = 0; isort1 < 5; isort1++){
      if(DijunctionInfo1.at(dijuncfin1).at(isort1) == -1){
        Tempsorting2.push_back(Dijunction1.at(dijuncfin1).at(isort1));
      }
    }
    Tempsorting1.push_back(Tempsorting2);
    Tempsorting2.clear();

  } //first looping to add fake partons is Finished, now sort the order of the legs to be identity 1,1,0,-1,-1

  //dignostic measure{Dijunction1}
  //std::cout <<endl;
  //std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Dijunction1.size(); icheck1++){
    //std::cout <<" Dijunction : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Dijunction1.at(icheck1).size(); icheck2++){
      //std::cout <<" Leg "<<icheck2<<" : "<<"{ identity : "<< DijunctionInfo1.at(icheck1).at(icheck2)<<" }  ";
      for(int icheck3 = 0; icheck3 < Dijunction1.at(icheck1).at(icheck2).size(); icheck3++){
        //std::cout <<" ( "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).col()<<" , "<<Dijunction1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      //std::cout <<endl;
    }
  }
  //dignostic measure{Tempsorting}
  //std::cout <<endl;
  //std::cout <<" List of Dijunction Structure"<<endl;
  for(int icheck1 = 0; icheck1 < Tempsorting1.size(); icheck1++){
    //std::cout <<" Tempsorted Dijunction1 : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Tempsorting1.at(icheck1).size(); icheck2++){
      for(int icheck3 = 0; icheck3 < Tempsorting1.at(icheck1).at(icheck2).size(); icheck3++){
        //std::cout <<" ( "<<Tempsorting1.at(icheck1).at(icheck2).at(icheck3).col()<<" , "<<Tempsorting1.at(icheck1).at(icheck2).at(icheck3).acol()<< " )  ";
      }
      //std::cout <<endl;
    }
  }
// adding fake partons and sorting them with right order{Junctionleg1,2 shared leg, antijunctionleg1,2} is finished so far, Now form a fake mothers to prepare for PYTHIA
 for(int itag1 = 0; itag1 < Tempsorting1.size(); itag1++){ /// we need to link Tempsorting1 with Transitdijunction1,
//First, put fake B,Bbar into 1st and 2nd posittion of Transitdijunction2

     parton_collection FakeBaryonElements; // array of quarks in fake B
     parton_collection FakeAntibaryonElements; // array of anti quarks in fake Bbar

     FakeBaryonElements.add(Tempsorting1.at(itag1).at(0).back()); // quarks to form fake FakeBaryon
     FakeBaryonElements.add(Tempsorting1.at(itag1).at(1).back());
     FakeBaryonElements.add(Tempsorting1.at(itag1).at(2).at(0));
// for the corrspondence with pdg particle id, sorting the these ids based on absolute value.

    for(int iswap = 0; iswap < 3; iswap++){
    if(abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[1].id())){ std::swap(FakeBaryonElements[2], FakeBaryonElements[1]) ; }
    if(abs(FakeBaryonElements[1].id()) > abs(FakeBaryonElements[0].id())){ std::swap(FakeBaryonElements[1], FakeBaryonElements[0]) ; }
    if(abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[0].id())){ std::swap(FakeBaryonElements[2], FakeBaryonElements[0]) ; }
    }


     FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(2).at(1)); //anti quarks to form FakeAntiBaryon
     FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(3).back());
     FakeAntibaryonElements.add(Tempsorting1.at(itag1).at(4).back());

     for(int iswap = 0; iswap < 3; iswap++){
     if(abs(FakeAntibaryonElements[2].id()) > abs(FakeAntibaryonElements[1].id())){ std::swap(FakeAntibaryonElements[2], FakeAntibaryonElements[1]) ; }
     if(abs(FakeAntibaryonElements[1].id()) > abs(FakeAntibaryonElements[0].id())){ std::swap(FakeAntibaryonElements[1], FakeAntibaryonElements[0]) ; }
     if(abs(FakeAntibaryonElements[2].id()) > abs(FakeAntibaryonElements[0].id())){ std::swap(FakeAntibaryonElements[2], FakeAntibaryonElements[0]) ; }
     }

     int m_id1 = abs(FakeBaryonElements[0].id())*1000 + abs(FakeBaryonElements[1].id())*100 + abs(FakeBaryonElements[2].id())*10 + 4;
     int m_id2 = (abs(FakeAntibaryonElements[0].id())*1000 + abs(FakeAntibaryonElements[0].id())*100 + abs(FakeAntibaryonElements[0].id())*10 + 4)*-1;

     FakeBaryonElements.clear();
     FakeAntibaryonElements.clear();
     //std::cout <<endl<<" Trial : "<<itag1<<endl;
     //std::cout <<endl<<"Fake Baryon id is "<<m_id1<<endl;
     //std::cout <<"Fake AntiBaryon id is "<<m_id2<<endl;

     HHparton fakeB;
     HHparton fakeBbar;
     fakeB.id(m_id1); fakeB.PY_stat(-11); fakeB.mass(2*xmq); fakeB.e(2*xmq);
     fakeBbar.id(m_id2); fakeBbar.PY_stat(-11); fakeBbar.mass(2*xmq); fakeBbar.e(2*xmq);

     //so far, we formed two fake partons, which would be at the center of Junction and Antijunction, Now put them 0th and 1st position of Transirdijunction1 vector

     Transitdijunction2.push_back(fakeB);
     Transitdijunction2.push_back(fakeBbar);
     Transitdijunction1.push_back(Transitdijunction2);
     Transitdijunction2.clear();

 }
 //the next step is to set M,D Tags
 for(int iMD1 = 0; iMD1 < Tempsorting1.size(); iMD1++){ //after working in these Legs, we will return the partons to Transitdijunction1,
// so far, the order in the vectors{Dijunction1, DijunctionInfo1, Tempsorting, Transitdijunction1} are unified, but the difference is the value in each vectors
  int Intag = 0; // tag for internal M,D tagging process
     vector<HHparton> Leg1 = Tempsorting1.at(iMD1).at(0); // for convenience, reassigning legs in sorted dijunction structure.{identity 1,1,0,-1,-1}
     vector<HHparton> Leg2 = Tempsorting1.at(iMD1).at(1);
     vector<HHparton> Leg3 = Tempsorting1.at(iMD1).at(2);
     vector<HHparton> Leg4 = Tempsorting1.at(iMD1).at(3);
     vector<HHparton> Leg5 = Tempsorting1.at(iMD1).at(4);
//starting from Leg1{identity1}
  for(int ileg1 = 0; ileg1 < Leg1.size(); ileg1++){
    Leg1.at(ileg1).PY_par1(0);  // the mother of first leg1, 2 should be first fakebaryon{located 1st in the Transitdijunction2 vector!}
    Leg1.at(ileg1).PY_par2(0);
    Transitdijunction1.at(iMD1).push_back(Leg1.at(ileg1)); //toss Leg1 to the Transirdijunction1.at{iMD}, which would be the last step before WaitingLineforPY
  }
  for(int ileg2 = 0; ileg2 < Leg2.size(); ileg2++){
    Leg2.at(ileg2).PY_par1(0); // the mother of first leg1, 2 should be first fakebaryon{located 1st in the Transitdijunction2 vector!}
    Leg2.at(ileg2).PY_par2(0);
    Transitdijunction1.at(iMD1).push_back(Leg2.at(ileg2));
  }
//for leg3, because of first and end fake particles, it goes differently from other legs
  Intag = Intag + 2 + Leg1.size() + Leg2.size(); // Previous tag + two fakemothers + Leg1 + Leg2
  //First of all, fake parton pairs are located for the convenience of MD tagging
  Leg3.at(0).PY_par1(0);//daughter of fakeB
  Leg3.at(0).PY_par2(0);
  Leg3.at(0).PY_dau1(Intag + Leg4.size() + Leg5.size() + 2); // daughter tags for gluons between two fake qqbar
  Leg3.at(0).PY_dau2(Intag + Leg4.size() + Leg5.size() + Leg3.size()- 1); // since
  Leg3.at(1).PY_par1(1);//daughter of fakeBbar
  Leg3.at(1).PY_par2(1);
  Leg3.at(1).PY_dau1(Intag + Leg4.size() + Leg5.size() + 2);
  Leg3.at(1).PY_dau2(Intag + Leg4.size() + Leg5.size() + Leg3.size()- 1);

  Transitdijunction1.at(iMD1).push_back(Leg3.at(0));
  Transitdijunction1.at(iMD1).push_back(Leg3.at(1)); // push first two fake partons and append remnants at the last, this is for convenience


  Transitdijunction1.at(iMD1).at(0).PY_dau1(2); // correcting daughter tags of fake B located at 0
  Transitdijunction1.at(iMD1).at(0).PY_dau2(Intag);
  Transitdijunction1.at(iMD1).at(1).PY_dau1(Intag + 1); //setting starting daughter tag of fakeBbar located at 1
  Transitdijunction1.at(iMD1).at(1).PY_dau2(Intag + 1 + Leg4.size() + Leg5.size() );

  for(int ileg4 = 0 ; ileg4 < Leg4.size(); ileg4++){
    Leg4.at(ileg4).PY_par1(1);
    Leg4.at(ileg4).PY_par2(1);
    Transitdijunction1.at(iMD1).push_back(Leg4.at(ileg4));
  }
  for(int ileg5 = 0 ; ileg5 < Leg5.size(); ileg5++){
    Leg5.at(ileg5).PY_par1(1);
    Leg5.at(ileg5).PY_par2(1);
    Transitdijunction1.at(iMD1).push_back(Leg5.at(ileg5));
  }

  for(int ileg3 = 2; ileg3 < Leg3.size(); ileg3++){ // mother tags for the gluons between two fake qqbar pair
  Leg3.at(ileg3).PY_par1(Intag);
  Leg3.at(ileg3).PY_par2(Intag + 1);
  Transitdijunction1.at(iMD1).push_back(Leg3.at(ileg3));
  }

  Intag = 0;//reseting the tag

} // Inner tagging loop is finished for dijunction system,{That is, Information in Tempsorting1 is transfered to Transitdijunction1}

//dignostic measure{Transitdijunction1}
//std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2 >> dau1 >> dau2 "<< endl;
for(int icheck1 = 0; icheck1 <  Transitdijunction1.size(); icheck1++ ){
  for(int i = 0; i <  Transitdijunction1.at(icheck1).size(); i++){
    vector<HHparton> temp = Transitdijunction1.at(icheck1);
    //std::cout  << i <<"   "<< temp.at(i).id() <<"   "<< temp.at(i).PY_stat() <<"   "<< temp.at(i).PY_par1() <<"   "<< temp.at(i).PY_par2() <<"   "
    //<< temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<  temp.at(i).col()<< " , " << temp.at(i).acol() << " ) " <<endl;
  }
}

//TODO: Now,start woring on Single Junction System.
for(int iSJ = 0; iSJ < Singlejunction1.size(); iSJ++){ // this process is for trasferring information from Singlejunction1 into Transitsinglejunction1, and finally they will be located at WaitingLineforPY

  //based on the partons at the end of each Leg, add fake mother baryon first.
   parton_collection FakeBaryonElements; // array of quarks in fake mother

   //HHparton q1 = Singlejunction1.at(iSJ).at(0).at(0);
   //HHparton q2 = Singlejunction1.at(iSJ).at(1).at(0);
   //HHparton q3 = Singlejunction1.at(iSJ).at(2).at(0);

   HHparton q1 = Singlejunction1.at(iSJ).at(0).back();
   HHparton q2 = Singlejunction1.at(iSJ).at(1).back();
   HHparton q3 = Singlejunction1.at(iSJ).at(2).back();

   HHparton q1fin = q1;
   HHparton q2fin = q2;
   HHparton q3fin = q3;

   FakeBaryonElements.add(q1fin); //add first particle in 1st leg in iSJ_st singlejunction
   FakeBaryonElements.add(q2fin); //add first particle in 2nd leg in iSJ_st singlejunction
   FakeBaryonElements.add(q3fin); //add first particle in 3rd leg in iSJ_st singlejunction

   // for the corrspondence with pdg particle id, sorting the these ids based on absolute value

   for(int iswap = 0; iswap < 3; iswap++){
   if(abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[1].id())){ std::swap(FakeBaryonElements[2], FakeBaryonElements[1]) ; }
   if(abs(FakeBaryonElements[1].id()) > abs(FakeBaryonElements[0].id())){ std::swap(FakeBaryonElements[1], FakeBaryonElements[0]) ; }
   if(abs(FakeBaryonElements[2].id()) > abs(FakeBaryonElements[0].id())){ std::swap(FakeBaryonElements[2], FakeBaryonElements[0]) ; }
   }



   int m_id = abs(FakeBaryonElements[0].id())*1000 + abs(FakeBaryonElements[1].id())*100 + abs(FakeBaryonElements[2].id())*10 + 4;
   if(q1.id() < 0){
     m_id = -1*m_id;
   }

   HHparton fakeB;
   fakeB.id(m_id); fakeB.PY_stat(-11); fakeB.mass(2*xmq); fakeB.e(2*xmq);
   Transitsinglejunction2.push_back(fakeB);
   Transitsinglejunction1.push_back(Transitsinglejunction2);
   Transitsinglejunction2.clear();
} // so farfake mother is added to all single junction vectors{vectors of partons}, so we need to remember all execution number is incremented by one because of this fake baryons

for(int iSJ = 0; iSJ < Singlejunction1.size(); iSJ++){

  vector<HHparton> Leg1 = Singlejunction1.at(iSJ).at(0);
  vector<HHparton> Leg2 = Singlejunction1.at(iSJ).at(1);
  vector<HHparton> Leg3 = Singlejunction1.at(iSJ).at(2);

  for(int ileg1 = 0; ileg1 < Leg1.size(); ileg1++ ){
    Leg1.at(ileg1).PY_par1(0); // setting the mother tag as zero, which indicates the first particle in the Transit
    Leg1.at(ileg1).PY_par2(0);
    Transitsinglejunction1.at(iSJ).push_back(Leg1.at(ileg1));
  }

  for(int ileg2 = 0; ileg2 < Leg2.size(); ileg2++ ){
    Leg2.at(ileg2).PY_par1(0); // setting the mother tag as zero, which indicates the first particle in the Transit
    Leg2.at(ileg2).PY_par2(0);
    Transitsinglejunction1.at(iSJ).push_back(Leg2.at(ileg2));
  }

  for(int ileg3 = 0; ileg3 < Leg3.size(); ileg3++ ){
    Leg3.at(ileg3).PY_par1(0); // setting the mother tag as zero, which indicates the first particle in the Transit
    Leg3.at(ileg3).PY_par2(0);
    Transitsinglejunction1.at(iSJ).push_back(Leg3.at(ileg3));
  }

  Transitsinglejunction1.at(iSJ).at(0).PY_dau1(1);//dau tag starting from 1 {after zero, which means fake mother baryon}
  Transitsinglejunction1.at(iSJ).at(0).PY_dau2(Leg1.size() + Leg2.size() + Leg3.size());//dau tag starting from 1 {after zero, which means fake mother baryon}
}// At last, iSJ_st Singlejunction1 partons are moved to  iSJ_st vector component in Transitsinglejunction1

//dignostic measure{Transitsinglejunction1}
//std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2 >> dau1 >> dau2 "<< endl;
for(int icheck1 = 0; icheck1 <  Transitsinglejunction1.size(); icheck1++ ){
  //std::cout <<endl<<"SingleJunction : "<<icheck1<<endl;
  for(int i = 0; i <  Transitsinglejunction1.at(icheck1).size(); i++){
    vector<HHparton> temp = Transitsinglejunction1.at(icheck1);
    //std::cout  << i <<"   "<< temp.at(i).id() <<"   "<< temp.at(i).PY_stat() <<"   "<< temp.at(i).PY_par1() <<"   "<< temp.at(i).PY_par2() <<"   "
    //<< temp.at(i).PY_dau1() <<"   "<< temp.at(i).PY_dau2() << "  ( "<<  temp.at(i).col()<< " , " << temp.at(i).acol() << " ) " <<endl;
  }
}


//

  /*
  //dignostic measure{Dijunction1}
  std::cout <<endl;
  std::cout <<" List of Tempsorting for Dijunction"<<endl;
  for(int icheck1 = 0; icheck1 < Tempsorting1.size(); icheck1++){
    std::cout <<" Sorted Leg : "<<icheck1<<endl;
    for(int icheck2 = 0 ; icheck2 < Tempsorting1.at(icheck1).size(); icheck2++){
        std::cout <<" ( "<<Tempsorting1.at(icheck1).at(icheck2).col()<<" , "<<Tempsorting1.at(icheck1).at(icheck2).acol()<< " )  ";
    }
    std::cout <<endl;
  }
  */

//So far, M,D tags for each junction and dijunction structure are designated. and EP_conservation checking will be done through Transirdijunction1 and Transitsinglejunction1 vectors
for(int itd1 = 0; itd1 < Transitdijunction1.size(); itd1++){
    bool EP_conserved = false;
    while(!EP_conserved){
      EP_conserved = true;
      for(int i=0; i<Transitdijunction1.at(itd1).size(); ++i){
        if(Transitdijunction1[itd1][i].PY_stat() >= 0){continue;}
        FourVector P_new(0.,0.,0.,0.); FourVector pos_new(0.,0.,0.,0.);
        int jmax = (Transitdijunction1[itd1][i].PY_dau2() > Transitdijunction1[itd1][i].PY_dau1()) ? Transitdijunction1[itd1][i].PY_dau2() : Transitdijunction1[itd1][i].PY_dau1();
        for(int j=Transitdijunction1[itd1][i].PY_dau1(); j<jmax+1; ++j){
          double n = double(j-Transitdijunction1[itd1][i].PY_dau1())+1.;
          P_new.Set(P_new.x()+Transitdijunction1[itd1][j].px(),P_new.y()+Transitdijunction1[itd1][j].py(),P_new.z()+Transitdijunction1[itd1][j].pz(),P_new.t()+Transitdijunction1[itd1][j].e());
          pos_new.Set(pos_new.x()+(Transitdijunction1[itd1][j].x()-pos_new.x())/n,pos_new.y()+(Transitdijunction1[itd1][j].y()-pos_new.y())/n,pos_new.z()+(Transitdijunction1[itd1][j].z()-pos_new.z())/n,pos_new.t()+(Transitdijunction1[itd1][j].x_t()-pos_new.t())/n);
        }
        if((dif2(P_new,Transitdijunction1[itd1][i].P())+(P_new.t()-Transitdijunction1[itd1][i].e())*(P_new.t()-Transitdijunction1[itd1][i].e()) > 0.00000001/*0.0001^2*/) ||
          (dif2(pos_new,Transitdijunction1[itd1][i].pos())+(pos_new.t()-Transitdijunction1[itd1][i].x_t())*(pos_new.t()-Transitdijunction1[itd1][i].x_t()) > 0.00000001)){
          Transitdijunction1[itd1][i].P(P_new); Transitdijunction1[itd1][i].pos(pos_new);
          Transitdijunction1[itd1][i].mass( Transitdijunction1[itd1][i].e()*Transitdijunction1[itd1][i].e()
            - Transitdijunction1[itd1][i].px()*Transitdijunction1[itd1][i].px() - Transitdijunction1[itd1][i].py()*Transitdijunction1[itd1][i].py() - Transitdijunction1[itd1][i].pz()*Transitdijunction1[itd1][i].pz() );
          Transitdijunction1[itd1][i].mass( (Transitdijunction1[itd1][i].mass() >= 0.) ? sqrt(Transitdijunction1[itd1][i].mass()) : -sqrt(-Transitdijunction1[itd1][i].mass()) ); //don't need mass >0 for reco now
          EP_conserved = false;
        }
      }
    }
    for(int i=0; i<Transitdijunction1.at(itd1).size(); ++i){if((Transitdijunction1[itd1][i].PY_stat() == -21) && (Transitdijunction1[itd1][i].PY_dau2() > Transitdijunction1[itd1][i].PY_dau1())){
      Transitdijunction1[itd1][i].px(Transitdijunction1[itd1][i].px()/2.);
      Transitdijunction1[itd1][i].py(Transitdijunction1[itd1][i].py()/2.);
      Transitdijunction1[itd1][i].pz(Transitdijunction1[itd1][i].pz()/2.);
      Transitdijunction1[itd1][i].e( Transitdijunction1[itd1][i].e() /2.);
      Transitdijunction1[itd1][i].mass(Transitdijunction1[itd1][i].mass()/2.);
    }}
}

//for transitsinglejunction1
for(int its1 = 0; its1 < Transitsinglejunction1.size(); its1++){
    bool EP_conserved = false;
    while(!EP_conserved){
      EP_conserved = true;
      for(int i=0; i<Transitsinglejunction1.at(its1).size(); ++i){
        if(Transitsinglejunction1[its1][i].PY_stat() >= 0){continue;}
        FourVector P_new(0.,0.,0.,0.); FourVector pos_new(0.,0.,0.,0.);
        int jmax = (Transitsinglejunction1[its1][i].PY_dau2() > Transitsinglejunction1[its1][i].PY_dau1()) ? Transitsinglejunction1[its1][i].PY_dau2() : Transitsinglejunction1[its1][i].PY_dau1();
        for(int j=Transitsinglejunction1[its1][i].PY_dau1(); j<jmax+1; ++j){
          double n = double(j-Transitsinglejunction1[its1][i].PY_dau1())+1.;
          P_new.Set(P_new.x()+Transitsinglejunction1[its1][j].px(),P_new.y()+Transitsinglejunction1[its1][j].py(),P_new.z()+Transitsinglejunction1[its1][j].pz(),P_new.t()+Transitsinglejunction1[its1][j].e());
          pos_new.Set(pos_new.x()+(Transitsinglejunction1[its1][j].x()-pos_new.x())/n,pos_new.y()+(Transitsinglejunction1[its1][j].y()-pos_new.y())/n,pos_new.z()+(Transitsinglejunction1[its1][j].z()-pos_new.z())/n,pos_new.t()+(Transitsinglejunction1[its1][j].x_t()-pos_new.t())/n);
        }
        if((dif2(P_new,Transitsinglejunction1[its1][i].P())+(P_new.t()-Transitsinglejunction1[its1][i].e())*(P_new.t()-Transitsinglejunction1[its1][i].e()) > 0.00000001/*0.0001^2*/) ||
          (dif2(pos_new,Transitsinglejunction1[its1][i].pos())+(pos_new.t()-Transitsinglejunction1[its1][i].x_t())*(pos_new.t()-Transitsinglejunction1[its1][i].x_t()) > 0.00000001)){
          Transitsinglejunction1[its1][i].P(P_new); Transitsinglejunction1[its1][i].pos(pos_new);
          Transitsinglejunction1[its1][i].mass( Transitsinglejunction1[its1][i].e()*Transitsinglejunction1[its1][i].e()
            - Transitsinglejunction1[its1][i].px()*Transitsinglejunction1[its1][i].px() - Transitsinglejunction1[its1][i].py()*Transitsinglejunction1[its1][i].py() - Transitsinglejunction1[its1][i].pz()*Transitsinglejunction1[its1][i].pz() );
          Transitsinglejunction1[its1][i].mass( (Transitsinglejunction1[its1][i].mass() >= 0.) ? sqrt(Transitsinglejunction1[its1][i].mass()) : -sqrt(-Transitsinglejunction1[its1][i].mass()) ); //don't need mass >0 for reco now
          EP_conserved = false;
        }
      }
    }
    for(int i=0; i<Transitsinglejunction1.at(its1).size(); ++i){if((Transitsinglejunction1[its1][i].PY_stat() == -21) && (Transitsinglejunction1[its1][i].PY_dau2() > Transitsinglejunction1[its1][i].PY_dau1())){
      Transitsinglejunction1[its1][i].px(Transitsinglejunction1[its1][i].px()/2.);
      Transitsinglejunction1[its1][i].py(Transitsinglejunction1[its1][i].py()/2.);
      Transitsinglejunction1[its1][i].pz(Transitsinglejunction1[its1][i].pz()/2.);
      Transitsinglejunction1[its1][i].e( Transitsinglejunction1[its1][i].e() /2.);
      Transitsinglejunction1[its1][i].mass(Transitsinglejunction1[its1][i].mass()/2.);
    }}
}
//EP_conservation checking through dijunction and single junction is finished!!, and Mother Daughter tag incrementation done from below.{starting from Transitdijunction1 to Transitsinglejunction1}
int standard = 0; //integer for tracking number of particles attached.


for(int itagdij1 = 0; itagdij1 < Transitdijunction1.size(); itagdij1++ ){
  for(int itagdij2 = 0; itagdij2 < Transitdijunction1.at(itagdij1).size(); itagdij2++){
  if(Transitdijunction1.at(itagdij1).at(itagdij2).PY_par1() >= 0){
  Transitdijunction1.at(itagdij1).at(itagdij2).PY_par1(Transitdijunction1.at(itagdij1).at(itagdij2).PY_par1() + standard);
  }
  if(Transitdijunction1.at(itagdij1).at(itagdij2).PY_par2() >= 0){
  Transitdijunction1.at(itagdij1).at(itagdij2).PY_par2(Transitdijunction1.at(itagdij1).at(itagdij2).PY_par2() + standard);
  }
  if(Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau1() >= 0){
  Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau1(Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau1() + standard);
  }
  if(Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau2() >= 0){
  Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau2(Transitdijunction1.at(itagdij1).at(itagdij2).PY_dau2() + standard);
  }

  WaitingLineforPY.push_back(Transitdijunction1.at(itagdij1).at(itagdij2));
  }
  standard = standard + Transitdijunction1.at(itagdij1).size();
}

for(int itagsj1 = 0; itagsj1 < Transitsinglejunction1.size(); itagsj1++ ){
  for(int itagsj2 = 0; itagsj2 < Transitsinglejunction1.at(itagsj1).size(); itagsj2++){
  if(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par1() >= 0){
  Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par1(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par1() + standard);
  }
  if(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par2() >= 0){
  Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par2(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_par2() + standard);
  }
  if(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau1() >= 0){
  Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau1(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau1() + standard);
  }
  if(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau2() >= 0){
  Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau2(Transitsinglejunction1.at(itagsj1).at(itagsj2).PY_dau2() + standard);
  }

  WaitingLineforPY.push_back(Transitsinglejunction1.at(itagsj1).at(itagsj2));
  }
  standard = standard + Transitsinglejunction1.at(itagsj1).size();
}

for(int ifins = 0; ifins < finalstring.size(); ifins++){
  WaitingLineforPY.push_back(finalstring.at(ifins));
}

//std::cout <<endl<<" Let's check final entry before PY!!"<<endl;
//std::cout <<endl<< " the order of data is eNum >> id >> stat >> par1 >> par2 >> dau1 >> dau2 "<< endl;
for(int fincheck = 0; fincheck < WaitingLineforPY.size(); fincheck++){
  vector<HHparton> temp = WaitingLineforPY;
  //std::cout  << fincheck <<"   "<< temp.at(fincheck).id() <<"   "<< temp.at(fincheck).PY_stat() <<"   "<< temp.at(fincheck).PY_par1() <<"   "<< temp.at(fincheck).PY_par2() <<"   "
  //<< temp.at(fincheck).PY_dau1() <<"   "<< temp.at(fincheck).PY_dau2() << "  ( "<<  temp.at(fincheck).col()<< " , " << temp.at(fincheck).acol() << " )   " << temp.at(fincheck).x_t()<<endl;
}

for(int ifin = 0; ifin < WaitingLineforPY.size(); ifin++){
SP_prepremn.add(WaitingLineforPY[ifin]);
}

  Tempjunctions.clear(); // clear these informations for next running
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
		//const double mm_to_fm = 100000000000.0; const double fm_to_mm = 1./mm_to_fm;
		const double mm_to_fm = 1.; const double fm_to_mm = 1./mm_to_fm;
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
			//for junction mother, adding junction to list manually for color tracing appendJunction(int kind, int col0, int col1, int col2);
			if(std::abs(event[event.size()-1].id()) > 1114){
				event.appendJunction(((event[event.size()-1].id()>0) ? 1 : 2), HH_pyremn[i].PY_tag1(), HH_pyremn[i].PY_tag2(), HH_pyremn[i].PY_tag3());
			}
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
			//hadout.x(event[i].xProd()*mm_to_fm); hadout.y(event[i].yProd()*mm_to_fm); hadout.z(event[i].zProd()*mm_to_fm); hadout.x_t(event[i].tProd()*mm_to_fm);
			
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
			
			
			
			
			//***************************************************************************************************//
			int had_col = 0; int had_acol = 0; had_col = event[i].col(); had_acol = event[i].acol();
			bool info_found = false;
			
			//3 main cases: col!=acol, col==acol!=0, col==acol==0
			if(had_col!=had_acol){//col!=acol -> this hadron should trace back to a single gluon; find it and grab its spacetime info
				for(int irem=0; irem<HH_pyremn.num(); ++irem){
					if(((HH_pyremn[irem].col()==had_col) && (HH_pyremn[irem].acol()==had_acol)) || ((HH_pyremn[irem].col()==had_acol) && (HH_pyremn[irem].acol()==had_col))){//alter for glu loops
						if(HH_pyremn[irem].PY_stat()<=0){continue;}
						hadout.x(HH_pyremn[irem].x()); hadout.y(HH_pyremn[irem].y()); hadout.z(HH_pyremn[irem].z()); hadout.x_t(HH_pyremn[irem].x_t());
						info_found = true;
					}
				}
				if(!info_found){ //this is bad; could either be a hadron from *multiple* segments, or still from a single gluon under color reconnections - OR BOTH!
					//this is a first-order handling to find all the intermediate partons using the input partons -
						//it really should be done using pythia's history, but should suffice for most cases, for now...
					//will apply Dijkstra to the partons in the current string to find the shortest path from had_col to had_acol
					
					//start by finding the partons with the col/acol tags
					int ptn1 = -1; int ptn2 = -1;
					for(int irem=0; irem<HH_pyremn.num(); ++irem){
						if((HH_pyremn[irem].col() ==had_col ) && (HH_pyremn[irem].PY_stat()>0)){ptn1=irem;}
						if((HH_pyremn[irem].acol()==had_acol) && (HH_pyremn[irem].PY_stat()>0)){ptn2=irem;}
						if((ptn1>=0) && (ptn2>=0)){break;}
					}
					
					//if ptn1 or ptn2 < 0, then we're giving up on this hadron (needs pythia history to properly reconstruct the mother parton(s))
					//otherwise, we'll trace along the string from one parton end to the other (using Dijkstra) to find all mother partons
					if((ptn1>=0) && (ptn2>=0)){
						//for 0th order approx, we'll just average the positions of the known parton ends
						double pos_x, pos_y, pos_z, pos_t; pos_x=0.; pos_y=0.; pos_z=0.; pos_t=0.;
						pos_x+=HH_pyremn[ptn1].x(); pos_y+=HH_pyremn[ptn1].y(); pos_z+=HH_pyremn[ptn1].z(); pos_t+=HH_pyremn[ptn1].x_t();
						pos_x+=HH_pyremn[ptn2].x(); pos_y+=HH_pyremn[ptn2].y(); pos_z+=HH_pyremn[ptn2].z(); pos_t+=HH_pyremn[ptn2].x_t();
						pos_x/=2.; pos_y/=2.; pos_z/=2.; pos_t/=2.;
						hadout.x(pos_x); hadout.y(pos_y); hadout.z(pos_z); hadout.x_t(pos_t);
						info_found = true;
					}
					
				}
			}
			else if((had_col==had_acol) && (had_col!=0)){//col==acol -> this hadron traces back to a string segment; find the two partons and interpolate position
				//there are 2 cases here - an easy case where the color tags match the original partons, and a hard case where we need to trace history
				int ptn1 = 0; int ptn2 = 0; int ptn3 = 0; bool col_found = false; bool acol_found = false; bool ptn3_found = false;
				while(ptn1<HH_pyremn.num()){if((HH_pyremn[ptn1].col() ==had_col) && (HH_pyremn[ptn1].PY_stat()>0)){col_found =true; break;} ++ptn1;}
				while(ptn2<HH_pyremn.num()){if((HH_pyremn[ptn2].acol()==had_acol)&& (HH_pyremn[ptn2].PY_stat()>0)){acol_found=true; break;} ++ptn2;}
				
				//if we don't find color/anticolor tags in input partons, we'll need to trace along pythia event/history to find either/both
				//we can grab what pythia claims are the mothers, then search that for the color tag, which should return ONE candidate
				//repeat this until event_i <= HH_pyremn.num() - this will be our input
				//would it be better to start at the beginning, grab everything that doesn't have a terminal color tag in a hadron, and trace those to respective hadron(s)?
				
				//since these *only* should happen in junction systems at the junction, use the junction list to find the other 2 color tags, then find those in orig. ptns!
				if(!col_found){
					//find the junction with matching color, grab the other 2
					int coll[2] = {0,0};
					for(int iJ=0;iJ<pythia.event.sizeJunction();++iJ){for(int iC=0;iC<3;++iC){
						if(pythia.event.colJunction(iJ,iC)==had_col){
							if(iC==0){coll[0]=pythia.event.colJunction(iJ,1); coll[1]=pythia.event.colJunction(iJ,2);}
							else if(iC==1){coll[0]=pythia.event.colJunction(iJ,0); coll[1]=pythia.event.colJunction(iJ,2);}
							else{coll[0]=pythia.event.colJunction(iJ,0); coll[1]=pythia.event.colJunction(iJ,1);}
							col_found=true; break;
						}
					} if(col_found){break;}}
					//since pythia refuses to keep track of initial junctions given to it, we'll do it manually
					if(!col_found){
						for(int irem=0;irem<HH_pyremn.num();++irem){if(std::abs(HH_pyremn[irem].id()) > 1112){
							if(     HH_pyremn[irem].PY_tag1()==had_col){coll[0]=HH_pyremn[irem].PY_tag2(); coll[1]=HH_pyremn[irem].PY_tag3(); col_found=true; break;}
							else if(HH_pyremn[irem].PY_tag2()==had_col){coll[0]=HH_pyremn[irem].PY_tag1(); coll[1]=HH_pyremn[irem].PY_tag3(); col_found=true; break;}
							else if(HH_pyremn[irem].PY_tag3()==had_col){coll[0]=HH_pyremn[irem].PY_tag1(); coll[1]=HH_pyremn[irem].PY_tag2(); col_found=true; break;}
						}}
					}
					//now we know the anti-color tags needed, we look for them.
					if(col_found){
						bool found_tags[2] = {0,0}; ptn1 = 0;
						while(ptn1<HH_pyremn.num()){if((HH_pyremn[ptn1].acol()==coll[0]) && (HH_pyremn[ptn1].PY_stat()>0)){found_tags[0]=true; break;} ++ptn1;}
						while(ptn3<HH_pyremn.num()){if((HH_pyremn[ptn3].acol()==coll[1]) && (HH_pyremn[ptn2].PY_stat()>0)){found_tags[1]=true; break;} ++ptn3;}
						if(found_tags[0] && found_tags[1]){ptn3_found=true;}
					}
				}
				if(!acol_found){
					//check to see if there's already 3 partons - if so, then we need a warning/error!
					if(ptn3_found){JSWARN << "A hadron was found with more than 3 partonic parents for space-time info!";}
					//otherwise, works just like above...
					int coll[2] = {0,0};
					for(int iJ=0;iJ<pythia.event.sizeJunction();++iJ){for(int iC=0;iC<3;++iC){
						if(pythia.event.colJunction(iJ,iC)==had_acol){
							if(iC==0){coll[0]=pythia.event.colJunction(iJ,1); coll[1]=pythia.event.colJunction(iJ,2);}
							else if(iC==1){coll[0]=pythia.event.colJunction(iJ,0); coll[1]=pythia.event.colJunction(iJ,2);}
							else{coll[0]=pythia.event.colJunction(iJ,0); coll[1]=pythia.event.colJunction(iJ,1);}
							acol_found=true; break;
						}
					} if(acol_found){break;}}
					//since pythia refuses to keep track of initial junctions given to it, we'll do it manually
					if(!acol_found){
						for(int irem=0;irem<HH_pyremn.num();++irem){if(std::abs(HH_pyremn[irem].id()) > 1112){
							if(     HH_pyremn[irem].PY_tag1()==had_col){coll[0]=HH_pyremn[irem].PY_tag2(); coll[1]=HH_pyremn[irem].PY_tag3(); acol_found=true; break;}
							else if(HH_pyremn[irem].PY_tag2()==had_col){coll[0]=HH_pyremn[irem].PY_tag1(); coll[1]=HH_pyremn[irem].PY_tag3(); acol_found=true; break;}
							else if(HH_pyremn[irem].PY_tag3()==had_col){coll[0]=HH_pyremn[irem].PY_tag1(); coll[1]=HH_pyremn[irem].PY_tag2(); acol_found=true; break;}
						}}
					}
					//now we know the anti-color tags needed, we look for them.
					if(acol_found){
						bool found_tags[2] = {0,0}; ptn2 = 0;
						while(ptn2<HH_pyremn.num()){if((HH_pyremn[ptn2].col()==coll[0]) && (HH_pyremn[ptn1].PY_stat()>0)){found_tags[0]=true; break;} ++ptn2;}
						while(ptn3<HH_pyremn.num()){if((HH_pyremn[ptn3].col()==coll[1]) && (HH_pyremn[ptn2].PY_stat()>0)){found_tags[1]=true; break;} ++ptn3;}
						if(found_tags[0] && found_tags[1]){ptn3_found=true;}
					}
				}
				
				//turning off warning iff things work
				if(col_found && acol_found){info_found=true;}
				
				//setting position from found partons
				double pos_x, pos_y, pos_z, pos_t; pos_x=0.; pos_y=0.; pos_z=0.; pos_t=0.;
				pos_x+=HH_pyremn[ptn1].x(); pos_y+=HH_pyremn[ptn1].y(); pos_z+=HH_pyremn[ptn1].z(); pos_t+=HH_pyremn[ptn1].x_t();
				pos_x+=HH_pyremn[ptn2].x(); pos_y+=HH_pyremn[ptn2].y(); pos_z+=HH_pyremn[ptn2].z(); pos_t+=HH_pyremn[ptn2].x_t();
				if(ptn3_found){
					pos_x+=HH_pyremn[ptn3].x(); pos_y+=HH_pyremn[ptn3].y(); pos_z+=HH_pyremn[ptn3].z(); pos_t+=HH_pyremn[ptn3].x_t();
					pos_x/=3.; pos_y/=3.; pos_z/=3.; pos_t/=3.;
				}
				else{pos_x/=2.; pos_y/=2.; pos_z/=2.; pos_t/=2.;}
				hadout.x(pos_x); hadout.y(pos_y); hadout.z(pos_z); hadout.x_t(pos_t);
			}
			else if((had_col==had_acol) && (had_col==0)){//col==acol==0 -> same as previous, but pythia didn't give color tag since this was from a short string q-qbar
				double avg_x, avg_y, avg_z, avg_t; avg_x=0.; avg_y=0.; avg_z=0.; avg_t=0.;
				int n_ptns = 0; //should come out to 2 in the end, but keeping track just in case.
				for(int imot=0; imot<mothers.size(); ++imot){ if(mothers[imot] <= HH_pyremn.num()){ info_found = true; //if no mothers, then no position, and throw error!
					avg_x+=HH_pyremn[mothers[imot]-1].x(); avg_y+=HH_pyremn[mothers[imot]-1].y(); avg_z+=HH_pyremn[mothers[imot]-1].z(); avg_t+=HH_pyremn[mothers[imot]-1].x_t();
					++n_ptns;
				}}
				if(n_ptns>0){avg_x/=double(n_ptns); avg_y/=double(n_ptns); avg_z/=double(n_ptns); avg_t/=double(n_ptns);}
				hadout.x(avg_x); hadout.y(avg_y); hadout.z(avg_z); hadout.x_t(avg_t);
			}
			
			//if(!info_found){JSWARN << "  Could not find spacetime information for hadron in event!" /*<< col_found << ":" << acol_found*/;
			if(!info_found){JSWARN << "  Could not find spacetime information for hadron in event!" /*<< col_found << ":" << acol_found*/;
			hadout.x(0.); hadout.y(0.); hadout.z(0.); hadout.x_t(0.);}
			//***************************************************************************************************//
			
			//the mother procedure might skip some partons if there are junctions involved
			//this can be 'repaired' by taking a 'mother' parton, then checking over all the partons in its string! (both adding to parents / checking if thermal)
			//this is done in hadronization calling function, after this invoke_py function is finished
			HH_hadrons.add(hadout);
		}}
		
/*		event.list();
		std::cout << "\n\n";
		for(int icount=0;icount<HH_pyremn.num();++icount){
			std::cout << HH_pyremn[icount].id() << ", " << HH_pyremn[icount].col() << ":" << HH_pyremn[icount].acol() << ", " <<  
			  HH_pyremn[icount].x() << ", " << HH_pyremn[icount].y() << ", " << HH_pyremn[icount].z() << ", " << HH_pyremn[icount].x_t() << "\n";
		}
		std::cout << "\n";
		for(int icount=0;icount<HH_hadrons.num();++icount){if(HH_hadrons[icount].is_final()){continue;}
			std::cout << HH_hadrons[icount].id() << ", " << HH_hadrons[icount].x() << ", " << HH_hadrons[icount].y() << ", " << HH_hadrons[icount].z() << ", " << HH_hadrons[icount].x_t() << "\n";
		}
		std::cout << "\n";
		for(int icount=0;icount<HH_hadrons.num();++icount){if(!HH_hadrons[icount].is_final()){continue;}
			std::cout << HH_hadrons[icount].id() << ", " << HH_hadrons[icount].x() << ", " << HH_hadrons[icount].y() << ", " << HH_hadrons[icount].z() << ", " << HH_hadrons[icount].x_t() << "\n";
		}
		std::cout << "\n==================================================================================\n";
*/		
		need_hadronization = false;
	}
	
	return success;
}
