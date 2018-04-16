/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
/** Create a pythia collision at a specified point and return the two inital hard partons
    @group JetScape (modular/task) based framework
    @version Revision 0.1
    @date Jun 29, 2017
*/

#include "PythiaGun.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;


PythiaGun::~PythiaGun()
{
  VERBOSE(8);
}

void PythiaGun::InitTask()
{

  JSDEBUG<<"Initialize PythiaGun"; 
  VERBOSE(8);

  // Show initialization at INFO level
  readString("Init:showProcesses = off");
  readString("Init:showChangedSettings = off");
  readString("Init:showMultipartonInteractions = off");
  readString("Init:showChangedParticleData = off");
  if ( JetScapeLogger::Instance()->GetInfo() ) {
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
  readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration  
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  readString("PartonLevel:ISR = on");
  readString("PartonLevel:MPI = on");
  readString("PartonLevel:FSR = off");
  readString("PromptPhoton:all=off");
  readString("WeakSingleBoson:all=off");
  readString("WeakDoubleBoson:all=off");

  tinyxml2::XMLElement *PythiaXmlDescription=GetHardXML()->FirstChildElement("PythiaGun");
  tinyxml2::XMLElement *xmle; 

  string s;
  // For parsing text
  stringstream numbf(stringstream::app|stringstream::in|stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);       numbf.setf(ios::showpoint);       numbf.precision(1);
  stringstream numbi(stringstream::app|stringstream::in|stringstream::out);
    
  if ( !PythiaXmlDescription ) {
    WARN << "Cannot initialize Pythia Gun";
    throw std::runtime_error("Cannot initialize Pythia Gun");
  }

  xmle = PythiaXmlDescription->FirstChildElement( "name" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  s = xmle->GetText();
  SetId(s);
  // cout << s << endl;
  
  xmle = PythiaXmlDescription->FirstChildElement( "pTHatMin" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&pTHatMin);
  xmle = PythiaXmlDescription->FirstChildElement( "pTHatMax" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&pTHatMax);
  
  VERBOSE(7) <<"Pythia Gun with "<< pTHatMin << " < pTHat < " << pTHatMax ;
  
  numbf.str("PhaseSpace:pTHatMin = "); numbf << pTHatMin;
  readString ( numbf.str() );
  numbf.str("PhaseSpace:pTHatMax = "); numbf << pTHatMax;
  readString ( numbf.str() );
  
  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
  readString("Random:setSeed = on");
  numbi.str("Random:seed = ");
  unsigned int seed = 0;
  if ( RandomXmlDescription ){
    xmle = RandomXmlDescription->FirstChildElement( "seed" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
    xmle->QueryUnsignedText(&seed);
  } else {
    WARN << "No <Random> element found in xml, seeding to 0";
  }
  VERBOSE(7) <<"Seeding pythia to "<< seed ;
  numbi << seed;
  readString( numbi.str() );
  
  // Species
  readString("Beams:idA = 2212");
  readString("Beams:idB = 2212");
  
  // Energy
  xmle = PythiaXmlDescription->FirstChildElement( "eCM" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  xmle->QueryDoubleText(&eCM);   
  numbf.str("Beams:eCM = "); numbf << eCM;
  readString ( numbf.str() );
  
  xmle = PythiaXmlDescription->FirstChildElement( "LinesToRead" );
  if ( xmle ){
    std::stringstream lines;
    lines << xmle->GetText();
    int i=0;
    while(std::getline(lines,s,'\n')){
      if( s.find_first_not_of (" \t\v\f\r") == s.npos ) continue; // skip empty lines
      VERBOSE(7) <<  "Also reading in: " << s;
      readString (s);
    }
  }
  
  // And initialize
  init(); // Pythia>8.1
  
}

void PythiaGun::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();
  //Reading vir_factor from xml for MATTER
   tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss ) {
    WARN << "Couldn't find tag Eloss";
    throw std::runtime_error ("Couldn't find tag Eloss");    
  }
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
  if ( !matter ) {
    WARN << "Couldn't find tag Eloss -> Matter";
    throw std::runtime_error ("Couldn't find tag Eloss -> Matter");
  }
  double vir_factor;
  matter->FirstChildElement("vir_factor")->QueryDoubleText(&vir_factor);

  bool flag62=false;
  vector<Pythia8::Particle> p62;
  

  // sort by pt
  struct greater_than_pt {
    inline bool operator() (const Pythia8::Particle& p1, const Pythia8::Particle& p2) {
      return ( p1.pT() > p2.pT());
    }
  };

  do{
    next();
    p62.clear();
    
    // pTarr[0]=0.0; pTarr[1]=0.0;
    // pindexarr[0]=0; pindexarr[1]=0;

    for(int parid=0; parid<event.size(); parid++){
      if ( parid<3 )continue;      // 0, 1, 2: total event and beams      
      Pythia8::Particle& particle = event[parid];

      // only accept particles after MPI
      if ( particle.status()!=62 ) continue;
       
      // only accept gluons and quarks
      if ( fabs( particle.id() ) > 3 && particle.id() !=21 ) continue;

      // reject rare cases of very soft particles that don't have enough e to get
      // reasonable virtuality
      if ( particle.pT() < 1.0/sqrt(vir_factor) ) continue;
  
      // accept
      p62.push_back( particle );

    }

    // if you want at least 2
    if ( p62.size() < 2 ) continue;
    //if ( p62.size() < 1 ) continue;
    
    // Now have all candidates, sort them
    // sort by pt
    std::sort( p62.begin(), p62.end(), greater_than_pt() );
    // // check...
    // for (auto& p : p62 ) cout << p.pT() << endl;
    
    flag62=true;

  }while(!flag62);



  //next();
  //JSDEBUG<<"Number of Pythia partons = "<<event.size();


  // stringstream numbi(stringstream::app|stringstream::in|stringstream::out);
  // numbi.str("delme");
  // numbi << GetCurrentEvent();
  // numbi << ".csv";
  
  // ofstream delme (numbi.str());
  // for (int i=0;i<100000;++i){
  //   auto idx = dist(engine);
  //   auto coord = ini->CoordFromIdx( idx );
  //   // cout << " =================================> " <<  idx
  //   // 	 << " x = " << get<0>(coord) 
  //   // 	 << " y = " << get<1>(coord) 
  //   // 	 << " z = " << get<2>(coord)
  //   // 	 << std::endl;

  //   delme << get<0>(coord) <<","<<get<1>(coord) << endl;
  // }
  // delme.close();

  double p[4], xLoc[4];

  // This location should come from an initial state
  for (int i=0;i<=3; i++) {
    xLoc[i] = 0.0;
  };

  // // Roll for a starting point
  // // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  // std::random_device device;
  // std::mt19937 engine(device()); // Seed the random number engine

  
  if (!ini) {
      WARN << "No initial state module, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    auto num_bin_coll = ini->GetNumOfBinaryCollisions();
    if ( num_bin_coll.size()==0 ){
      WARN << "num_of_binary_collisions is empty, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
    } else {	 
      std::discrete_distribution<> dist( begin(num_bin_coll),end(num_bin_coll) ); // Create the distribution
    
      // Now generate values
      auto idx = dist( *GetMt19937Generator() );
      auto coord = ini->CoordFromIdx( idx );
      xLoc[1] = get<0>( coord );
      xLoc[2] = get<1>( coord );
    }
  }
    

  // Loop through particles

  // Only top two
  //for(int np = 0; np<2; ++np){

  // Accept them all
  for(int np = 0; np<p62.size(); ++np){
    Pythia8::Particle& particle = p62.at( np );

    VERBOSE(7)<<"Adding particle with pid = " << particle.id()
	      <<" at x=" << xLoc[1]
	      <<", y=" << xLoc[2]
	      <<", z=" << xLoc[3];
    
    VERBOSE(7) <<"Adding particle with pid = " << particle.id()
	       << ", pT = " << particle.pT()
	       << ", y = " << particle.y()
	       << ", phi = " << particle.phi()
	       << ", e = " << particle.e();
    VERBOSE(7) <<" at x=" << xLoc[1]
	       <<", y=" << xLoc[2]
	       <<", z=" << xLoc[3];
    
    AddParton(make_shared<Parton>(0, particle.id(),0,particle.pT(),particle.y(),particle.phi(),particle.e(),xLoc) );

  }
  

 
  VERBOSE(8)<<GetNHardPartons();
}
