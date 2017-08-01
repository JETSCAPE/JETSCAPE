/** Create a pythia collision at a specified point and return the two inital hard partons

    Based on PGun and JSPythia

    @group JetScape (modular/task) based framework
    @author Kolja Kauder
    @version Revision 0.1
    @date Jun 29, 2017
*/


#include "PythiaGun.hpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>

PythiaGun::~PythiaGun()
{
  VERBOSE(8);
}

void PythiaGun::InitTask()
{

  DEBUG<<"Initialize PythiaGun"; 
  VERBOSE(8);
  
  // No event record printout.
  readString("Next:numberShowInfo = 0"); 
  readString("Next:numberShowProcess = 0"); 
  readString("Next:numberShowEvent = 0"); 

  // Standard settings
  readString("HardQCD:all = on"); // will repeat this line in the xml for demonstration  
  readString("HadronLevel:Decay = off");
  readString("HadronLevel:all = off");
  readString("PartonLevel:ISR = on");
  readString("PartonLevel:FSR = off");
  
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
  cout << s << endl;
  
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
  // xml limits us to unsigned int :-/
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
    std::stringstream lines; lines << xmle->GetText();
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

  next();
  DEBUG<<"Number of Pythia partons = "<<event.size();


  // stringstream numbi(stringstream::app|stringstream::in|stringstream::out);
  // numbi.str("delme");
  // numbi << GetCurrentEvent();
  // numbi << ".csv";
  
  // ofstream delme (numbi.str());
  // for (int i=0;i<100000;++i){
  //   auto idx = dist(engine);
  //   auto coord = ini->coord_from_idx( idx );
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

  // Roll for a starting point
  // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
  std::random_device device;
  std::mt19937 engine(device()); // Seed the random number engine

  if (!ini) {
      WARN << "No initial state module, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    auto num_bin_coll = ini->get_num_of_binary_collisions();
    if ( num_bin_coll.size()==0 ){
      WARN << "num_of_binary_collisions is empty, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
    } else {	 
      std::discrete_distribution<> dist( begin(num_bin_coll),end(num_bin_coll) ); // Create the distribution
    
      // Now generate values
      auto idx = dist(engine);
      auto coord = ini->coord_from_idx( idx );
      xLoc[1] = get<0>( coord );
      xLoc[2] = get<1>( coord );
    }
  }
    

  // Loop through particles
  for ( int nP = 0; nP<event.size() ; ++nP ){
    if ( nP<3 )continue;      // 0, 1, 2: total event and beams      
    Pythia8::Particle& particle = event[nP];
    if ( particle.status()==-23 ){
      // cout << "particle.id()=" << particle.id() << endl;
      VERBOSE(7)<<"Adding particle with pid = " << particle.id()
		<<" at x=" << xLoc[1]
		<<", y=" << xLoc[2]
		<<", z=" << xLoc[3];

      AddParton(make_shared<Parton>(0, particle.id(),0,particle.pT(),particle.y(),particle.phi(),particle.e(),xLoc) );
    }
  }

 
  VERBOSE(8)<<GetNHardPartons();
}
