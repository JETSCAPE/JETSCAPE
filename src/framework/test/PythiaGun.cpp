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
    
  if ( PythiaXmlDescription ) {
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

    //random seed
    readString("Random:setSeed = on");
    numbi.str("Random:seed = ");
    numbi << 0;
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
    
    // // DEBUG
    // next();
    // cout << event.size() << endl;

  } else {
    WARN << "Cannot initialize Pythia Gun";
    throw std::runtime_error("Cannot initialize Pythia Gun");
  }
}
 
void PythiaGun::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

  next();
  DEBUG<<"Number of Pythia partons = "<<event.size();

  double p[4], xLoc[4];

  /// @TODO: This location should come from an initial state
  for (int i=0;i<=3; i++) {
    xLoc[i] = 0.0;
  };

  // Loop through particles
  for ( int nP = 0; nP<event.size() ; ++nP ){
    if ( nP<3 )continue;      // 0, 1, 2: total event and beams      
    Pythia8::Particle& particle = event[nP];
    if ( particle.status()==-23 ){
      AddParton(make_shared<Parton>(0, particle.id(),0,particle.pT(),particle.y(),particle.phi(),particle.e(),xLoc) );
    }
  }

 
  VERBOSE(8)<<GetNHardPartons();
}
