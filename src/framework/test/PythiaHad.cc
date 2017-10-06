#include "PythiaHad.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

using namespace Jetscape;
using namespace Pythia8;

Pythia pythia;
Event& event      = pythia.event;
ParticleData& pdt = pythia.particleData;

PythiaHad::PythiaHad()
{
  SetId("PythiaHad");
  VERBOSE(8);
}

PythiaHad::~PythiaHad()
{
  VERBOSE(8);
}

void PythiaHad::Init()
{
  DEBUG<<"Initialize PythiaHad";
  VERBOSE(8);

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Standard settings
  pythia.readString("ProcessLevel:all = off");
  

  // ************
  // CANNOT MAKE THE XML WORK, GIVES SEG FAULT WHEN EXECUTING...
  // ************
  /*

  // XML settings
  tinyxml2::XMLElement *PythiaXmlDescription=GetHadronXML()->FirstChildElement("PythiaHad");
  tinyxml2::XMLElement *xmle;  

  string s;
  // For parsing text
  stringstream numbf(stringstream::app|stringstream::in|stringstream::out);
  numbf.setf(ios::fixed, ios::floatfield);       numbf.setf(ios::showpoint);       numbf.precision(1);
  stringstream numbi(stringstream::app|stringstream::in|stringstream::out);

  if ( !PythiaXmlDescription ) {
    WARN << "Cannot initialize Pythia Had";
    throw std::runtime_error("Cannot initialize Pythia Had");
  }

  xmle = PythiaXmlDescription->FirstChildElement( "name" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  s = xmle->GetText();
  SetId(s);
  cout << s << endl;

  //To use in future when choosing lund string parameters
  //xmle = PythiaXmlDescription->FirstChildElement( "pTHatMin" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  //xmle->QueryDoubleText(&pTHatMin);
  //xmle = PythiaXmlDescription->FirstChildElement( "pTHatMax" ); if ( !xmle ) throw std::runtime_error("Cannot parse xml");
  //xmle->QueryDoubleText(&pTHatMax);
  
  //VERBOSE(7) <<"Pythia Gun with "<< pTHatMin << " < pTHat < " << pTHatMax ;
  
  //numbf.str("PhaseSpace:pTHatMin = "); numbf << pTHatMin;
  //readString ( numbf.str() );
  //numbf.str("PhaseSpace:pTHatMax = "); numbf << pTHatMax;
  //readString ( numbf.str() );

  // random seed
  // xml limits us to unsigned int :-/ -- but so does 32 bits Mersenne Twist
  tinyxml2::XMLElement *RandomXmlDescription=JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Random" );
  pythia.readString("Random:setSeed = on");
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
  pythia.readString( numbi.str() );

  xmle = PythiaXmlDescription->FirstChildElement( "LinesToRead" );
  if ( xmle ){
    std::stringstream lines; lines << xmle->GetText();
    int i=0;
    while(std::getline(lines,s,'\n')){
      if( s.find_first_not_of (" \t\v\f\r") == s.npos ) continue; // skip empty lines
      VERBOSE(7) <<  "Also reading in: " << s;
      pythia.readString (s);
    }
  }
  
  */

  // And initialize
  pythia.init();

}

void PythiaHad::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   w.lock()->WriteComment("Hadronization Module : "+GetId());
   w.lock()->WriteComment("Hadronization to be implemented accordingly ...");
}

void PythiaHad::DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  INFO<<"Start Hadronizing using Lund string model...";
  event.reset();
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {  
    int ide=pIn[ipart]->pid();
    double px=pIn[ipart]->px();
    double py=pIn[ipart]->py();
    double pz=pIn[ipart]->pz();
    double ee=pIn[ipart]->e();
    double mm=pdt.m0(int(ide));
    ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
    //Need to provide color singlet object
    //Dummy setup that provides an antiquark for every quark
    int col, acol;
    if (ide>0) col=102, acol=0;
    else col=0; acol=102;
    if (ide==22) col=0, acol=0;
    if (ide==21) ide=22, col=0, acol=0; 
    event.append(int(ide),23,int(col),int(acol),px,py,pz,ee,mm);
    event.append(-int(ide),23,int(acol),int(col),-px,-py,-pz,ee,mm); 
  }
  pythia.next();
  for (unsigned int ipart=0; ipart < event.size(); ++ipart)
  {
    if (event[ipart].isFinal())
    {
      int ide=pythia.event[ipart].id();
      FourVector p(pythia.event[ipart].px(),pythia.event[ipart].py(),pythia.event[ipart].pz(),pythia.event[ipart].e());
      FourVector x;
      hOut.push_back(std::make_shared<Hadron> (Hadron (0,ide,0,p,x)));
      INFO << "Produced Hadron has id = " << pythia.event[ipart].id();
    }
  } 
  INFO<<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";

}
