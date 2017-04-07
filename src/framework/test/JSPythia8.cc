//Pythia8 Test ...
#include "JSPythia8.h"

JSPythia8::JSPythia8() : HardProcess() , Pythia()
{
  SetId("JSPythia8");
  VERBOSE(8);
}

JSPythia8::~JSPythia8()
{
  VERBOSE(8);
}

void JSPythia8::InitTask()
{
  
  DEBUG<<"Initialize Pythia8 Brick (Test) ...";
  VERBOSE(8);

  /*
  tinyxml2::XMLElement *py8=GetHardXML()->FirstChildElement("Pythia8");

  if (py8)
    {
      string s = pgun->FirstChildElement( "name" )->GetText();
      DEBUG << s << " to be initilizied ...";
      
      pgun->FirstChildElement("pT")->QueryDoubleText(&fixed_pT);

      DEBUG << s << " with fixed pT = "<<fixed_pT;
      INFO<<"Parton Gun with fixed pT = "<<fixed_pT;
      
    }
  else
    {
      WARN << " : Pythia8 not properly initialized in XML file ...";
      exit(-1);
    }
  */

  // Read in from our xml file ...
  // to be implemented (see above)...
  
  // Generator. Shorthand for event.
  //char* ptHatMin="100.0";
  //char* ptHatMax="200.0";

  // Process selection.
  readString("HardQCD:all = on");
  readString("HadronLevel:Decay = off");
  readString("PartonLevel:ISR = on");
  readString("PartonLevel:FSR = off");

  string pmin="PhaseSpace:pTHatMin = ";
  string pmax="PhaseSpace:pTHatMax = ";
  pmin.append("100");
  pmax.append("200");

  readString(pmin);
  readString(pmax);

  // No event record printout.
  readString("Next:numberShowInfo = 0");
  readString("Next:numberShowProcess = 0");
  readString("Next:numberShowEvent = 0");

  //pythia.readString("Tune:pp = 5");
  //pythia.readString("SigmaProcess:Kfactor = 0.7");

  //random seed
  //readString("Random:setSeed = on");
  //readString("Random:seed = 0");

  //settings.listAll();
  
  //initialization ...
  //init(2212, 2212,5020.);  // Pythia < 8.2
  init(); // Pythia>8.1
}
 
void JSPythia8::Exec()
{
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

  next();

  DEBUG<<"Number of Pythia partons = "<<event.size();

  //DUMMY, to be added accordingly, later ...
  
  for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
      // Add the partons needed (5 and 6 or others !???)
      
    }
  
  
  // DEBUG: Test ...
  
  double pAssign[4], xLoc[4];
  for (int i=0;i<=3; i++) {
     xLoc[i] = 0.0;
   };

   pAssign[0] = 11.0; 
   pAssign[3] = 10.0;
   pAssign[1] = pAssign[2] = 1.0;

   for (int i=0;i<2;i++)
     {
       AddParton(make_shared<Parton>(1,21,0,pAssign,xLoc));
     }
  
  
  VERBOSE(8)<<GetNHardPartons();
  
  if (GetNHardPartons()<1)
    {WARN<<"Pythia 8 initial hard partons has to be implemented accordingly ..."; exit(-1);}
}
