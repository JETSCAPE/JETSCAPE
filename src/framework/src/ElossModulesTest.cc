// quick and dirty test class implementation for Eloss modules ...

#include "ElossModulesTest.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#define MAGENTA "\033[35m"

Matter::Matter()
{
  SetId("Matter");
  VERBOSE(8);
}

Matter::~Matter()
{
  VERBOSE(8);
}

void Matter::Init()
{
  INFO<<"Intialize Matter ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );  
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
 
  if (matter)
    {
      //cout<<matter->FirstChildElement( "name" )->GetText() <<endl;
      //DEBUG << __PRETTY_FUNCTION__ <<":"<<__LINE__ << " ";//<<matter->FirstChildElement( "name" )->GetText();
      
      string s = matter->FirstChildElement( "name" )->GetText();
    
      DEBUG << s << " to be initilizied ...";

      double m_qhat=-99.99;
      matter->FirstChildElement("qhat")->QueryDoubleText(&m_qhat);
      SetQhat(m_qhat);
      
      DEBUG  << s << " with qhat = "<<GetQhat();      	
    }
  else
    {
      WARN << " : Matter not properly initialized in XML file ...";
      exit(-1);
    }
}

void Matter::Exec()
{
   INFO<<"Run Matter ...";
   DEBUG<<"Qhat = "<<GetQhat();
   DEBUG<<"Emit Signal: jetSignal(10,20.3)";
   jetSignal(10,20.3);
   double edensity=-1;
   edensitySignal(1,edensity);
   DEBUG<< MAGENTA<<"Received edensity = "<<edensity<<" for t="<<1;
}

// ----------------------

Martini::Martini()
{
  SetId("Martini");
  VERBOSE(8);
}

Martini::~Martini()
{
  VERBOSE(8);
}

void Martini::Init()
{
  INFO<<"Intialize Martini ...";
}

void Martini::Exec()
{
  INFO<<"Run Martini ...";
  DEBUG<<"Qhat = "<<GetQhat();
  DEBUG<<"Emit Signal: jetSignal(100,200.3)";
  //cout<<jetSignal.is_connected()<<endl;
  jetSignal(100,200.3);
}
