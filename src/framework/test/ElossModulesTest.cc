// quick and dirty test class implementation for Eloss modules ...

#include "ElossModulesTest.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#include "fluid_dynamics.h"

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


void Matter::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   w.lock()->Write("ElossModule Parton List: "+GetId());
   w.lock()->Write("Energy loss to be implemented accordingly ...");
}

void Matter::Exec()
{
   INFO<<"Run Matter ...";
   DEBUG<<"Qhat = "<<GetQhat();
   
   /*
   DEBUG<<"Emit Signal: jetSignal(10,20.3)";
   jetSignal(10,20.3);
   double edensity=-1;
   edensitySignal(1,edensity);
   DEBUG<< MAGENTA<<"Received edensity = "<<edensity<<" for t="<<1;
   */

   if (GetShowerInitiatingParton())
     {
       //cout<<shared_from_this().get()<<endl;
       cout<< *GetShowerInitiatingParton()<<endl;
       //PrintShowerInitiatingParton();
     }
   
   FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
   GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);
   
   DEBUG<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
   
   //check_fluid_info_ptr->Print();
   
   delete check_fluid_info_ptr;
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
  //DEBUG<<"Emit Signal: jetSignal(100,200.3)";
  //cout<<jetSignal.is_connected()<<endl;
  //jetSignal(100,200.3);

  //cout<<shared_from_this().get()<<endl;
  // Logger not working with overloaded << from Parton class ...
  // check and resolve ...
   if (GetShowerInitiatingParton())
     {
       cout<< *GetShowerInitiatingParton()<<endl;
       //PrintShowerInitiatingParton();
     }
   
   FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
   GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);  
   DEBUG<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
   
   //check_fluid_info_ptr->Print();
   
   delete check_fluid_info_ptr;
}
