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

#include "ColoredHadronization.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "tinyxml2.h"


using namespace Jetscape;
using namespace Pythia8;

Pythia8::Pythia ColoredHadronization::pythia ("IntentionallyEmpty",false);


ColoredHadronization::ColoredHadronization()
{
  SetId("MyHadroTest");
  VERBOSE(8);
}

ColoredHadronization::~ColoredHadronization()
{
  VERBOSE(8);
}

void ColoredHadronization::Init()
{
  tinyxml2::XMLElement *hadronization= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("JetHadronization" );
    
  if ( !hadronization ) {
    JSWARN << "Couldn't find tag Jet Hadronization";
    throw std::runtime_error ("Couldn't find tag Jet Hadronization");
  }
  if (hadronization) {
    string s = hadronization->FirstChildElement( "name" )->GetText();
    JSDEBUG << s << " to be initializied ...";
        
    double p_read_xml = 10000 ;
        
    hadronization->FirstChildElement("eCMforHadronization")->QueryDoubleText(&p_read_xml);
    p_fake = p_read_xml;
        
    VERBOSE(2)<<"Start Hadronizing using the PYTHIA module...";  
        
    // Show initialization at DEBUG or high verbose level
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("Init:showChangedParticleData = off");
    if ( JetScapeLogger::Instance()->GetDebug() || JetScapeLogger::Instance()->GetVerboseLevel()>2 ) {
      pythia.readString("Init:showProcesses = on");
      pythia.readString("Init:showChangedSettings = on");
      pythia.readString("Init:showMultipartonInteractions = on");
      pythia.readString("Init:showChangedParticleData = on");
    }
        
    // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    if ( JetScapeLogger::Instance()->GetDebug() || JetScapeLogger::Instance()->GetVerboseLevel()>2 ) {
      pythia.readString("Next:numberShowInfo = 1");
      pythia.readString("Next:numberShowProcess = 1");
      pythia.readString("Next:numberShowEvent = 1");
    }
        
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("PartonLevel:FSR=off");
    pythia.readString("HadronLevel:Decay = off");
    pythia.init();
        
  }
    
}

void ColoredHadronization::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  auto f = w.lock();
  if ( !f ) return;
  f->WriteComment("Hadronization Module : "+GetId());
  f->WriteComment("Hadronization to be implemented accordingly ...");
}

void ColoredHadronization::DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{

  Event& event = pythia.event;
  event.reset();
  double pz = p_fake;
   
  JSDEBUG << "&&&&&&&&&&&&&&&&&&& the number of showers are: " << shower.size();
  for(unsigned int ishower=0; ishower <  shower.size(); ++ishower)
  {
      JSDEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size() << " partons in the shower number " << ishower;
      for(unsigned int ipart=0; ipart <  shower.at(ishower).size(); ++ipart)
      {
          double onshellE = pow(pow(shower.at(ishower).at(ipart)->px(),2) + pow(shower.at(ishower).at(ipart)->py(),2) + pow(shower.at(ishower).at(ipart)->pz(),2) ,0.5 ) ;
          
            if ( shower.at(ishower).at(ipart)->pid()==22 )
            {
                
                double blurb;
                JSINFO << BOLDYELLOW << " photon found with " ;
                JSINFO << BOLDYELLOW << "px = " << shower.at(ishower).at(ipart)->px();
                //cin >> blurb;
            }
          event.append(shower.at(ishower).at(ipart)->pid(),23,shower.at(ishower).at(ipart)->color(),shower.at(ishower).at(ipart)->anti_color(),
		       shower.at(ishower).at(ipart)->px(),shower.at(ishower).at(ipart)->py(),shower.at(ishower).at(ipart)->pz(),onshellE);
      }
      unsigned int color, anti_color;
      int pid;
  
      anti_color = shower.at(ishower).at(0)->min_anti_color();
      color = shower.at(ishower).at(0)->min_color();
  
      if ((color>100)&&(anti_color>100))
      {
          pid = 21;
      }
      else if ((color>100)&&(anti_color<100))
      {
          pid = -1;
      }
      else
      {
          pid = 1;
      }
  
      pz = -1*pz;
      event.append(pid, 23, anti_color, color, 0.2, 0.2, pz, sqrt(pz*pz + 0.08));       
    
      VERBOSE(2) <<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";
    }
  
  pythia.next();
  // event.list();
  
  unsigned int ip=hOut.size();
  for (unsigned int i=0; i<event.size(); ++i){
    if ( !event[i].isFinal() )   continue;
//    if ( !event[i].isHadron() )  continue;
    if(fabs(event[i].eta())>20)  continue; //To prevent "nan" from propagating, very rare though
    
    double x[4] = {0,0,0,0};
    hOut.push_back(make_shared<Hadron>(ip,event[i].id(),event[i].status(),event[i].pT(),event[i].eta(), event[i].phi(), event[i].e(), x));
    ++ip;    
  }
  
  shower.clear();
    
}
