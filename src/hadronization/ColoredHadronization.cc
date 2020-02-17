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

// Register the module with the base class
RegisterJetScapeModule<ColoredHadronization> ColoredHadronization::reg("ColoredHadronization");

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
  
  std::string s = GetXMLElementText({"JetHadronization", "name"});
  JSDEBUG << s << " to be initializied ...";
  
  double p_read_xml = GetXMLElementDouble({"JetHadronization", "eCMforHadronization"});
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
  pythia.readString("HadronLevel:Decay = on");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10.0");
  pythia.init();
    
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
                
                JSINFO << BOLDYELLOW << " photon found in colored hadronization with " ;
                JSINFO << BOLDYELLOW << "px = " << shower.at(ishower).at(ipart)->px();
                //cin >> blurb;
            }
          event.append(shower.at(ishower).at(ipart)->pid(),23,shower.at(ishower).at(ipart)->color(),shower.at(ishower).at(ipart)->anti_color(),
		       shower.at(ishower).at(ipart)->px(),shower.at(ishower).at(ipart)->py(),shower.at(ishower).at(ipart)->pz(),onshellE);
      }

	  //first, find unpaired color and anticolor tags.
	  std::vector<int> cols; std::vector<int> acols;
	  for(unsigned int ipart=0; ipart <  shower.at(ishower).size(); ++ipart){
		  if(shower.at(ishower).at(ipart)->pid()        == 22){continue;}
		  if(shower.at(ishower).at(ipart)->color()      != 0 ){ cols.push_back(shower.at(ishower).at(ipart)->color());}
		  if(shower.at(ishower).at(ipart)->anti_color() != 0 ){acols.push_back(shower.at(ishower).at(ipart)->anti_color());}
	  }
	  //the outcomes are: 1-unpaired color tag, 2-unpaired anticolor tag, 3-both an unpaired color & anticolor tag, 4-no unpaired tags
	  //1-add an antiquark, 2-add a quark, 3-add a gluon, 4-add nothing (possibly photon only event)
	  int icol=0;
	  while(icol<cols.size()){
		  bool foundpair = false;
		  for(int iacol=0;iacol<acols.size();++iacol){if(cols[icol]==acols[iacol]){cols.erase(cols.begin()+icol); acols.erase(acols.begin()+iacol); foundpair=true; continue;}}
		  if(!foundpair){++icol;}
	  }
	  
	  int pid=0; int color=0; int anti_color=0;
	  if(     (cols.size()>0 ) && (acols.size()>0 )){pid=21; color=cols[0]; anti_color=acols[0];}
	  else if((cols.size()>0 ) && (acols.size()==0)){pid=-1; color=cols[0]; anti_color=0       ;}
	  else if((cols.size()==0) && (acols.size()>0 )){pid= 1; color=0      ; anti_color=acols[0];}
	  
	  if(pid != 0){pz = -1*pz; event.append(pid, 23, anti_color, color, 0.2, 0.2, pz, sqrt(pz*pz + 0.08));}    
    
    VERBOSE(2) <<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";
  }
	
  //there still may be color tag duplicates - will SegFault if color_reconnections is ever invoked.
  //this should be fixed *here*, before pythia.next() below, if that's ever a concern.
  //scan over list of color & anticolor tags - if we find a duplicate, set it to a new value >0, <101 (will fail for more than 100 duplicates)
  std::vector<std::vector<int>> col_instances;
  for(unsigned int i=0; i<event.size(); ++i){
	  bool newcol = true; bool newacol = true;
	  for(int icols=0;icols<col_instances.size();++icols){
		  if(event[i].col()  == col_instances[icols][0]){++col_instances[icols][1]; newcol =false;}
		  if(event[i].acol() == col_instances[icols][0]){++col_instances[icols][1]; newacol=false;}
	  }
	  if(newcol &&  (event[i].col()  != 0)){std::vector<int> tmpcol; tmpcol.push_back(event[i].col());  tmpcol.push_back(1); col_instances.push_back(tmpcol);}
	  if(newacol && (event[i].acol() != 0)){std::vector<int> tmpcol; tmpcol.push_back(event[i].acol()); tmpcol.push_back(1); col_instances.push_back(tmpcol);}
  }
  col_instances.erase(std::remove_if(col_instances.begin(), col_instances.end(),[](const std::vector<int>& val){return val[1] <= 2;}), col_instances.end());
  int updcol = 1;
  while(col_instances.size() > 0){
	  int nupd = 2;
	  for(unsigned int i=event.size()-1; i>=0; --i){
		  if(col_instances[0][0]==event[i].col() ){event[i].col( updcol); --nupd;}
		  if(col_instances[0][0]==event[i].acol()){event[i].acol(updcol); --nupd;}
		  if(nupd == 0){break;}
	  }
	  ++updcol; col_instances[0][1] -= 2;
	  col_instances.erase(std::remove_if(col_instances.begin(), col_instances.end(),[](const std::vector<int>& val){return val[1] <= 2;}), col_instances.end());
  }
  
  pythia.next();
  // event.list();
  
  unsigned int ip=hOut.size();
  for (unsigned int i=0; i<event.size(); ++i){
    if ( !event[i].isFinal() )   continue;
    //if ( !event[i].isHadron() )  continue;
    if(fabs(event[i].eta())>20)  continue; //To prevent "nan" from propagating, very rare though
    
    double x[4] = {0,0,0,0};
    hOut.push_back(make_shared<Hadron>(ip,event[i].id(),event[i].status(),event[i].pT(),event[i].eta(), event[i].phi(), event[i].e(), x));
    ++ip;    
  }
  
  shower.clear();
    
}
