#include "HadronizationModuleTest.h"
#include "JetScapeLogger.h"


using namespace Jetscape;
using namespace Pythia8;


HadronizationModuleTest::HadronizationModuleTest()
{
  SetId("MyHadroTest");
  VERBOSE(8);
}

HadronizationModuleTest::~HadronizationModuleTest()
{
  VERBOSE(8);
}

void HadronizationModuleTest::Init()
{

}

void HadronizationModuleTest::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   w.lock()->WriteComment("Hadronization Module : "+GetId());
   w.lock()->WriteComment("Hadronization to be implemented accordingly ...");
}

void HadronizationModuleTest::DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  VERBOSE(2)<<"Start Hadronizing using the PYTHIA module...";
  // cmldir will be detected automatically and override this setting anyway.
  // Only set because it's needed to suppress the banner
  string xmlDir = "DONTUSETHIS";
  bool printBanner = false;
  if ( JetScapeLogger::Instance()->GetVerboseLevel()>7 ) printBanner = true;
  Pythia pythia (xmlDir,printBanner);
    
  Event& event      = pythia.event;

  // Show initialization at DEBUG or high verbose level
  pythia.readString("Init:showProcesses = off");
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showMultipartonInteractions = off");
  pythia.readString("Init:showChangedParticleData = off");
  if ( JetScapeLogger::Instance()->GetDebug() || JetScapeLogger::Instance()->GetVerboseLevel()>7 ) {
    pythia.readString("Init:showProcesses = on");
    pythia.readString("Init:showChangedSettings = on");
    pythia.readString("Init:showMultipartonInteractions = on");
    pythia.readString("Init:showChangedParticleData = on");
  }

  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0"); 
  pythia.readString("Next:numberShowProcess = 0"); 
  pythia.readString("Next:numberShowEvent = 0"); 

  pythia.readString("ProcessLevel:all = off");
  pythia.readString("PartonLevel:FSR=on");
  pythia.init();
  event.reset();
  double pz = 100.0;
//    event.append( 21, 23, 101, 103, 0., 0., 100.0, 100 );
//    event.append( 21, 23, 103, 101, 0., 0., -100.0, 100);
    
  DEBUG << "&&&&&&&&&&&&&&&&&&& the number of showers are: " << shower.size();
  for(unsigned int ishower=0; ishower <  shower.size(); ++ishower)  
{
    


  DEBUG << "&&&&&&&&&&&&&&&&&&& there are " << shower.at(ishower).size() << " partons in the shower number " << ishower;
  for(unsigned int ipart=0; ipart <  shower.at(ishower).size(); ++ipart)
  {
      double onshellE = pow(pow(shower.at(ishower).at(ipart)->px(),2) + pow(shower.at(ishower).at(ipart)->py(),2) + pow(shower.at(ishower).at(ipart)->pz(),2) ,0.5 ) ;
      event.append(shower.at(ishower).at(ipart)->pid(),23,shower.at(ishower).at(ipart)->color(),shower.at(ishower).at(ipart)->anti_color(),
                   shower.at(ishower).at(ipart)->px(),shower.at(ishower).at(ipart)->py(),shower.at(ishower).at(ipart)->pz(),onshellE);
  }
    unsigned int color, anti_color;
    int pid;
    
    anti_color = shower.at(ishower).at(0)->min_anti_color();
    color = shower.at(ishower).at(0)->min_color();
    
    if ((color>100)&&(anti_color>100)){
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
    event.append(pid, 23, anti_color, color, 0.2, 0.2, pz, 100.0004);
    
    
      
/*    if(ipart%2==0)
    {
      hOut.push_back(std::dynamic_pointer_cast<Hadron> (pIn.at(ipart)));
    }
    else
    {
      pOut.push_back(pIn.at(ipart));
    }
  }
 */

   // event.list();

/*   for (unsigned int i=0; i <  event.size(); ++i)
    {
        if (event[i].status()==91)
        {
            double x[4];
            x[0]=x[1]=x[2]=x[3]=0.0;
            
            
            
//            Hadron hadron(i,event[i].id(),91,event[i].pT(),event[i].eta(), event[i].phi(), event[i].e(), x);
            
            hOut.push_back(make_shared<Hadron>(i,event[i].id(),91,event[i].pT(),event[i].eta(), event[i].phi(), event[i].e(), x));
        }
        
    }
*/
    
    
   VERBOSE(2) <<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";
}
        pythia.next();
    
    shower.clear();
    
}
