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

void HadronizationModuleTest::DoHadronization(vector<shared_ptr<Parton>>& pIn, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut)
{
  INFO<<"Start Hadronizing using the PYTHIA module...";
    Pythia pythia;
    Event& event      = pythia.event;
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("PartonLevel:FSR=on");
    pythia.init();
    event.reset();
//    event.append( 21, 23, 101, 103, 0., 0., 100.0, 100 );
//    event.append( 21, 23, 103, 101, 0., 0., -100.0, 100);
    
    
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {
      double onshellE = pow(pow(pIn.at(ipart)->px(),2) + pow(pIn.at(ipart)->py(),2) + pow(pIn.at(ipart)->pz(),2) ,0.5 ) ;
      event.append(pIn.at(ipart)->pid(),23,pIn.at(ipart)->color(),pIn.at(ipart)->anti_color(),
                   pIn.at(ipart)->px(),pIn.at(ipart)->py(),pIn.at(ipart)->pz(),onshellE);
  }
    unsigned int color, anti_color;
    int pid;
    
    anti_color = pIn.at(0)->min_anti_color();
    color = pIn.at(0)->min_color();
    
    if ((color>100)&&(anti_color>100)){
        pid = 21;
    }
    else if ((color>100)&&(anti_color<100))
    {
        pid = 1;
    }
    else
    {
        pid = -1;
    }
    
   event.append(pid, 23, anti_color, color, 0.2, 0.2, 100.0, 100.0004);
    
    
      
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
    pythia.next();
   // event.list();

/*    for (unsigned int i=0; i <  event.size(); ++i)
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
    
    
  INFO<<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";

}
