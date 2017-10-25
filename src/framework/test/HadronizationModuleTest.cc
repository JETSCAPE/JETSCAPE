#include "HadronizationModuleTest.h"
#include "JetScapeLogger.h"


using namespace Jetscape;


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
  INFO<<"Start Hadronizing using my test module...";
  for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
  {
    if(ipart%2==0)
    {
      hOut.push_back(std::dynamic_pointer_cast<Hadron> (pIn.at(ipart)));
    }
    else
    {
      pOut.push_back(pIn.at(ipart));
    }
  }
  INFO<<"There are " << hOut.size() << " Hadrons and " << pOut.size() << " partons after Hadronization";

}
