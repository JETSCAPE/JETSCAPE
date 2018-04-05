
#include "Hadronization.h"
#include "JetScapeLogger.h"
#include <string>
#include <vector>
#include <iostream>
#include "JetScapeSignalManager.h"
#include "JetScapeWriterStream.h"

using namespace std;

namespace Jetscape {

Hadronization::Hadronization()
{
  TransformPartonsConnected = false;
}

Hadronization::~Hadronization()
{
}

void Hadronization::Clear()
{
  VERBOSESHOWER(8);
  outHadrons.clear();  
}

void Hadronization::Init()
{
  JetScapeModuleBase::Init();
  
  // May need to read some configuration info from XML file here
  

  if (GetNumberOfTasks()<1)
  {
    WARN << " : No valid Hadronization modules found ...";
    exit(-1);
  }
 
  INFO<<"Found "<<GetNumberOfTasks()<<" Hadronization Tasks/Modules Initialize them ... ";
  JetScapeTask::InitTasks();

}

void Hadronization::DoHadronize()
{
  VERBOSE(2)<<"Get Recombination Partons...";

  if(inPartons.size()>0)
  {
    VERBOSE(2)<<"There are Partons ready for Recombination...";
    TransformPartons(inPartons, outHadrons, outPartons);    
  }
  else
  {
    VERBOSE(2)<<"There is no Parton ready for Recombination...";
  }
}


void Hadronization::Exec()
{
  VERBOSE(2)<<"Run Hadronization Exec...";
  VERBOSE(2)<<"Found "<<GetNumberOfTasks()<<" Hadronization Tasks/Modules Execute them ... ";

  //this->outHadrons = make_shared<vector<shared_ptr<Hadron>>>();
  //this->outPartons = make_shared<vector<shared_ptr<Parton>>>();

  DoHadronize();
 
}

void Hadronization::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(4)<<"In Hadronization::WriteTask";
  auto f = w.lock();
  if ( !f ) return;
  
  f->WriteComment("Hadronization module: "+GetId());

  if(GetHadrons().size()>0) {
    f->WriteComment("Final State Hadrons");
    int i=0;
    for ( auto& h : GetHadrons() ){
      f->WriteWhiteSpace("["+to_string(i)+"] H");
      f->Write( h );
      ++i;
    }
  }  else   {
    f->WriteComment("There are no Hadrons");
  }

}








}























