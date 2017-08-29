#include<iostream>

#include "PartonPrinter.h"
#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.hpp"
#include "JetScapeLogger.h"

using namespace std;

namespace Jetscape {

PartonPrinter::PartonPrinter()
{
}

PartonPrinter::~PartonPrinter()
{
}

void PartonPrinter::Init()
{
  this->SetId("Printer");
}

void PartonPrinter::Exec()
{
  INFO<<"$$$$$Run PartonPrinter: print shower from event # "<<GetCurrentEvent()<<" ...";
}

void PartonPrinter::GetFinalPartons(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons)
{
  if(pShower)
  {
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size(); ++ipart)
    {
      fPartons.push_back( pShower.get()->GetFinalPartons().at(ipart));
      cout << "############### FINAL PARTON IN THE VECTOR NUMBER : " << ipart << " is " << fPartons.at(ipart) << "###################\n";
    }

  }
}

void PartonPrinter::Clear()
{
 // Clear();
}

void PartonPrinter::GetPartonsAtTime(shared_ptr<PartonShower> pShower,  vector<shared_ptr<Parton>>& fPartons, double time)
{
}

void PartonPrinter::PrintFinalPartons(shared_ptr<PartonShower> pShower)
{
 if(pShower)
  {
    cout << "@@@@@@@@@@@@@@@@ Number of FINAL Partons: " << pShower.get()->GetFinalPartons().size() << "@@@@@@@@@@@@@@@@@@@@@@@\n";
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size(); ++ipart)
    {
      cout << "@@@@@@@@@@@@@@@@ FINAL PARTON NUMBER : " << ipart << " is " <<  pShower.get()->GetFinalPartons().at(ipart) << "@@@@@@@@@@@\n";
    }
  }
  else
  {
     cout << "@@@@@@@@@@@@@@@@@ THERE IS NO SHOWER NOW \n";
  }
 
}

} // end namespace Jetscape




























