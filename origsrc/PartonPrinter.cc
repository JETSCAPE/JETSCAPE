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
  VERBOSE(2) <<"Run PartonPrinter: print shower from event # "<<GetCurrentEvent()<<" ...";
}

void PartonPrinter::GetFinalPartons(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons)
{
  if(pShower)
  {
    vector<shared_ptr<Parton>> vPin;
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size(); ++ipart)
    {
      fPartons.push_back( pShower.get()->GetFinalPartons().at(ipart));
        dist_output << ipart << " " << fPartons.at(ipart)->pid() << " " << fPartons.at(ipart)->e() << " " << fPartons.at(ipart)->px() << " " << fPartons.at(ipart)->py() << " " << fPartons.at(ipart)->pz() << endl;
      vPin.push_back( pShower.get()->GetFinalPartons().at(ipart));	
    }
    this->pFinals.push_back(vPin);
  }
}

void PartonPrinter::GetFinalPartons2(shared_ptr<PartonShower> pShower)
{
  if(pShower)
  {
/*
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size(); ++ipart)
    {
      this->pFinals.push_back( pShower.get()->GetFinalPartons());
      //cout << "############### FINAL PARTON IN THE VECTOR NUMBER : " << ipart << " is " << this->pFinals.at(ipart) << "###################\n";
    }
*/
    //this->pFinals.clear();
    this->pFinals.push_back( pShower.get()->GetFinalPartons());
  }
}

void PartonPrinter::Clear()
{
    this->pFinals.clear();
}

void PartonPrinter::GetPartonsAtTime(shared_ptr<PartonShower> pShower,  vector<shared_ptr<Parton>>& fPartons, double time)
{
}

} // end namespace Jetscape




























