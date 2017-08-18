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

vector<shared_ptr<Parton>> PartonPrinter::GetFinalPartons(shared_ptr<PartonShower> pShower)
{
  return pShower->GetFinalPartons();
}

void PartonPrinter::Clear()
{
 // Clear();
}

vector<shared_ptr<Parton>> PartonPrinter::GetPartonsAtTime(shared_ptr<PartonShower> pShower, double time)
{
  if(partons.size()==0)
  {
    graph::edge_iterator eIt, eEnd;
    for(eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt!=eEnd; ++eIt)
    {
      //if(pShower->GetParton(*eIt).get()->form_time()==time)
      {
	this->partons.push_back(pShower->GetParton(*eIt));
      }
    }

/*    for(unsigned int ipart=0; ipart<this->pShower->GetNumberOfPartons(); ipart++)
    {
      if(this->pShower->GetPartons()[ipart].form_time()==time)
      {
	this->partons.push_back(this->pShower->GetPartons()[ipart]);
      }	
    }
*/
    return this->partons;
  }
  else
  {
    return this->partons;
  }	
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




























