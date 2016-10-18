#include "interface/Pythia8.h"

Pythia8::Pythia8(Parameter parameter_list):
  pset = parameter_list
{}

Pythia8::~Pythia8()
{}

void Pythia8::generatePartons(HepMC::GenEvent & event)
{
  // set up Pythia first with parameters

  // then run 1st step of Pythia

  // then add the hard partons to the list of particles to the HepMC event.
}
