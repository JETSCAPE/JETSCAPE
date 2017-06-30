// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "InitialState.h"
#include "JetScapeWriter.h"
#include <iostream>

namespace Jetscape {

InitialState::InitialState(){
    collEnergy = 13000;
    temperature = 175; //MeV? Get the units from HepMC?
    FourVector v(0,0,0,0);
    initialVtx.set_location(v);
}

InitialState::InitialState(double collE, Vertex vtx)
{
  collEnergy = collE;
  initialVtx = vtx;

}

InitialState::~InitialState()
{
}

void InitialState::Init(){
    
    //Set Pythia (or PGun) collision Energy?
    
}

void InitialState::Exec(){
    // Do whatever is needed to figure out the internal temp...
    
}

void InitialState::Clear(){
}

void InitialState::Write(weak_ptr<JetScapeWriter> w){

    //Write out the original vertex so the writer can keep track of it...
    w.lock()->Write(make_shared<Vertex>(initialVtx));
}

} // end namespace Jetscape
