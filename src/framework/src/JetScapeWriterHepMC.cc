// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeWriterHepMC.h"
#include "JetScapeLogger.h"
#include "HardProcess.h"
#include "JetScapeSignalManager.h"
#include "GTL/node.h"

JetScapeWriterHepMC::~JetScapeWriterHepMC()
{
  if (GetActive())
      Close();
}

void JetScapeWriterHepMC::WriteEvent()
{
  INFO<< GetCurrentEvent() << " in HepMC ... ";
  
  //This just writes an empty event since it doesn't have access to anything interesting...
  //GenEvent evt(Units::GEV,Units::MM);
  //Print::content(evt);
  //write_event(evt);
}

void JetScapeWriterHepMC::WriteEvent(weak_ptr<PartonShower> ps){

    shared_ptr<PartonShower> pShower = ps.lock();
    GenEvent evt(Units::GEV,Units::MM);

    //do all the vertices and link all the particles
    for(unsigned int ivtx=0; ivtx<pShower->GetNumberOfVertices(); ivtx++){
        node vtxNode = pShower->GetNodeAt(ivtx);
        HepMC::GenVertex *vtx = castVtxToHepMC(pShower->GetVertex(vtxNode));
        node::in_edges_iterator inPartonsIt;
        for(inPartonsIt=vtxNode.in_edges_begin(); inPartonsIt!=vtxNode.in_edges_end(); inPartonsIt++){
            HepMC::GenParticle *p1 = castParticleToHepMC(pShower->GetParton(*inPartonsIt));
            vtx->add_particle_in(p1);
        }
        node::out_edges_iterator outPartonsIt;
        for(outPartonsIt=vtxNode.out_edges_begin(); outPartonsIt!=vtxNode.out_edges_end(); outPartonsIt++){
            HepMC::GenParticle *p2 = castParticleToHepMC(pShower->GetParton(*outPartonsIt));
            vtx->add_particle_out(p2);
        }
        evt.add_vertex(vtx);
    }
    //Print::content(evt);
    write_event(evt);
}

void JetScapeWriterHepMC::Init()
{
   if (GetActive())
     {
       INFO<<"JetScape HepMC Writer initialized with output file = "<<GetOutputFileName();
     }
}

void JetScapeWriterHepMC::Exec()
{
  if (GetActive())
    WriteEvent();
}
