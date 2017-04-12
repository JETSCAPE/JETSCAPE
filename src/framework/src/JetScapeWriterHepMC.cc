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

  //This function dumps the particles in the JetScapeEvent container (that everyone has access to)
  GenEvent evt(Units::GEV,Units::MM);
  HepMC::GenVertex *vtx = new HepMC::GenVertex(HepMC::FourVector(0,0,0,0));
  //if(!(this->getParentTask()->getEvent())){
  //  INFO << "JetScape Event not-initialized - exiting...";
  //}
  //vector<Parton> pList = this->getParentTask()->getEvent()->getPartonCollection();
  
  weak_ptr<HardProcess> hproc = JetScapeSignalManager::Instance()->GetHardProcessPointer();
  vector<shared_ptr<Parton> > pList = hproc.lock()->GetPartonList();

  INFO << "Parton list size to dump: " << pList.size();
  for(unsigned int ipart=0; ipart<pList.size(); ipart++){
      HepMC::GenParticle *p1 = castParticleToHepMC(pList.at(ipart));
      INFO << "   parton " << ipart << " pt " << p1->momentum().perp();
      vtx->add_particle_out(p1);
  }
  evt.add_vertex(vtx);
  write_event(evt);
}

void JetScapeWriterHepMC::WriteEvent(weak_ptr<PartonShower> ps){

    //This function dumps the particles in a specific parton shower to a file
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
