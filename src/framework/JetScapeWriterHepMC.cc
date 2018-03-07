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

namespace Jetscape {

JetScapeWriterHepMC::~JetScapeWriterHepMC()
{
  if (GetActive())
      Close();
}

void JetScapeWriterHepMC::Write(weak_ptr<Vertex> v){ 
    //Setting vertex from initial state...
    vertices.push_back(castVtxToHepMC(v.lock()));
    vertexFlag = true;
}


void JetScapeWriterHepMC::WriteEvent()
{
  INFO<< GetCurrentEvent() << " in HepMC ... ";

  //This function dumps the particles in the JetScapeEvent container (that everyone has access to)
  GenEvent evt(Units::GEV,Units::MM);
  
  INFO << " found " << vertices.size() << " vertices in the list...";
  
  for(unsigned int ivtx=0; ivtx<vertices.size(); ivtx++){
    evt.add_vertex(vertices.at(ivtx));
  }
  /*HepMC::GenVertex *vtx;
  if(vertexFlag) vtx = &initStateVtx;
  else vtx = new HepMC::GenVertex(HepMC::FourVector(0,0,0,0));
  
  weak_ptr<HardProcess> hproc = JetScapeSignalManager::Instance()->GetHardProcessPointer();
  vector<shared_ptr<Parton> > pList = hproc.lock()->GetPartonList();

  INFO << "Parton list size to dump: " << pList.size();
  for(unsigned int ipart=0; ipart<pList.size(); ipart++){
      HepMC::GenParticle *p1 = castParticleToHepMC(pList.at(ipart));
      INFO << "   parton " << ipart << " pt " << p1->momentum().perp();
      vtx->add_particle_out(p1);
  }
  evt.add_vertex(vtx);*/
  write_event(evt);
  vertices.clear();
}

void JetScapeWriterHepMC::Write(weak_ptr<PartonShower> ps){

    //This function dumps the particles in a specific parton shower to a file
    shared_ptr<PartonShower> pShower = ps.lock();

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
        //evt.add_vertex(vtx);
        vertices.push_back(vtx);
    }
    //Print::content(evt);
    //write_event(evt);
}

void JetScapeWriterHepMC::Init()
{
    vertexFlag = false;
    if (GetActive())
     {
       INFO<<"JetScape HepMC Writer initialized with output file = "<<GetOutputFileName();
     }
}

void JetScapeWriterHepMC::Exec()
{
}

void JetScapeWriterHepMC::WriteTask(weak_ptr<JetScapeWriter> w){
    //redirect - do the writing in Write, not in exec...
    WriteEvent();
}


} // end namespace Jetscape
