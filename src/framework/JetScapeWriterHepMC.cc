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
  
  JetScapeWriterHepMC::~JetScapeWriterHepMC(){
    if (GetActive())      Close();
  }
  
  // // We should never accept anything other than a full shower
  // void JetScapeWriterHepMC::Write(weak_ptr<Vertex> v){
  //   //Setting vertex from initial state
  //   vertices.push_back(castVtxToHepMC(v.lock()));
  // }
  
  void JetScapeWriterHepMC::WriteHeaderToFile() {
    // Create event here - not actually writing
    evt = GenEvent(Units::GEV,Units::MM);
    
    // Expects pb, pythia delivers mb
    auto xsec = make_shared<HepMC::GenCrossSection>();
    xsec->set_cross_section( GetHeader().GetSigmaGen() * 1e9, 0);
    xsec->set_cross_section( GetHeader().GetSigmaGen() * 1e9, GetHeader().GetSigmaErr() * 1e9);
    evt.set_cross_section( xsec );
    evt.weights().push_back( GetHeader().GetEventWeight() );
  }

  void JetScapeWriterHepMC::WriteEvent() {
    
    INFO<< GetCurrentEvent() << " in HepMC ... ";
    
    // Now write the shower
    INFO << " found " << vertices.size() << " vertices in the list...";
    
    // Have collected all vertices now, add them to the event
    for( auto v : vertices )      evt.add_vertex( v );
    
    write_event(evt);
    vertices.clear();
  }
  
  //This function dumps the particles in a specific parton shower to the event
  void JetScapeWriterHepMC::Write(weak_ptr<PartonShower> ps){
    shared_ptr<PartonShower> pShower = ps.lock();
    if ( !pShower ) return;
    
    // Collect all vertices, remember their origin
    map <int, HepMC::GenVertex*> IdVertexMap;
    PartonShower::node_iterator nIt,nEnd;
    for ( nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd; ++nIt){
      IdVertexMap [nIt->id() ] = castVtxToHepMC( pShower->GetVertex( *nIt ) ) ;
      IdVertexMap [nIt->id() ]->set_status( nIt->id() ); // Also remember it in the output
    }

    // Connect them via partons
    PartonShower::edge_iterator eIt,eEnd;
    for ( eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd; ++eIt){
      // cast edge
      HepMC::GenParticle* p = castParticleToHepMC( pShower->GetParton(*eIt));

      // where it's from
      IdVertexMap [eIt->source().id()]->add_particle_out( p );
      // where it points to
      IdVertexMap [eIt->target().id()]->add_particle_in( p );
    }

    // I _think_ every vertex needs an incoming and an outgoing particle, let's attach those
    for ( auto& indval : IdVertexMap ){
      auto& vp = indval.second;
      // no incoming, copy outgoing
      if ( vp->particles_in_size() == 0 ){
	// this should make a copy
	// Normally, only shallow copies are made which messes everything up...
	auto pp = *vp->particles_out_const_begin();
	auto p = make_shared<HepMC::GenParticle>( pp->momentum(), pp->pdg_id(), pp->status() );
	vp->add_particle_in( p );
      }
      // // no outgoing, copy incoming
      if ( vp->particles_out_size() == 0 ){
	// this should make a copy
	// Normally, only shallow copies are made which messes everything up...
	auto pp = *vp->particles_in_const_begin();
	auto p = make_shared<HepMC::GenParticle>( pp->momentum(), pp->pdg_id(), pp->status() );
      	vp->add_particle_out( p );
      }
    }

    // copy to list of vertices
    // there's room for optimization - but note that id is not unique across showers
    // vertices.clear();    // DEBUG ONLY -- keep only one shower
    for ( auto& indval : IdVertexMap ){
      vertices.push_back ( indval.second );
    }

    // //do all the vertices and link all the particles
    // for(unsigned int ivtx=0; ivtx<pShower->GetNumberOfVertices(); ivtx++){
    //   node vtxNode = pShower->GetNodeAt(ivtx);
    //   HepMC::GenVertex *vtx = castVtxToHepMC(pShower->GetVertex(vtxNode));
    //   node::in_edges_iterator inPartonsIt;
    //   for(inPartonsIt=vtxNode.in_edges_begin(); inPartonsIt!=vtxNode.in_edges_end(); inPartonsIt++){
    // 	HepMC::GenParticle *p1 = castParticleToHepMC(pShower->GetParton(*inPartonsIt));
    // 	vtx->add_particle_in(p1);
    //   }
    //   node::out_edges_iterator outPartonsIt;
    //   for(outPartonsIt=vtxNode.out_edges_begin(); outPartonsIt!=vtxNode.out_edges_end(); outPartonsIt++){
    // 	HepMC::GenParticle *p2 = castParticleToHepMC(pShower->GetParton(*outPartonsIt));
    // 	vtx->add_particle_out(p2);
    //   }
    //   //evt.add_vertex(vtx); // not sure why this doesn't work
    //   vertices.push_back(vtx);
    // }

    //Print::content(evt);
  }

  void JetScapeWriterHepMC::Init()
  {
    if (GetActive()) {
      INFO<<"JetScape HepMC Writer initialized with output file = "<<GetOutputFileName();
    }
  }
  
  void JetScapeWriterHepMC::Exec(){
    // Too early for HepMC, do the writing in Write, WriteTask
    // Create basic event structure
    // WriteEvent();    
  }
  
  void JetScapeWriterHepMC::WriteTask(weak_ptr<JetScapeWriter> w){
    //redirect - do the writing in WriteEvent, not in exec...
    WriteEvent();
  }


} // end namespace Jetscape
