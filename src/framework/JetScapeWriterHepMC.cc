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
#include <GTL/topsort.h>

using HepMC::Units;

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
    // TODO: GeV seems right, but I don't think we actually measure in mm
    // Should multiply all lengths by 1e-12 probably
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
    
    // The following works after a fashion
    // I.e. it gives a healthy event as read by HepMC3_fileIO_example.exe
    // However, it also doubles and triples partons, so try a different approach
#if FALSE
    // Collect all vertices, remember their origin
    map <int, HepMC::GenVertex*> IdVertexMap;
    PartonShower::node_iterator nIt,nEnd;
    for ( nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd; ++nIt){
      IdVertexMap [ nIt->id() ] = castVtxToHepMC( pShower->GetVertex( *nIt ) ) ;
      IdVertexMap [ nIt->id() ]->set_status( nIt->id() ); // Also remember it in the output
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

    // Every vertex needs an incoming and an outgoing particle, let's attach those
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
#endif

    // Need topological order, see
    // https://hepmc.web.cern.ch/hepmc/differences.html
    // That means if parton p1 comes into vertex v, and p2 goes out of v,
    // then p1 has to be created (bestowed an id) before p2
    // Take inspiration from hepmc3.0.0/interfaces/pythia8/src/Pythia8ToHepMC3.cc
    // But pythia showers are different from our existing graph structure,
    // So instead try to modify the first attempt to respect top. order
    // and don't create vertices and particles more than once
    // What we have:
    // - every node knows it's parent(s) and daughter(s)
    // - HepMC vertices don't need to be created in order
    // However:
    // - even though our graph is directed, we don't know for certain
    //   that the gtl node iterator delivers nodes in a topological order
    // But we can at least catch problems and throw errors
    // --> Strategy
    // For every node:
    // 1. check whether it has a mother. Every vertex needs an incoming and an outgoing particle.
    //    a. if not create a dummy one (clone of an outgoer?)
    //    b. for all mothers:    
    //         i) has it already

    // Let's try GTL's topsort
    // 1. Check that our graph is sane
    if ( !pShower->is_acyclic() ) throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Graph is not acyclic.");

    // 2.
    topsort topsortsearch;
    topsortsearch.scan_whole_graph(true);
    topsortsearch.start_node();// defaulted to first node
    topsortsearch.run(*pShower);

    // topsort::topsort_iterator nEnd = topsortsearch.top_order_end();
    auto nEnd = topsortsearch.top_order_end(); // this is a topsort::topsort_iterator

    // Need to keep track of already created ones
    map< int , GenParticlePtr > CreatedPartons;
      
    // probably unnecessary, used for consistency checks
    map< int , bool > vused;
    bool foundRoot=false;
    for ( auto nIt = topsortsearch.top_order_begin(); nIt != nEnd; ++nIt){
      // cout << *nIt << "  " << nIt->indeg() << "  " << nIt->outdeg() << endl;

      // Should be the only time we see this node.
      if ( vused.find ( nIt->id() ) != vused.end() ) throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Reusing a vertex.");
      vused[ nIt->id() ]=true;

      cout << pShower->GetVertex( *nIt )->x_in().t() << endl;

      // 1. Create a new vertex.
      // --------------------------------------------
      auto v = castVtxToHepMC( pShower->GetVertex( *nIt ) );

      // 2. Incoming edges?
      // --------------------------------------------
      // 2.1: No. Need to create one (thanks, HepMC)
      if ( nIt->indeg() == 0 ){
	if (foundRoot){ // Should only happen once unless we merged showers
	  WARN << "Found a second root! Should only happen if we merged showers. "
	       << " If that's the case, please comment out the following throw and recompile";
	  throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Found a second root.");
	}
	foundRoot=true;

	// Some choices here. In general, could just attach to a dummy.
	// Maybe better: The root should only have exactly one daughter, let's check and clone
	if ( nIt->outdeg() != 1 ){
	  throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Root has illegal number of daughters");
	}
	auto out   = pShower->GetParton( *(nIt->out_edges_begin()) );
	auto hepin = castParticleToHepMC( out );
	v->add_particle_in( hepin );
      }	      

      // 2.2: Have incoming edges.
      if ( nIt->indeg() > 0 ){
	//  In the current framework, it should only be one.
	//  So we will catch anything more but provide a mechanism that should work anyway.
	if (nIt->indeg() > 1 ){
	  WARN << "Found more than one mother parton! Should only happen if we added medium particles. "
	       << "The code should work, but proceed with caution";
	}

	auto inIt  = nIt->in_edges_begin();
	auto inEnd = nIt->in_edges_end();
	for ( /* nop */; inIt != inEnd; ++inIt ){

	  auto phepin = CreatedPartons.find ( inIt->id() );
	  if ( phepin != CreatedPartons.end() ){
	    // We should already have one!
	    v->add_particle_in( phepin->second );
	  } else {
	    WARN << "Incoming particle out of nowhere. This could maybe happen if we pick up medium particles "
		 << " but is probably a topsort problem. Try using the code after this throw but be very careful.";
	    throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Incoming particle out of nowhere.");

	    auto in   = pShower->GetParton( *inIt );
	    auto hepin = castParticleToHepMC( in );
	    CreatedPartons [ inIt->id() ] = hepin;
	    v->add_particle_in( hepin );	    
	  }
	  // auto out   = pShower->GetParton( *(nIt->out_edges_begin()) );
	  // auto hepin = castParticleToHepMC( out );
	  // v->add_particle_in( hepin );
	}
      }	      

      // 3. Outgoing edges?
      // --------------------------------------------
      // 3.1: No. Need to create one (thanks, HepMC)
      if ( nIt->outdeg() == 0 ){
      	// Some choices here. In general, could just attach to a dummy.
	// That may be the best option, but it is hard to imagine a physics scenario
	// where a final vertex should have more than one parent
	// So let's clone the incoming guy
      	if ( nIt->indeg() != 1 ){
      	  throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Need exactly one parent to clone final state partons.");
      	}
      	auto in     = pShower->GetParton( *(nIt->in_edges_begin()) );
      	auto hepout = castParticleToHepMC( in );
      	v->add_particle_out( hepout );
	// Note that we do not register this particle. Since it's pointing nowhere it can never be reused.
      }	      

      // 3.2: Otherwise use it and register it
      if ( nIt->outdeg() > 0 ){
	auto outIt  = nIt->out_edges_begin();
	auto outEnd = nIt->out_edges_end();
	for ( /* nop */; outIt != outEnd; ++outIt ){
	  if ( CreatedPartons.find ( outIt->id() ) != CreatedPartons.end() ){
	    throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Trying to recreate a preexisting GenParticle.");
	  }
	  auto out = pShower->GetParton( *outIt );
	  auto hepout = castParticleToHepMC( out );		  
	  CreatedPartons [ outIt->id() ] = hepout;
	  v->add_particle_out( hepout );
	}
      }	      


      // if ( nIt->indeg() == 0 )
      vertices.push_back ( v );
      cout.precision(3);
      cout << v->position().t() << endl;
      
    }
    // for ( auto v : vused ) cout << v.first << "  " << v.second << endl;
    // throw std::runtime_error("Done");
	
    // // 1. Collect all vertices, remember their id
    // // Creation order of vertices is irrelevant to HepMC
    // map <int, HepMC::GenVertexPtr> IdVertexMap;
    // map <int, HepMC::GenParticlePtr> IdParticleMap;
    // PartonShower::node_iterator nIt,nEnd;
    // for ( nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd; ++nIt){
    //   IdVertexMap [ nIt->id() ] = castVtxToHepMC( pShower->GetVertex( *nIt ) ) ;
    //   IdVertexMap [ nIt->id() ]->set_status( nIt->id() ); // Also remember id in the output

    //   // Now every particle
    // }


    // // 1. Fill particle information 
    // std::vector<GenParticlePtr> hepevt_particles;
    // PartonShower::edge_iterator eIt,eEnd;
    // for ( eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd; ++eIt){
    //   hepevt_particles.push_back( castParticleToHepMC( pShower->GetParton(*eIt) ) );

    //   // not sure how to deal with mass
    //   // hepevt_particles.back()->set_generated_mass( pyev[i].m() );
    // }


    // // 2. Fill vertex information
    // std::vector<GenVertexPtr> vertex_cache;
    

    // // Create a particle instance for each entry and fill a map, and 
    // // a vector which maps from the particle index to the GenParticle address.
    // std::vector<GenParticle*> hepevt_particles;
    // PartonShower::edge_iterator eIt,eEnd;
    // for ( eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd; ++eIt){      
    //   // Fill the particle.
    //   HepMC::GenParticle* p = castParticleToHepMC( pShower->GetParton(*eIt) );
    //   hepevt_particles.push_back( p );
    //   // hepevt_particles[i]->suggest_barcode(i);
    //   // hepevt_particles[i]->set_generated_mass( momFac * pyev[i].m() );

    //   // // Colour flow uses index 1 and 2.
    //   // int colType = pyev[i].colType();
    //   // if (colType ==  1 || colType == 2)
    //   // 	hepevt_particles[i]->set_flow(1, pyev[i].col());
    //   // if (colType == -1 || colType == 2)
    //   // 	hepevt_particles[i]->set_flow(2, pyev[i].acol());
    // }

    // // Could set beam particles here - not trivial in AA?
    // // evt->set_beam_particles( hepevt_particles[1], hepevt_particles[2] );

    // // 3. Loop over particles AGAIN, this time creating vertices from our own vertices.
    // // The HEPEVT pointers are bi-directional, so gives decay vertices as well.
    
    
    
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
