/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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

    auto heavyion = make_shared<HepMC::GenHeavyIon>();
    // see https://gitlab.cern.ch/hepmc/HepMC3/blob/master/include/HepMC/GenHeavyIon.h
    if ( GetHeader().GetNpart() > -1 ){
      // Not clear what the difference is...
      heavyion->Ncoll_hard = GetHeader().GetNcoll();
      heavyion->Ncoll = GetHeader().GetNcoll();
    }
    if ( GetHeader().GetNcoll() > -1 ){
      // Hepmc separates into target and projectile.
      // Set one? Which? Both? half to each? setting projectile for now.
      // setting both might lead to weird problems when they get added up
      heavyion->Npart_proj = GetHeader().GetNpart();
    }
    if ( GetHeader().GetTotalEntropy() > -1 ){
      // nothing good in the HepMC standard. Something related to mulitplicity would work
    }

    if ( GetHeader().GetEventPlaneAngle() > -999 ){
      heavyion->event_plane_angle = GetHeader().GetEventPlaneAngle();
    }

    evt.set_heavy_ion( heavyion );  
  }
  
  void JetScapeWriterHepMC::WriteEvent() {    
    VERBOSE(1)<<"Run JetScapeWriterHepMC: Write event # "<<GetCurrentEvent();

    // Have collected all vertices now, add them to the event
    for( auto v : vertices )      evt.add_vertex( v );
    INFO << " found " << vertices.size() << " vertices in the list";
    
    write_event(evt);
    vertices.clear();
    hadronizationvertex=0;
  }
  
  //This function dumps the particles in a specific parton shower to the event
  void JetScapeWriterHepMC::Write(weak_ptr<PartonShower> ps){
    shared_ptr<PartonShower> pShower = ps.lock();
    if ( !pShower ) return;
    
    // Need topological order, see
    // https://hepmc.web.cern.ch/hepmc/differences.html
    // That means if parton p1 comes into vertex v, and p2 goes out of v,
    // then p1 has to be created (bestowed an id) before p2
    // Take inspiration from hepmc3.0.0/interfaces/pythia8/src/Pythia8ToHepMC3.cc
    // But pythia showers are different from our existing graph structure,
    // So instead try to modify the first attempt to respect top. order
    // and don't create vertices and particles more than once

    // Using GTL's topsort
    // 1. Check that our graph is sane
    if ( !pShower->is_acyclic() ) throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Graph is not acyclic.");

    // 2.
    topsort topsortsearch;
    topsortsearch.scan_whole_graph(true);
    topsortsearch.start_node(); // defaults to first node
    topsortsearch.run(*pShower);
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
	auto hepin = castPartonToHepMC( out );
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
		 << " but is probably a topsort problem. Try using the code after this throw() but be very careful.";
	    throw std::runtime_error("PROBLEM in JetScapeWriterHepMC: Incoming particle out of nowhere.");

	    auto in   = pShower->GetParton( *inIt );
	    auto hepin = castPartonToHepMC( in );
	    CreatedPartons [ inIt->id() ] = hepin;
	    v->add_particle_in( hepin );	    
	  }
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
      	auto hepout = castPartonToHepMC( in );
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
	  auto hepout = castPartonToHepMC( out );
	  CreatedPartons [ outIt->id() ] = hepout;
	  v->add_particle_out( hepout );
	}
      }	      
      
      vertices.push_back ( v );      
    }
  }

  void JetScapeWriterHepMC::Write(weak_ptr<Hadron> h)
  {
    auto hadron = h.lock();
    if ( !hadron ) return;

    // No clear source for most hadrons
    // Also, a graph with e.g. recombination hadrons would have loops,
    // (though the direction should still make it acyclic?)
    // Not sure how this is supposed to be done in HepMC3
    // Our solution: Attach all hadrons to one dedicated hadronization vertex.
    // Future option: Have separate shower and bulk vertices?

    // Create if it doesn't exist yet
    if ( !hadronizationvertex ) {
      // dummy position
      HepMC::FourVector vtxPosition( 0,0,0, 100 ); // set it to a late time...
      hadronizationvertex =  make_shared<GenVertex>(vtxPosition);

      // dummy mother -- could also maybe use the first/hardest shower initiator
      HepMC::FourVector pmom(0, 0, 0, 0);
      make_shared<GenParticle> (pmom, 0, 0);
      hadronizationvertex->add_particle_in( make_shared<GenParticle> (pmom, 0, 0) );

      vertices.push_back ( hadronizationvertex );
    }

    // now attach
    hadronizationvertex->add_particle_out( castHadronToHepMC( hadron ) );
  }

  void JetScapeWriterHepMC::Init()
  {
    if (GetActive()) {
      INFO<<"JetScape HepMC Writer initialized with output file = "<<GetOutputFileName();
    }
  }
  
  void JetScapeWriterHepMC::Exec(){
    // Nothing to do
  }

  // // NEVER use this!
  // // Can work with only one writer, but with a second one it gets called twice
  // void JetScapeWriterHepMC::WriteTask(weak_ptr<JetScapeWriter> w){
  //   //redirect - do the writing in WriteEvent, not in exec...
  //   cerr << "GETTING CALLED"<<endl;
  //   // WriteEvent();
  // }


} // end namespace Jetscape
