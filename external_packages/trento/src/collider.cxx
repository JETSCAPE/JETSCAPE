// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include "collider.h"

#include <cmath>
#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"
#include "nucleus.h"
#include <iostream>

namespace trento {

namespace {

// Helper functions for Collider ctor.

// Create one nucleus from the configuration.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  const auto& nucleon_dmin = var_map["nucleon-min-dist"].as<double>();
  const auto& nucleon_width = var_map["nucleon-width"].as<double>();
  return Nucleus::create(species, nucleon_width, nucleon_dmin);
}

// Determine the maximum impact parameter.  If the configuration contains a
// non-negative value for bmax, use it; otherwise, fall back to the minimum-bias
// default.
double determine_bmax(const VarMap& var_map,
    const Nucleus& A, const Nucleus& B, const NucleonProfile& profile) {
  auto bmax = var_map["b-max"].as<double>();
  if (bmax < 0.)
    bmax = A.radius() + B.radius() + profile.max_impact();
  return bmax;
}

// Determine the asymmetry parameter (Collider::asymmetry_) for a pair of
// nuclei.  It's just rA/(rA+rB), falling back to 1/2 if both radii are zero
// (i.e. for proton-proton).
double determine_asym(const Nucleus& A, const Nucleus& B) {
  double rA = A.radius();
  double rB = B.radius();
  double sum = rA + rB;
  if (sum < 0.1)
    return 0.5;
  else
    return rA/sum;
}

}  // unnamed namespace

// Lots of members to initialize...
// Several helper functions are defined above.
Collider::Collider(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_profile_(var_map),
      nevents_(var_map["number-events"].as<int>()),
      ntrys_(0),
      bmin_(var_map["b-min"].as<double>()),
      bmax_(determine_bmax(var_map, *nucleusA_, *nucleusB_, nucleon_profile_)),
      npartmin_(var_map["npart-min"].as<int>()),
      npartmax_(var_map["npart-max"].as<int>()),
      stotmin_(var_map["s-min"].as<double>()),
      stotmax_(var_map["s-max"].as<double>()),
      asymmetry_(determine_asym(*nucleusA_, *nucleusB_)),
      event_(var_map),
      output_(var_map),
      with_ncoll_(var_map["ncoll"].as<bool>())
{
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

// See header for explanation.
Collider::~Collider() = default;

void Collider::run_events() {
  // The main event loop.
  for (int n = 0; n < nevents_; ++n) {
	if (n%1000 == 0 && n!=0) std::cout<< "# "<< n << " events generated" << std::endl; 
    // Sampling the impact parameter also implicitly prepares the nuclei for
    // event computation, i.e. by sampling nucleon positions and participants.

    // WK: an extra do-while loop, sample events, until it meets the Npart or 
    // Entropy cut provided from command lines
    bool fullfil_Npart_cut=false, fullfil_Entropy_cut=false;
    double b;
    do{
    	b = sample_impact_param();
    	// Pass the prepared nuclei to the Event.  It computes the entropy profile
    	// (thickness grid) and other event observables.
    	event_.compute(*nucleusA_, *nucleusB_, nucleon_profile_);
        fullfil_Npart_cut = (npartmin_ < event_.npart()) 
								&& (event_.npart() <= npartmax_);
        fullfil_Entropy_cut = (stotmin_ < event_.multiplicity()) 
								&& (event_.multiplicity() <= stotmax_); 
	}while( (!fullfil_Npart_cut) || (!fullfil_Entropy_cut) );
    // Write event data.
    output_(n, b, event_);
	records this_event;
	this_event.i = n;
	this_event.b = b;
	this_event.npart = event_.npart();
	this_event.mult = event_.multiplicity();
	all_records_.push_back(this_event);
  }
  double cross_section = nevents_*M_PI*(bmax_*bmax_ - bmin_*bmin_)/ntrys_;
  double cross_section_err = cross_section/std::sqrt(1.*nevents_);
  //std::cout << "# cross-section = " << cross_section
  //			<< " +/- " << cross_section_err <<" [fm^2]" << std::endl; 
}

double Collider::sample_impact_param() {
  // Sample impact parameters until at least one nucleon-nucleon pair
  // participates.  The bool 'collision' keeps track -- it is effectively a
  // logical OR over all possible participant pairs.
  double b;
  bool collision = false;

  do {
    // Sample b from P(b)db = 2*pi*b.
    b = bmin_ + (bmax_ - bmin_) * std::sqrt(random::canonical<double>());

    // Offset each nucleus depending on the asymmetry parameter (see header).
    nucleusA_->sample_nucleons(asymmetry_ * b);
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);

    // Check each nucleon-nucleon pair.
    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
		bool AB_collide = nucleon_profile_.participate(A, B);

        if (with_ncoll_) {
			// WK: only init Ncoll and Ncoll density at the first binary collision:
			if (AB_collide && (!collision) ) event_.clear_TAB();
			// WK: to calculate binary collision denstiy, each collision 
			// contribute independently its Tpp. Therefore, if one pair collide, 
			// it calls the event object to accumulate Tpp to the Ncoll density
			// Ncoll density = Sum Tpp		
			if (AB_collide) event_.accumulate_TAB(A, B, nucleon_profile_);
		}

		// update collision flag
        collision = AB_collide || collision;
      }
    }
    ntrys_ ++;
  } while (!collision);

  return b;
}

}  // namespace trento
