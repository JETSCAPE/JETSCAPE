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

#include "TrentoCollision.h"

#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include "nucleus.h"

namespace Jetscape {

// Helper functions for Collider ctor.
namespace {
// Create one nucleus from the configuration.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  const auto& nucleon_width = var_map["nucleon-width"] .as<double>();
  return Nucleus::create(species, nucleon_width);
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

// void write_stream(std::ostream& os, int width,
//     int num, double impact_param, const Event& event) {
//   using std::fixed;
//   using std::setprecision;
//   using std::setw;
//   using std::scientific;

//   // Write a nicely-formatted line of event properties.
//   os << setprecision(10)
//      << setw(width)            << num
//      << setw(15) << fixed      << impact_param
//      << setw(5)                << event.npart()
//      << setw(18) << scientific << event.multiplicity()
//      << fixed;

//   for (const auto& ecc : event.eccentricity())
//     os << setw(14)             << ecc.second;

//   for (const auto& phi_n : event.participant_plane())
//     os << setw(14)             << phi_n.second;

//   // Write the mass center (x, y)
//   os << setw(14) << event.mass_center_index().first;
//   os << setw(14) << event.mass_center_index().second;

//   os << '\n';
// }


} // end unnamedspace


TrentoCollision::TrentoCollision(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_profile_(var_map),
      nevents_(var_map["number-events"].as<int>()),
      bmin_(var_map["b-min"].as<double>()),
      bmax_(determine_bmax(var_map, *nucleusA_, *nucleusB_, nucleon_profile_)),
      asymmetry_(determine_asym(*nucleusA_, *nucleusB_)),
      event_(var_map),
      output_(var_map) {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

// See header for explanation.
TrentoCollision::~TrentoCollision() = default;


void TrentoCollision::sample_(double smin, double smax){
    double b;
    while (true) {
        b = sample_impact_param();
        // Pass the prepared nuclei to the Event.  It computes the entropy profile
        // (thickness grid) and other event observables.
        event_.compute(*nucleusA_, *nucleusB_, nucleon_profile_);

        double total_entropy = event_.multiplicity();
        if (total_entropy >= smin && total_entropy <= smax) break;
    }

    info_.impact_parameter = b;
    info_.num_participant = event_.npart();
    info_.total_entropy = event_.multiplicity();
    info_.ecc = event_.eccentricity();
    info_.psi = event_.participant_plane();

    // Write the mass center (x, y)
    info_.xmid = event_.mass_center_index().first;
    info_.ymid = event_.mass_center_index().second;

    // write_stream(std::cout, 2, 1, b, event_);
}

double TrentoCollision::sample_impact_param() {
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
        collision = nucleon_profile_.participate(A, B) || collision;
      }
    }
  } while (!collision);

  return b;
}

} // end namespace Jetscape
