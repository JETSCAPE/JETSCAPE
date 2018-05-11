// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef COLLIDER_H
#define COLLIDER_H

#include <memory>

#include "fwd_decl.h"
#include "event.h"
#include "nucleon.h"
#include "output.h"

namespace trento {

/// \rst
/// Orchestrates event generation and output.  Owns instances of several other
/// TRENTO classes and knows how they work together.  Responsible for sampling
/// impact parameters.  After instantiation, call ``run_events()`` to do
/// everything.
///
/// Example::
///
///   Collider collider{var_map};
///   collider.run_events();
///
/// \endrst
class Collider {
 public:
  /// Instantiate from the configuration.
  explicit Collider(const VarMap& var_map);

  /// Declare a destructor to properly handle the std::unique_ptr<Nucleus>
  /// members.  At this point in the code, Nucleus is an incomplete type so the
  /// compiler does not know how to delete it.  Therefore the destructor is
  /// defined (as default) in the implementation file, at which point Nucleus is
  /// fully defined.  See e.g. Item 22 of Effective Modern C++ by Scott Meyers
  /// and this stackoverflow answer: http://stackoverflow.com/a/6089065.
  ~Collider();

  /// Run events and output.
  void run_events();

 private:
  // Most of these are pretty self-explanatory...

  /// Sample a min-bias impact parameter within the set range.
  double sample_impact_param();

  /// Pair of nucleus projectiles.
  std::unique_ptr<Nucleus> nucleusA_, nucleusB_;

  /// The nucleon profile instance.
  NucleonProfile nucleon_profile_;

  /// Number of events to run.
  const int nevents_;

  /// Minimum and maximum impact parameter.
  const double bmin_, bmax_;

  /// Parameterizes the degree of asymmetry between the two projectiles.  Used
  /// to apportion the total impact parameter to each projectile so that the
  /// resulting overlap is approximately centered.  Given the two nuclear radii
  /// rA and rB (see Nucleus::radius()), the parameter is
  ///
  ///   asymmetry = rA / (rA + rB)
  ///
  /// and the impact parameter b is apportioned
  ///
  ///   offset_A = asymmetry * b
  ///   offset_B = (asymmetry-1) * b
  ///
  /// So for example, two identical projectiles would have asymmetry = 0.5, and
  /// each would be offset by half the impact parameter.  A proton-lead system
  /// would have asymmetry = 1, meaning the lead nucleus would be offset by the
  /// entire impact parameter and the proton would not be offset at all.
  const double asymmetry_;

  /// The event instance.
  Event event_;

  /// The output instance.
  Output output_;
};

}  // namespace trento

#endif  // COLLIDER_H
