// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include <boost/math/constants/constants.hpp>

#include "fast_exp.h"
#include "fwd_decl.h"
#include "random.h"

namespace trento {

class Nucleon;

/// \rst
/// Encapsulates properties shared by all nucleons: transverse thickness
/// profile, cross section, fluctuations.  Responsible for sampling
/// nucleon-nucleon participation with given `\sigma_{NN}`.
/// \endrst
class NucleonProfile {
 public:
  /// Instantiate from the configuration.
  explicit NucleonProfile(const VarMap& var_map);

  /// The radius at which the nucleon profile is truncated.
  double radius() const;

  /// The maximum impact parameter for participation.
  double max_impact() const;

  /// Randomly fluctuate the profile.  Should be called prior to evaluating the
  /// thickness function for a new nucleon.
  void fluctuate();

  /// Compute the thickness function at a (squared) distance from the profile
  /// center.
  double thickness(double distance_sqr) const;

  /// WK: same as above, but without the Gamma fluctuation, 
  /// used in the calculation of binary collision density 
  double deterministic_thickness(double distance_sqr) const;

  /// WK: return Tpp given bpp^2
  double norm_Tpp(double bpp_sqr) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(Nucleon& A, Nucleon& B) const;

 private:
  /// Width of Gaussian thickness function.
  const double width_sqr_;

  /// Truncate the Gaussian at this radius.
  const double trunc_radius_sqr_;

  /// Maximum impact parameter for participants.
  const double max_impact_sqr_;

  /// Cache (-1/2w^2) for use in the thickness function exponential.
  /// Yes, this actually makes a speed difference...
  const double neg_one_div_two_width_sqr_;

  /// WK 1/4w^2
  const double neg_one_div_four_width_sqr_;

  /// WK 1/4pi
  const double one_div_four_pi_;

  /// Dimensionless parameter set to reproduce the inelastic nucleon-nucleon
  /// cross section \sigma_{NN}.  Calculated in constructor.
  const double cross_sec_param_;

  /// Fast exponential for calculating the thickness profile.
  const FastExp<double> fast_exp_;

  /// Fluctuation distribution.
  std::gamma_distribution<double> fluct_dist_;

  /// Thickness function prefactor = fluct/(2*pi*w^2)
  double prefactor_;
 
  /// bool variable to calcualte Ncoll
  bool with_ncoll_;
};

/// \rst
/// Represents a single nucleon.  Stores its position and whether or not it's a
/// participant.  These properties are globally readable, but can only be set
/// through ``Nucleus`` and ``NucleonProfile``.
/// \endrst
class Nucleon {
 public:
  /// Only a default constructor is necessary\---the class is designed to be
  /// constructed once and repeatedly updated.
  Nucleon() = default;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// The longitudinal \em z position.
  double z() const;

  /// Whether or not this nucleon is a participant.
  bool is_participant() const;

 private:
  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;

  /// The NucleonProfile samples participants so must be able to set
  /// participation status.
  friend bool NucleonProfile::participate(Nucleon&, Nucleon&) const;

  /// Set the position and reset participant status to false.
  void set_position(double x, double y, double z);

  /// Mark as a participant.
  void set_participant();

  /// Internal storage of the position.
  double x_, y_, z_;

  /// Internal storage of participant status.
  bool participant_;
};

// These functions are short, called very often, and account for a large
// fraction of the total computation time, so request inlining.

// Nucleon inline member functions

inline double Nucleon::x() const {
  return x_;
}

inline double Nucleon::y() const {
  return y_;
}

inline double Nucleon::z() const {
  return z_;
}

inline bool Nucleon::is_participant() const {
  return participant_;
}

inline void Nucleon::set_position(double x, double y, double z) {
  x_ = x;
  y_ = y;
  z_ = z;
  participant_ = false;
}

inline void Nucleon::set_participant() {
  participant_ = true;
}

// NucleonProfile inline member functions

inline double NucleonProfile::radius() const {
  return std::sqrt(trunc_radius_sqr_);
}

inline double NucleonProfile::max_impact() const {
  return std::sqrt(max_impact_sqr_);
}

inline void NucleonProfile::fluctuate() {
  prefactor_ = fluct_dist_(random::engine) *
     math::double_constants::one_div_two_pi / width_sqr_;
}

inline double NucleonProfile::thickness(double distance_sqr) const {
  if (distance_sqr > trunc_radius_sqr_)
    return 0.;
  return prefactor_ * fast_exp_(neg_one_div_two_width_sqr_*distance_sqr);
}

// WK
inline double NucleonProfile::deterministic_thickness(double distance_sqr) const {
  if (distance_sqr > trunc_radius_sqr_)
    return 0.;
  return math::double_constants::one_div_two_pi / width_sqr_ 
		* fast_exp_(neg_one_div_two_width_sqr_*distance_sqr);
}

// WK
inline double NucleonProfile::norm_Tpp(double bpp_sqr) const  {
  return one_div_four_pi_ / width_sqr_ 
		* fast_exp_(neg_one_div_four_width_sqr_*bpp_sqr);
}

inline bool NucleonProfile::participate(Nucleon& A, Nucleon& B) const {
  // If both nucleons are already participants, there's nothing to do, unless
  // in Ncoll mode
  if (A.is_participant() && B.is_participant() && (! with_ncoll_))
    return true;

  double dx = A.x() - B.x();
  double dy = A.y() - B.y();
  double distance_sqr = dx*dx + dy*dy;

  // Check if nucleons are out of range.
  if (distance_sqr > max_impact_sqr_)
    return false;

  // The probability is
  //   P = 1 - exp(...)
  // which we could sample as
  //   P > U
  // where U is a standard uniform (0, 1) random number.  We can also compute
  //   1 - P = exp(...)
  // and then sample
  //   (1 - P) > (1 - U)
  // or equivalently
  //   (1 - P) < U
  auto one_minus_prob = std::exp(
      -std::exp(cross_sec_param_ - .25*distance_sqr/width_sqr_));

  // Sample one random number and decide if this pair participates.
  if (one_minus_prob < random::canonical<double>()) {
    A.set_participant();
    B.set_participant();
    return true;
  }

  return false;
}

}  // namespace trento

#endif  // NUCLEON_H
