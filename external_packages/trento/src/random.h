// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef RANDOM_H
#define RANDOM_H

#include <limits>
#include <random>

#include <boost/math/constants/constants.hpp>

#include "fwd_decl.h"

namespace trento { namespace random {

/// Mersenne Twister engine, 64-bit preset.
using Engine = std::mt19937_64;

/// Global variable defined in \c random.cxx.
extern Engine engine;

/// Helper function to easily generate random numbers in [0, 1).
template <typename RealType = double>
inline RealType canonical() {
  return std::generate_canonical
           <RealType, std::numeric_limits<RealType>::digits>
           (engine);
}

/// Sample a spherical polar angle cos(theta).
template <typename RealType = double>
inline double cos_theta() {
  return 2 * canonical<RealType>() - 1;
}

/// Sample a spherical azimuthal angle phi.
template <typename RealType = double>
inline double phi() {
  return math::constants::two_pi<RealType>() * canonical<RealType>();
}

}}  // namespace trento::random

#endif  // RANDOM_H
