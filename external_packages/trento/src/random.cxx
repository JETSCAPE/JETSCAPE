// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include "random.h"

namespace trento { namespace random {

// Seed random number generator from hardware device.
Engine engine{std::random_device{}()};

}}  // namespace trento::random
