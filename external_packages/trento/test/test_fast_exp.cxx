// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/fast_exp.h"

#include "catch.hpp"

#include "../src/random.h"

using namespace trento;

TEST_CASE( "fast exponential" ) {
  double xmin = -3.2;
  double xmax = 2.7;
  std::size_t nsteps = 1000;
  double tolerance = .25*std::pow((xmax - xmin)/nsteps, 2);

  FastExp<double> fast_exp{xmin, xmax, nsteps};

  std::uniform_real_distribution<double> dist{xmin, xmax};

  double worst_err = 0.;
  for (int i = 0; i < 1000; ++i) {
    auto x = dist(random::engine);
    auto approx = fast_exp(x);
    auto exact = std::exp(x);
    auto err = std::fabs(approx-exact)/exact;
    worst_err = std::max(worst_err, err);
  }

  CHECK( worst_err < tolerance );

#ifndef NDEBUG
  CHECK_THROWS_AS( fast_exp(xmin - 1), std::out_of_range );
#endif
}
