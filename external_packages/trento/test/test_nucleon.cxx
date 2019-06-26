// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleon.h"

#include <cmath>

#include "catch.hpp"
#include "util.h"

#include "../src/nucleus.h"
#include "../src/random.h"

using namespace trento;

TEST_CASE( "nucleon" ) {
  auto fluct = 1. + .5*random::canonical<>();
  auto xsec = 4. + 3.*random::canonical<>();
  auto width = .5 + .2*random::canonical<>();
  auto wsq = width*width;

  auto var_map = make_var_map({
      {"fluctuation",   fluct},
      {"cross-section", xsec},
      {"nucleon-width", width},
  });

  NucleonProfile profile{var_map};
  profile.fluctuate();

  // truncation radius
  auto R = profile.radius();
  CHECK( R == Approx(5*width) );

  // thickness function
  // check relative to zero
  auto tzero = profile.thickness(0.);
  CHECK( profile.thickness(wsq) == Approx(tzero*std::exp(-.5)) );

  // random point inside radius
  auto dsq = std::pow(R*random::canonical<>(), 2);
  CHECK( profile.thickness(dsq) == Approx(tzero*std::exp(-.5*dsq/wsq)).epsilon(1e-5).margin(1e-5) );

  // random point outside radius
  dsq = std::pow(R*(1+random::canonical<>()), 2);
  CHECK( profile.thickness(dsq) == 0. );

  // fluctuations
  // just check they have unit mean -- the rest is handled by the C++ impl.
  auto total = 0.;
  auto n = 1e6;
  for (auto i = 0; i < static_cast<int>(n); ++i) {
    profile.fluctuate();
    total += profile.thickness(0.) * (2*M_PI*wsq);
  }

  auto mean = total/n;
  CHECK( mean == Approx(1.).epsilon(.003) );

  // must use a Nucleus to set Nucleon position
  // Proton conveniently sets a deterministic position
  // a mock class would be better but this works fine
  Proton A{}, B{};
  A.sample_nucleons(0.);
  B.sample_nucleons(0.);
  auto& nA = *A.begin();
  auto& nB = *B.begin();
  CHECK( nA.x() == 0. );
  CHECK( nA.y() == 0. );
  CHECK( nA.z() == 0. );
  CHECK( !nA.is_participant() );

  // wait until the nucleons participate
  while (!profile.participate(nA, nB)) {}
  CHECK( nA.is_participant() );
  CHECK( nB.is_participant() );

  // resampling nucleons resets participant state
  A.sample_nucleons(0.);
  CHECK( !nA.is_participant() );

  // test cross section
  // min-bias impact params
  auto bmax = profile.max_impact();
  CHECK( bmax == Approx(6*width) );

  auto nev = 1e6;
  auto count = 0;
  for (auto i = 0; i < static_cast<int>(nev); ++i) {
    auto b = bmax * std::sqrt(random::canonical<>());
    A.sample_nucleons(.5*b);
    B.sample_nucleons(-.5*b);
    if (profile.participate(nA, nB))
      ++count;
  }

  auto xsec_mc = M_PI*bmax*bmax * static_cast<double>(count)/nev;

  // precision is better than this, but let's be conservative
  CHECK( xsec_mc == Approx(xsec).epsilon(.02) );

  // impact larger than max should never participate
  auto b = bmax + random::canonical<>();
  A.sample_nucleons(.5*b);
  B.sample_nucleons(-.5*b);
  CHECK( !profile.participate(nA, nB) );

  // very large fluctuation parameters mean no fluctuations
  auto no_fluct_var_map = make_var_map({
      {"fluctuation",   1e12},
      {"cross-section", xsec},
      {"nucleon-width", width},
  });

  NucleonProfile no_fluct_profile{no_fluct_var_map};
  no_fluct_profile.fluctuate();
  CHECK( no_fluct_profile.thickness(0) == Approx(1/(2*M_PI*wsq)) );

  CHECK_THROWS_AS([]() {
    // nucleon width too small
    auto bad_var_map = make_var_map({
        {"fluctuation",   1.},
        {"cross-section", 5.},
        {"nucleon-width", .1},
    });
    NucleonProfile bad_profile{bad_var_map};
  }(),
  std::domain_error);
}
