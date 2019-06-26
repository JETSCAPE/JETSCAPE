// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/collider.h"

#include "catch.hpp"
#include "util.h"

#include <iostream>

#include "../src/nucleus.h"

using namespace trento;

TEST_CASE( "collider" ) {
  constexpr auto N = 5;

  auto var_map = make_var_map({
    {"number-events", N},
    {"quiet", false},
    {"random-seed", static_cast<int64_t>(-1)},
    {"projectile", std::vector<std::string>{"Pb", "Pb"}},
    {"b-min", 0.},
    {"b-max", -1.},
    {"normalization", 1.},
    {"reduced-thickness", 0.},
    {"grid-max", 9.},
    {"grid-step", 0.3},
    {"fluctuation", 1.},
    {"cross-section", 6.4},
    {"nucleon-width", 0.5},
    {"nucleon-min-dist", 0.},
  });

  std::vector<int> nevent, npart;
  std::vector<double> impact, mult;

  // run collider normally and save output
  {
    capture_stdout capture;

    Collider collider{var_map};
    collider.run_events();

    std::string line;
    while (std::getline(capture.stream, line)) {
      nevent.push_back(0);
      impact.push_back(0);
      npart.push_back(0);
      mult.push_back(0);
      std::istringstream(line) >> nevent.back()
                               >> impact.back()
                               >> npart.back()
                               >> mult.back();
    }
  }

  // event numbers should be an integer sequence from zero
  std::vector<int> sequence(N);
  std::iota(sequence.begin(), sequence.end(), 0);
  CHECK( nevent == sequence );

  // verify impact parameters are within min-bias range
  auto impact_max = 2*Nucleus::create("Pb")->radius() + 6*.5;
  CHECK( impact >= std::vector<double>(N, 0.) );
  CHECK( impact <= std::vector<double>(N, impact_max) );

  // verify all events have at least 2 participants and at most 416 for PbPb
  CHECK( npart >= std::vector<int>(N, 2) );
  CHECK( npart <= std::vector<int>(N, 416) );

  // nonzero multiplicity
  CHECK( mult >= std::vector<double>(N, 0.) );
}

TEST_CASE( "fixed impact parameter" ) {
  constexpr auto N = 5;
  constexpr auto bfixed = 4.;

  auto var_map = make_var_map({
    {"number-events", N},
    {"quiet", false},
    {"random-seed", static_cast<int64_t>(-1)},
    {"projectile", std::vector<std::string>{"Au", "Au"}},
    {"b-min", bfixed},
    {"b-max", bfixed},
    {"normalization", 1.},
    {"reduced-thickness", 0.},
    {"grid-max", 9.},
    {"grid-step", 0.3},
    {"fluctuation", 1.},
    {"cross-section", 6.4},
    {"nucleon-width", 0.5},
    {"nucleon-min-dist", 0.2},
  });

  std::vector<double> impact;

  {
    Collider collider{var_map};

    capture_stdout capture;
    collider.run_events();

    std::string line;
    while (std::getline(capture.stream, line)) {
      double x;
      std::istringstream(line) >> x >> x;
      impact.push_back(x);
    }
  }

  // all events have the specified impact parameter
  CHECK( std::all_of(impact.cbegin(), impact.cend(),
    [&bfixed](double b) { return b == Approx(bfixed); }) );
}

TEST_CASE( "random seed" ) {
  std::vector<std::string> output(5);

  const auto seed = static_cast<int64_t>(std::random_device{}());

  // run several collider batches with the same seed
  std::generate(output.begin(), output.end(),
    [&seed]() {
      Collider collider{make_var_map({
        {"number-events", 3},
        {"quiet", false},
        {"random-seed", seed},
        {"projectile", std::vector<std::string>{"p", "U"}},
        {"b-min", 0.},
        {"b-max", -1.},
        {"normalization", 1.},
        {"reduced-thickness", 0.},
        {"grid-max", 9.},
        {"grid-step", 0.3},
        {"fluctuation", 1.},
        {"cross-section", 6.4},
        {"nucleon-width", 0.5},
        {"nucleon-min-dist", 0.4},
      })};

      capture_stdout capture;
      collider.run_events();

      return capture.stream.str();
    });

  // all collider batches are identical
  CHECK( std::all_of(output.cbegin(), output.cend(),
    [&output](const std::string& s) { return s == output.front(); }) );
}
