// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/event.h"

#include "catch.hpp"
#include "util.h"

#include "../src/nucleus.h"
#include "../src/random.h"

using namespace trento;

TEST_CASE( "event" ) {
  // Check Event class results against equivalent (but slower) methods.

  // Repeat for p == 0, p < 0, p > 0.
  auto pplus = .5 + .49*random::canonical<>();
  auto pminus = -.5 + .49*random::canonical<>();
  for (auto p : {0., pplus, pminus }) {
    // Random physical params.
    auto norm = 1. + .5*random::canonical<>();
    auto xsec = 4. + 3.*random::canonical<>();
    auto nucleon_width = .5 + .2*random::canonical<>();

    // Effectively disable fluctuations to deterministically compute thickness.
    auto fluct = 1e12;

    // Coarse-ish grid.
    auto grid_max = 9.;
    auto grid_step = 0.3;
    auto grid_nsteps = 60;

    auto var_map = make_var_map({
        {"normalization", norm},
        {"reduced-thickness", p},
        {"grid-max", grid_max},
        {"grid-step", grid_step},
        {"fluctuation",   fluct},
        {"cross-section", xsec},
        {"nucleon-width", nucleon_width},
    });

    Event event{var_map};
    NucleonProfile profile{var_map};

    CHECK( event.reduced_thickness_grid().num_dimensions() == 2 );
    CHECK( static_cast<int>(event.reduced_thickness_grid().shape()[0]) == grid_nsteps );
    CHECK( static_cast<int>(event.reduced_thickness_grid().shape()[1]) == grid_nsteps );

    auto nucleusA = Nucleus::create("Pb");
    auto nucleusB = Nucleus::create("Pb");

    // Sample impact param, nucleons, and participants.
    auto b = 4.*std::sqrt(random::canonical<>());
    nucleusA->sample_nucleons(+.5*b);
    nucleusB->sample_nucleons(-.5*b);

    for (auto&& A : *nucleusA)
      for (auto&& B : *nucleusB)
        profile.participate(A, B);

    // Run a normal Event.
    event.compute(*nucleusA, *nucleusB, profile);

    // Verify npart.
    auto count_part = [](const Nucleus& nucleus) {
      auto npart = 0;
      for (const auto& n : nucleus)
        if (n.is_participant())
          ++npart;
      return npart;
    };
    auto npart = count_part(*nucleusA) + count_part(*nucleusB);
    CHECK( npart == event.npart() );

    // Compute TR grid the slow way -- switch the order of grid and nucleon loops.
    boost::multi_array<double, 2> TR{boost::extents[grid_nsteps][grid_nsteps]};

    auto thickness = [&profile](const Nucleus& nucleus, double x, double y) {
      auto t = 0.;
      for (const auto& n : nucleus) {
        if (n.is_participant()) {
          auto dx = n.x() - x;
          auto dy = n.y() - y;
          t += profile.thickness(dx*dx + dy*dy);
        }
      }
      return t;
    };

    auto gen_mean = [p](double a, double b) {
      if (std::abs(p) < 1e-12)
        return std::sqrt(a*b);
      else if (p < 0. && (a < 1e-12 || b < 1e-12))
        return 0.;
      else
        return std::pow(.5*(std::pow(a, p) + std::pow(b, p)), 1./p);
    };

    for (auto iy = 0; iy < grid_nsteps; ++iy) {
      for (auto ix = 0; ix < grid_nsteps; ++ix) {
        auto x = (ix + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto y = (iy + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto TA = thickness(*nucleusA, x, y);
        auto TB = thickness(*nucleusB, x, y);
        TR[iy][ix] = norm * gen_mean(TA, TB);
      }
    }

    // Verify multiplicity.
    auto mult = std::pow(2*grid_max/grid_nsteps, 2) *
      std::accumulate(TR.origin(), TR.origin() + TR.num_elements(), 0.);
    CHECK( mult == Approx(event.multiplicity()) );

    // Verify each grid element.
    auto all_correct = true;
    for (const auto* t1 = TR.origin(),
         * t2 = event.reduced_thickness_grid().origin();
         t1 != TR.origin() + TR.num_elements();
         ++t1, ++t2) {
      if (*t1 != Approx(*t2).epsilon(1e-5)) {
        all_correct = false;
        WARN( "TR mismatch: " << *t1 << " != " << *t2 );
        break;
      }
    }
    CHECK( all_correct );

    // Verify eccentricity.
    auto sum = 0., xcm = 0., ycm = 0.;
    for (auto iy = 0; iy < grid_nsteps; ++iy) {
      for (auto ix = 0; ix < grid_nsteps; ++ix) {
        auto& t = TR[iy][ix];
        sum += t;
        xcm += t*ix;
        ycm += t*iy;
      }
    }
    xcm /= sum;
    ycm /= sum;

    for (auto n = 2; n <= 5; ++n) {
      auto real = 0., imag = 0., weight = 0.;
      for (auto iy = 0; iy < grid_nsteps; ++iy) {
        for (auto ix = 0; ix < grid_nsteps; ++ix) {
          auto x = ix - xcm;
          auto y = iy - ycm;
          auto w = TR[iy][ix] * std::pow(x*x + y*y, .5*n);
          // compute exp(i*n*phi) the naive way
          auto phi = std::atan2(y, x);
          real += w*std::cos(n*phi);
          imag += w*std::sin(n*phi);
          weight += w;
        }
      }
      auto ecc = std::sqrt(real*real + imag*imag) / weight;
      CHECK( ecc == Approx(event.eccentricity().at(n)).epsilon(1e-6).margin(1e-6) );
    }
  }

  // test grid size when step size does not evenly divide width
  auto var_map = make_var_map({
      {"normalization", 1.},
      {"reduced-thickness", 0.},
      {"grid-max", 10.},
      {"grid-step", 0.3}
  });

  Event event{var_map};

  CHECK( event.reduced_thickness_grid().num_dimensions() == 2 );
  CHECK( event.reduced_thickness_grid().shape()[0] == 67 );
  CHECK( event.reduced_thickness_grid().shape()[1] == 67 );
}
