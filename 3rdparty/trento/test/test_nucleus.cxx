// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleus.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>

#include "catch.hpp"
#include "util.h"

#include "../src/hdf5_utils.h"
#include "../src/random.h"

using namespace trento;

TEST_CASE( "proton" ) {
  auto nucleus = Nucleus::create("p");

  CHECK( dynamic_cast<Proton*>(nucleus.get()) != nullptr );

  // proton contains one nucleon
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == 1 );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 1 );

  // and has zero radius
  CHECK( nucleus->radius() == 0. );

  // sample position with random offset
  double offset = random::canonical<>();
  nucleus->sample_nucleons(offset);
  auto&& proton = *(nucleus->begin());

  // check correct position
  CHECK( proton.x() == offset );
  CHECK( proton.y() == 0. );

  // not a participant initially
  CHECK( !proton.is_participant() );
}

TEST_CASE( "deuteron" ) {
  auto nucleus = Nucleus::create("d");

  CHECK( dynamic_cast<Deuteron*>(nucleus.get()) != nullptr );

  CHECK( std::distance(nucleus->begin(), nucleus->end()) == 2 );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 2 );

  nucleus->sample_nucleons(0);

  CHECK( nucleus->cbegin()->x() == -std::next(nucleus->cbegin())->x() );
  CHECK( nucleus->cbegin()->y() == -std::next(nucleus->cbegin())->y() );

  CHECK( nucleus->radius() == Approx(-std::log(.01)/(2*0.457)) );
}

TEST_CASE( "lead nucleus" ) {
  auto nucleus = Nucleus::create("Pb");

  CHECK( dynamic_cast<WoodsSaxonNucleus*>(nucleus.get()) != nullptr );

  constexpr int A = 208;
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == A );

  CHECK( nucleus->radius() > 6. );

  double offset = nucleus->radius() * random::canonical<>();
  nucleus->sample_nucleons(offset);

  // average nucleon position
  double x = 0., y = 0.;
  for (const auto& nucleon : *nucleus) {
    x += nucleon.x();
    y += nucleon.y();
  }
  x /= A;
  y /= A;
  auto tolerance = .6;
  CHECK( std::abs(x - offset) < tolerance );
  CHECK( std::abs(y) < tolerance );

  // no initial participants
  bool initial_participants = std::any_of(
      nucleus->cbegin(), nucleus->cend(),
      [](decltype(*nucleus->cbegin())& n) {
        return n.is_participant();
      });
  CHECK( !initial_participants );
}

TEST_CASE( "copper nucleus" ) {
  auto nucleus = Nucleus::create("Cu");
  auto def_nucleus = Nucleus::create("Cu2");

  CHECK( dynamic_cast<WoodsSaxonNucleus*>(nucleus.get()) != nullptr );
  CHECK( dynamic_cast<DeformedWoodsSaxonNucleus*>(def_nucleus.get()) != nullptr );

  constexpr int A = 62;
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(def_nucleus->cbegin(), def_nucleus->cend()) == A );

  CHECK( nucleus->radius() > 4. );
  CHECK( def_nucleus->radius() > nucleus->radius() );
}

TEST_CASE( "gold nucleus" ) {
  auto nucleus = Nucleus::create("Au");
  auto def_nucleus = Nucleus::create("Au2");

  CHECK( dynamic_cast<WoodsSaxonNucleus*>(nucleus.get()) != nullptr );
  CHECK( dynamic_cast<DeformedWoodsSaxonNucleus*>(def_nucleus.get()) != nullptr );

  constexpr int A = 197;
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(def_nucleus->cbegin(), def_nucleus->cend()) == A );

  CHECK( nucleus->radius() > 6. );
  CHECK( def_nucleus->radius() > nucleus->radius() );
}

TEST_CASE( "uranium nucleus" ) {
  auto nucleus = Nucleus::create("U");

  CHECK( dynamic_cast<DeformedWoodsSaxonNucleus*>(nucleus.get()) != nullptr );

  constexpr int A = 238;
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == A );

  CHECK( nucleus->radius() > 8. );
}

#ifdef TRENTO_HDF5

TEST_CASE( "manual nucleus" ) {
  std::vector<float> positions = {
     0.,  0., 0.,
     1.,  0., 0.,
    -1.,  0., 0.,
     3.,  0., 0.,
     0.,  2., 0.,
     0., -2., 0.,
  };

  const auto A = positions.size() / 3;

  temporary_path temp{".hdf5"};

  {
    auto dataspace = hdf5::make_dataspace(std::array<hsize_t, 3>{1, A, 3});
    const auto& datatype = hdf5::type<decltype(positions)::value_type>();

    H5::H5File file{temp.path.string(), H5F_ACC_EXCL};

    auto dataset = file.createDataSet("test", datatype, dataspace);
    dataset.write(positions.data(), datatype);
  }

  auto nucleus = Nucleus::create(temp.path.string());

  CHECK( dynamic_cast<ManualNucleus*>(nucleus.get()) != nullptr );

  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == A );

  CHECK( nucleus->radius() == Approx(3.) );

  nucleus->sample_nucleons(0.);

  auto nucleon = nucleus->cbegin();

  CHECK( nucleon->x() == Approx(0.) );
  CHECK( nucleon->y() == Approx(0.) );

  std::advance(nucleon, 1);
  CHECK( nucleon->x() == Approx(-std::next(nucleon)->x()) );
  CHECK( nucleon->y() == Approx(-std::next(nucleon)->y()) );

  CHECK( nucleon->x() == Approx(std::next(nucleon, 2)->x()/3) );
  CHECK( nucleon->y() == Approx(std::next(nucleon, 2)->y()/3) );

  std::advance(nucleon, 3);
  CHECK( nucleon->x() == Approx(-std::next(nucleon)->x()) );
  CHECK( nucleon->y() == Approx(-std::next(nucleon)->y()) );

  CHECK_THROWS_AS( Nucleus::create("nonexistent.hdf"), std::invalid_argument );
}

#endif  // TRENTO_HDF5

TEST_CASE( "woods-saxon sampling" ) {
  int A = 200;
  double R = 6., a = .5;
  NucleusPtr nucleus{new WoodsSaxonNucleus{static_cast<std::size_t>(A), R, a}};

  // Test Woods-Saxon sampling.
  // This is honestly not a great test; while it does prove that the code
  // basically works as intended, it does not rigorously show that the generated
  // numbers are actually Woods-Saxon distributed.  Personally I am much more
  // convinced by plotting a histogram and visually comparing it to the smooth
  // curve.  The script 'plot-woods-saxon.py' in the 'woods-saxon' subdirectory
  // does this.

  // sample a bunch of nuclei and bin all the nucleon positions
  // using "p" for cylindrical radius to distinguish from spherical "r"
  auto dp = .5;
  auto nevents = 1000;
  auto nsamples = nevents * A;
  std::map<int, int> hist{};
  for (auto i = 0; i < nevents; ++i) {
    nucleus->sample_nucleons(0.);
    for (const auto& nucleon : *nucleus) {
      auto x = nucleon.x();
      auto y = nucleon.y();
      auto p = std::sqrt(x*x + y*y);
      ++hist[p/dp];
    }
  }

  // integrate the W-S dist from pmin to pmax and over all z
  double rmax = R + 10.*a;
  auto integrate_woods_saxon = [R, a, rmax](double pmin, double pmax) -> double {
    int nzbins = 1000;
    auto npbins = static_cast<int>(nzbins*(pmax-pmin)/rmax);
    double result = 0.;
    for (auto ip = 0; ip <= npbins; ++ip) {
      auto p = pmin + (ip*(pmax-pmin))/(npbins - 1);
      for (auto iz = 0; iz < nzbins; ++iz) {
        auto z = (iz*rmax)/(nzbins-1);
        result += p/(1.+std::exp((std::sqrt(p*p + z*z) - R)/a));
      }
    }
    return result;
  };

  double ws_norm = integrate_woods_saxon(0, rmax);

  // check all histogram bins agree with numerical integration
  auto all_bins_correct = true;
  std::ostringstream bin_output{};
  bin_output << std::fixed << std::boolalpha
             << "pmin     pmax     prob     cprob    ratio    pass\n";
  for (const auto& bin : hist) {
    auto pmin = bin.first * dp;
    auto pmax = pmin + dp;
    auto prob = static_cast<double>(bin.second) / nsamples;
    auto correct_prob = integrate_woods_saxon(pmin, pmax) / ws_norm;
    bool within_tol = prob == Approx(correct_prob).epsilon(.1);
    if (!within_tol)
      all_bins_correct = false;
    bin_output << pmin << ' '
               << pmax << ' '
               << prob << ' '
               << correct_prob << ' '
               << prob/correct_prob << ' '
               << within_tol << '\n';
  }
  INFO( bin_output.str() );
  CHECK( all_bins_correct );
}

TEST_CASE( "deformed woods-saxon sampling" ) {
  // This is tough to test because the sampled positions are randomly rotated
  // and the z-coordinate is discarded.  However, the authors have visually
  // checked the results by histogramming the samples (not done here, but easy
  // to reproduce).

  // As a not-very-stringent test, check that the mean ellipticity of a deformed
  // W-S nucleus is significantly larger than a symmetric nucleus.

  std::size_t A = 200;
  double R = 6., a = .5, beta2 = .3, beta4 = .1;

  auto nucleus_sym = WoodsSaxonNucleus{A, R, a};
  auto nucleus_def = DeformedWoodsSaxonNucleus{A, R, a, beta2, beta4};

  auto mean_ecc2 = [](Nucleus& nucleus) {
    double e2 = 0.;
    int nevents = 500;

    for (int n = 0; n < nevents; ++n) {
      nucleus.sample_nucleons(0.);
      double numer = 0.;
      double denom = 0.;
      for (const auto& nucleon : nucleus) {
        double x2 = std::pow(nucleon.x(), 2);
        double y2 = std::pow(nucleon.y(), 2);
        numer += x2 - y2;
        denom += x2 + y2;
      }

      e2 += std::fabs(numer) / denom;
    }

    return e2 / nevents;
  };

  CHECK( mean_ecc2(nucleus_def) > 2.*mean_ecc2(nucleus_sym) );
}

TEST_CASE( "nuclear radius" ) {
  constexpr auto R = 5., a = .5;
  constexpr auto radius = R + 3*a;

  WoodsSaxonNucleus nucleus{100, R, a};
  DeformedWoodsSaxonNucleus dummy_def_nucleus{100, R, a, 0., 0.};

  CHECK( nucleus.radius() == Approx(radius) );
  CHECK( dummy_def_nucleus.radius() == Approx(radius) );

  DeformedWoodsSaxonNucleus def_nucleus{100, R, a, .2, .1};
  CHECK( def_nucleus.radius() > radius );
}

TEST_CASE( "unknown nucleus species" ) {
  CHECK_THROWS_AS( Nucleus::create("hello"), std::invalid_argument );
}
