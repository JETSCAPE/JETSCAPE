/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "ThermPtnSampler.h"

#include <algorithm>
#include <omp.h>
#include <vector>

#include "gtest/gtest.h"

TEST(ThermPtnSampler, TEST_constructor) {
  unsigned int seed = 12345;
  double Tc = 0.155;
  ThermalPartonSampler sampler(seed, Tc);

  // expect no partons initially
  EXPECT_EQ(sampler.nTot(), 0);
  EXPECT_EQ(sampler.th_nL(), 0);
  EXPECT_EQ(sampler.th_nS(), 0);

  // set some brick parameters and check the expected outcome
  sampler.brick_length_width(8.0, 5.0);
  EXPECT_EQ(sampler.getL(), 20.0);    // 2*L+4
  EXPECT_EQ(sampler.getW(), 14.0);    // 2*W+4
  EXPECT_EQ(sampler.getTime(), 8.0);  // L

  sampler.brick_flow(0.2, 0.3, 0.4);
  EXPECT_EQ(sampler.getVx(), 0.2);
  EXPECT_EQ(sampler.getVy(), 0.3);
  EXPECT_EQ(sampler.getVz(), 0.4);
}

TEST(ThermPtnSampler, TEST_getClosestCachedTemp) {
  unsigned int seed = 12345;
  double Tc = 0.155;
  ThermalPartonSampler sampler(seed, Tc);

  // the for the function is in fm^{-1}, not in GeV
  // expected: 0,95×(0,155 GeV / (0,197327053 GeV*fm)) = 0,746223074
  // when choosing a value smaller than 0.746223074
  // expectation up to 1e-6
  EXPECT_NEAR(getClosestCachedTemp(sampler.getCacheCDFLight(), 0.7),
              0.746223074, 1e-6);
  EXPECT_NEAR(getClosestCachedTemp(sampler.getCacheCDFStrange(), 0.7),
              0.746223074, 1e-6);
}

TEST(ThermPtnSampler, TEST_BesselK2Function) {
  unsigned int seed = 12345;
  double Tc = 0.155;
  ThermalPartonSampler sampler(seed, Tc);

  EXPECT_NEAR(sampler.BesselK2function(2.0), 0.253759754566, 1e-6);
}

TEST(ThermPtnSampler, TEST_brick) {
  unsigned int seed = 42;
  double Tc = 0.16;
  ThermalPartonSampler sampler(seed, Tc);
  sampler.brick_length_width(8.0, 5.0);
  sampler.brick_flow(0.0, 0.0, 0.0);

  sampler.sample_brick();

  // check the number of partons
  EXPECT_EQ(sampler.nTot(), 3025);
  EXPECT_EQ(sampler.th_nL(), 2416);
  EXPECT_EQ(sampler.th_nS(), 609);
  EXPECT_EQ(sampler.nTot(), sampler.th_nL() + sampler.th_nS());
}

TEST(ThermPtnSampler, TEST_sample_2p1d) {
  unsigned int seed = 42;
  // in this test we artificially set Tc to 0.8 to increase the number of
  // partons otherwise in this cell there is nothing sampled
  double Tc = 0.8;  // used for the T cache in this case

  std::vector<double> surface_cell = {1.70523,   -2.48733,    -2.23816,    0,
                                      0.0145614, -0.00055864, -0.00199173, 0,
                                      0.8,       -0.205297,   -0.310245,   0};
  std::vector<std::vector<double>> surface = {surface_cell};

  double eta_max_boost_inv = 2.0;

  std::vector<int> nTotResults;
  std::vector<int> th_nLResults;
  std::vector<int> th_nSResults;
  std::vector<int> nTotExpected;
  int highest_power = 5;

  // Redirect stdout and stderr to /dev/null and check for success
  if (std::freopen("/dev/null", "w", stdout) == nullptr) {
    std::cerr << "Failed to redirect stdout to /dev/null\n";
  }
  if (std::freopen("/dev/null", "w", stderr) == nullptr) {
    std::cerr << "Failed to redirect stderr to /dev/null\n";
  }

  for (int i = 0; i <= highest_power; i++) {
    ThermalPartonSampler sampler(seed, Tc);
    sampler.set_hypersurface(surface);
    int num_threads = pow(2, i);
    omp_set_num_threads(num_threads);
    sampler.sample_2p1d(eta_max_boost_inv);

    nTotResults.push_back(sampler.nTot());
    th_nLResults.push_back(sampler.th_nL());
    th_nSResults.push_back(sampler.th_nS());
    nTotExpected.push_back(sampler.th_nL() + sampler.th_nS());
  }

  // Restore stdout and stderr
  if (std::freopen("/dev/tty", "w", stdout) == nullptr) {
    std::cerr << "Failed to restore stdout from /dev/tty\n";
  }
  if (std::freopen("/dev/tty", "w", stderr) == nullptr) {
    std::cerr << "Failed to restore stderr from /dev/tty\n";
  }

  auto check_all_equal = [](std::vector<int> vec) {
    return std::adjacent_find(vec.begin(), vec.end(),
                              std::not_equal_to<int>()) == vec.end();
  };

  // check that threads didn't introduce inconsistencies
  EXPECT_TRUE(check_all_equal(nTotResults));
  EXPECT_TRUE(check_all_equal(th_nLResults));
  EXPECT_TRUE(check_all_equal(th_nSResults));
  EXPECT_TRUE(check_all_equal(nTotExpected));

  // check the number of partons
  EXPECT_EQ(nTotResults[0], 76);
  EXPECT_EQ(th_nLResults[0], 54);
  EXPECT_EQ(th_nSResults[0], 22);
  EXPECT_EQ(nTotExpected[0], th_nLResults[0] + th_nSResults[0]);
}

TEST(ThermPtnSampler, TEST_sample_3p1d) {
  unsigned int seed = 42;
  // in this test we artificially set Tc to 0.8 to increase the number of
  // partons otherwise in this cell there is nothing sampled
  double Tc = 0.8;  // used for the T cache in this case

  std::vector<double> surface_cell = {1.70523,   -2.48733,    -2.23816,    0,
                                      0.0145614, -0.00055864, -0.00199173, 0,
                                      0.8,       -0.205297,   -0.310245,   0};
  std::vector<std::vector<double>> surface = {surface_cell};

  std::vector<int> nTotResults;
  std::vector<int> th_nLResults;
  std::vector<int> th_nSResults;
  std::vector<int> nTotExpected;
  int highest_power = 5;
  bool cartesian = false;

  for (int i = 0; i <= highest_power; i++) {
    ThermalPartonSampler sampler(seed, Tc);
    sampler.set_hypersurface(surface);
    sampler.sample_3p1d(cartesian);

    nTotResults.push_back(sampler.nTot());
    th_nLResults.push_back(sampler.th_nL());
    th_nSResults.push_back(sampler.th_nS());
    nTotExpected.push_back(sampler.th_nL() + sampler.th_nS());
  }

  auto check_all_equal = [](std::vector<int> vec) {
    return std::adjacent_find(vec.begin(), vec.end(),
                              std::not_equal_to<int>()) == vec.end();
  };

  // check that threads didn't introduce inconsistencies
  EXPECT_TRUE(check_all_equal(nTotResults));
  EXPECT_TRUE(check_all_equal(th_nLResults));
  EXPECT_TRUE(check_all_equal(th_nSResults));
  EXPECT_TRUE(check_all_equal(nTotExpected));

  // check the number of partons
  EXPECT_EQ(nTotResults[0], 3);
  EXPECT_EQ(th_nLResults[0], 2);
  EXPECT_EQ(th_nSResults[0], 1);
  EXPECT_EQ(nTotExpected[0], th_nLResults[0] + th_nSResults[0]);
}