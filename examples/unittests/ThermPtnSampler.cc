/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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
#include "gtest/gtest.h"

#include <vector>

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
    EXPECT_EQ(sampler.getL(), 20.0); // 2*L+4
    EXPECT_EQ(sampler.getW(), 14.0); // 2*W+4
    EXPECT_EQ(sampler.getTime(), 8.0); // L

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
    EXPECT_NEAR(getClosestCachedTemp(sampler.getCacheCDFLight() , 0.7), 0.746223074, 1e-6); 
    EXPECT_NEAR(getClosestCachedTemp(sampler.getCacheCDFStrange(), 0.7), 0.746223074, 1e-6);
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
    EXPECT_EQ(sampler.nTot(), 3079);
    EXPECT_EQ(sampler.th_nL(), 2444);
    EXPECT_EQ(sampler.th_nS(), 635);
    EXPECT_EQ(sampler.nTot(), sampler.th_nL() + sampler.th_nS());
}

TEST(ThermPtnSampler, TEST_sample_2p1d) {
    unsigned int seed = 42;
    // in this test we artificially set Tc to 0.8 to increase the number of partons
    // otherwise in this cell there is nothing sampled
    double Tc = 0.8; // used for the T cache in this case

    ThermalPartonSampler sampler(seed, Tc);

    std::vector<double> surface_cell = {1.70523, -2.48733, -2.23816, 0, 0.0145614, -0.00055864, -0.00199173, 0, 0.8, -0.205297, -0.310245, 0};
    std::vector<std::vector<double>> surface = {surface_cell};
    sampler.set_hypersurface(surface);

    double eta_max_boost_inv = 2.0;
    sampler.sample_2p1d(eta_max_boost_inv);

    // check the number of partons
    EXPECT_EQ(sampler.nTot(), 70);
    EXPECT_EQ(sampler.th_nL(), 47);
    EXPECT_EQ(sampler.th_nS(), 23);
    EXPECT_EQ(sampler.nTot(), sampler.th_nL() + sampler.th_nS());
}

TEST(ThermPtnSampler, TEST_sample_3p1d) {
    unsigned int seed = 42;
    // in this test we artificially set Tc to 0.8 to increase the number of partons
    // otherwise in this cell there is nothing sampled
    double Tc = 0.8; // used for the T cache in this case

    ThermalPartonSampler sampler(seed, Tc);

    std::vector<double> surface_cell = {1.70523, -2.48733, -2.23816, 0, 0.0145614, -0.00055864, -0.00199173, 0, 0.8, -0.205297, -0.310245, 0};
    std::vector<std::vector<double>> surface = {surface_cell};
    sampler.set_hypersurface(surface);

    bool cartesian = false;
    sampler.sample_3p1d(cartesian);

    // check the number of partons
    EXPECT_EQ(sampler.nTot(), 4);
    EXPECT_EQ(sampler.th_nL(), 3);
    EXPECT_EQ(sampler.th_nS(), 1);
    EXPECT_EQ(sampler.nTot(), sampler.th_nL() + sampler.th_nS());
}