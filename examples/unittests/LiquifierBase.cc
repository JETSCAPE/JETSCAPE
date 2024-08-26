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

#include <vector>

#include "LiquefierBase.h"
#include "gtest/gtest.h"

using namespace Jetscape;

// check function to add a droplet and Clear
TEST(LiquefierBaseTest, TEST_add_a_droplet_and_Clear) {
  LiquefierBase lqf;
  EXPECT_EQ(0, lqf.get_dropletlist_size());

  Droplet a_drop;
  lqf.add_a_droplet(a_drop);
  EXPECT_EQ(1, lqf.get_dropletlist_size());

  lqf.Clear();
  EXPECT_EQ(0, lqf.get_dropletlist_size());
}

TEST(LiquefierBaseTest, TEST_filter_partons) {
  LiquefierBase lqf;
  FourVector p_test1(1.0, 0.0, 0.0, 1.0);    // px, py, pz, E
  FourVector x_test1(0.0, 0.0, 0.0, 1.0);    // x, y, z, t
  Parton test1(0, 21, 0, p_test1, x_test1);  // 21 = gluon

  FourVector p_test2(3.0, 0.0, 0.0, 3.0);
  FourVector x_test2(0.0, 0.0, 0.0, 1.0);
  Parton test2(0, 21, 0, p_test2, x_test2);

  std::vector<Parton> pOut;

  pOut.push_back(test1);
  pOut.push_back(test2);
  lqf.filter_partons(pOut);
  EXPECT_EQ(2, pOut.size());
  EXPECT_EQ(lqf.get_drop_stat(), pOut[0].pstat());
  EXPECT_EQ(0, pOut[1].pstat());
  pOut.clear();

  // test a negative particle
  FourVector p_test3(3.0, 0.0, 0.0, 3.0);
  FourVector x_test3(0.0, 0.0, 0.0, 1.0);
  Parton test3(0, 21, -1, p_test3, x_test3);
  pOut.push_back(test3);
  lqf.filter_partons(pOut);
  EXPECT_EQ(1, pOut.size());
  for (const auto &iparton : pOut)
    EXPECT_EQ(lqf.get_drop_stat(), iparton.pstat());
  pOut.clear();

  // test a negative particle and a soft particle
  pOut.push_back(test1);
  pOut.push_back(test3);
  lqf.filter_partons(pOut);
  EXPECT_EQ(2, pOut.size());
  for (const auto &iparton : pOut)
    EXPECT_EQ(lqf.get_drop_stat(), iparton.pstat());
  pOut.clear();

  // test a negative particle and a hard particle
  pOut.push_back(test2);
  pOut.push_back(test3);
  lqf.filter_partons(pOut);
  EXPECT_EQ(2, pOut.size());
  EXPECT_EQ(0, pOut[0].pstat());
  EXPECT_EQ(lqf.get_drop_stat(), pOut[1].pstat());
  pOut.clear();
}

TEST(LiquefierBaseTest, TEST_check_energy_momentum_conservation) {
  LiquefierBase lqf;
  FourVector p_test1(1.0, 0.0, 0.0, 1.0);    // px, py, pz, E
  FourVector x_test1(0.0, 0.0, 0.0, 1.0);    // x, y, z, t
  Parton test1(0, 21, 0, p_test1, x_test1);  // 21 = gluon

  FourVector p_test2(3.0, 0.0, 0.0, 3.0);
  FourVector x_test2(0.0, 0.0, 0.0, 1.0);
  Parton test2(0, 21, 0, p_test2, x_test2);

  std::vector<Parton> pIn;
  std::vector<Parton> pOut;

  pIn.push_back(test1);
  pOut.push_back(test2);
  lqf.check_energy_momentum_conservation(pIn, pOut);
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(2, pOut.size());
  EXPECT_DOUBLE_EQ(pIn[0].p(0), pOut[0].p(0) + pOut[1].p(0));
  EXPECT_DOUBLE_EQ(pIn[0].p(1), pOut[0].p(1) + pOut[1].p(1));
  EXPECT_DOUBLE_EQ(pIn[0].p(2), pOut[0].p(2) + pOut[1].p(2));
  EXPECT_DOUBLE_EQ(pIn[0].p(3), pOut[0].p(3) + pOut[1].p(3));
  pIn.clear();
  pOut.clear();
}

// check function to add a hydro source term
TEST(LiquefierBaseTest, TEST_add_hydro_source) {
  LiquefierBase lqf;

  FourVector p_test1(1.0, 0.0, 0.0, 1.0);    // px, py, pz, E
  FourVector x_test1(0.0, 0.0, 0.0, 1.0);    // x, y, z, t
  Parton test1(0, 21, 0, p_test1, x_test1);  // 21 = gluon

  FourVector p_test2(3.0, 0.0, 0.0, 3.0);
  FourVector x_test2(0.0, 0.0, 0.0, 1.0);
  Parton test2(0, 21, 0, p_test2, x_test2);

  std::vector<Parton> pIn;
  std::vector<Parton> pOut;
  pIn.push_back(test1);
  pIn.push_back(test2);

  // test when pOut.size() == 0
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(0, lqf.get_dropletlist_size());
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test when pOut == pIn (hard parton)
  pIn.push_back(test2);
  pOut.push_back(test2);
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(0, lqf.get_dropletlist_size());
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(1, pOut.size());
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test when pOut == pIn (soft parton)
  pIn.push_back(test1);
  pOut.push_back(test1);
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(1, lqf.get_dropletlist_size());
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(1, pOut.size());
  EXPECT_DOUBLE_EQ(p_test1.t(), (lqf.get_a_droplet(0).get_pmu())[0]);
  EXPECT_DOUBLE_EQ(p_test1.x(), (lqf.get_a_droplet(0).get_pmu())[1]);
  EXPECT_DOUBLE_EQ(x_test1.t(), (lqf.get_a_droplet(0).get_xmu())[0]);
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test a negative particle
  FourVector p_test3(3.0, 0.0, 0.0, 3.0);
  FourVector x_test3(0.0, 0.0, 0.0, 1.0);
  Parton test3(0, 21, -1, p_test3, x_test3);
  pIn.push_back(test2);
  pOut.push_back(test3);
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(1, lqf.get_dropletlist_size());
  EXPECT_DOUBLE_EQ(p_test3.t(), (lqf.get_a_droplet(0).get_pmu())[0]);
  EXPECT_DOUBLE_EQ(p_test3.x(), (lqf.get_a_droplet(0).get_pmu())[1]);
  EXPECT_DOUBLE_EQ(x_test3.t(), (lqf.get_a_droplet(0).get_xmu())[0]);
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(2, pOut.size());
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test a negative particle and a soft particle
  pIn.push_back(test2);
  pOut.push_back(test1);
  pOut.push_back(test3);
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(1, lqf.get_dropletlist_size());
  EXPECT_DOUBLE_EQ(p_test2.t(), (lqf.get_a_droplet(0).get_pmu())[0]);
  EXPECT_DOUBLE_EQ(p_test2.x(), (lqf.get_a_droplet(0).get_pmu())[1]);
  EXPECT_DOUBLE_EQ(x_test2.t(), (lqf.get_a_droplet(0).get_xmu())[0]);
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(3, pOut.size());
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test a negative particle and a hard particle
  pIn.push_back(test1);
  pOut.push_back(test2);
  pOut.push_back(test3);
  lqf.add_hydro_sources(pIn, pOut);
  EXPECT_EQ(1, lqf.get_dropletlist_size());
  EXPECT_DOUBLE_EQ(p_test1.t() - p_test2.t(),
                   (lqf.get_a_droplet(0).get_pmu())[0]);
  EXPECT_DOUBLE_EQ(p_test1.x() - p_test2.x(),
                   (lqf.get_a_droplet(0).get_pmu())[1]);
  EXPECT_DOUBLE_EQ((x_test1.t() + x_test2.t()) / 2.,
                   (lqf.get_a_droplet(0).get_xmu())[0]);
  EXPECT_EQ(1, pIn.size());
  EXPECT_EQ(3, pOut.size());
  pIn.clear();
  pOut.clear();
  lqf.Clear();

  // test a soft particle in pIn only
  // pIn.push_back(test1);
  // lqf.add_hydro_sources(pIn, pOut);
  // EXPECT_EQ(1, lqf.get_dropletlist_size());
  // EXPECT_DOUBLE_EQ(p_test1.t(), (lqf.get_a_droplet(0).get_pmu())[0]);
  // EXPECT_DOUBLE_EQ(p_test1.x(), (lqf.get_a_droplet(0).get_pmu())[1]);
  // EXPECT_DOUBLE_EQ(x_test1.t(), (lqf.get_a_droplet(0).get_xmu())[0]);
  // EXPECT_EQ(0, pIn.size());
  // EXPECT_EQ(0, pOut.size());
  // pIn.clear();
  // pOut.clear();
  // lqf.Clear();
}
