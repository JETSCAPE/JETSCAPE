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

#include <iostream>

#include "CausalLiquefier.h"
#include "LiquefierBase.h"
#include "gtest/gtest.h"

using namespace Jetscape;

// check coordinate transformation
TEST(CausalLiquifierTest, TEST_COORDINATES) {
  CausalLiquefier lqf(0.3, 0.3, 0.3, 0.3);

  // for the transformation from tau-eta to t-z (configuration space)
  EXPECT_DOUBLE_EQ(0.0, lqf.get_t(0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_t(5.0, 0.0));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_z(0.0, 0.5));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_z(5.0, 0.0));

  for (double tau = 0.1; tau < 0.3; tau += 0.1) {
    for (double eta = 0.0; eta < 0.2; eta += 0.1) {
      double t = lqf.get_t(tau, eta);
      double z = lqf.get_z(tau, eta);
      EXPECT_DOUBLE_EQ(tau * tau, t * t - z * z);
      EXPECT_DOUBLE_EQ(exp(2.0 * eta), (t + z) / (t - z));
    }
  }

  // for the transformation from t-z to tau-eta (momentum)
  EXPECT_DOUBLE_EQ(0.0, lqf.get_ptau(0.0, 0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_ptau(5.0, 0.0, 0.0));
  EXPECT_DOUBLE_EQ(0.0, lqf.get_peta(0.0, 0.0, 0.5));
  EXPECT_DOUBLE_EQ(5.0, lqf.get_peta(0.0, 5.0, 0.0));
}

// check causality
TEST(CausalLiquifierTest, TEST_CAUSALITY) {
  CausalLiquefier lqf(0.3, 0.3, 0.3, 0.3);

  // check zero source outside of the cousal area
  double time = lqf.tau_delay;
  double r_bound = lqf.c_diff * time;
  for (int scale = 1.0; scale < 5.0; scale++) {
    double r_out = r_bound * scale;
    EXPECT_DOUBLE_EQ(0.0, lqf.rho_smooth(time, r_out));
    EXPECT_DOUBLE_EQ(0.0, lqf.rho_delta(time, r_out));
    EXPECT_DOUBLE_EQ(0.0, lqf.kernel_rho(time, r_out));
  }

  // check zero source before the deposition time
  time = 0.99 * (lqf.tau_delay - 0.5 * lqf.dtau);
  std::array<Jetscape::real, 4> jmu = {0.0, 0.0, 0.0, 0.0};
  std::array<Jetscape::real, 4> droplet_xmu = {0.0, 0.0, 0.0, 0.0};
  std::array<Jetscape::real, 4> droplet_pmu = {1.0, 1.0, 0.0, 0.0};
  Droplet drop_i(droplet_xmu, droplet_pmu);
  lqf.smearing_kernel(time, 0.0, 0.0, 0.0, drop_i, jmu);
  EXPECT_DOUBLE_EQ(0.0, jmu[0]);

  // check zero source after the deposition time
  time = 1.01 * (lqf.tau_delay + 0.5 * lqf.dtau);
  lqf.smearing_kernel(time, 0.0, 0.0, 0.0, drop_i, jmu);
  EXPECT_DOUBLE_EQ(0.0, jmu[0]);
}

// check conservation law
TEST(CausalLiquifierTest, TEST_CONSERVATION) {
  double dr = 0.005;  // in [fm]
  double dt = 0.005;  // in [fm]

  CausalLiquefier lqf(dt, dr, dr, dr);

  for (double t = 0.3; t < 3.0; t += dt) {
    double integrated_value = 0.0;
    for (double r = 0.5 * dr; r < 1.5 * t; r += dr) {
      integrated_value += 4.0 * M_PI * r * r * dr * lqf.kernel_rho(t, r);
    }
    // require conservation within 5%
    EXPECT_NEAR(1.0, integrated_value, 0.05);
  }
}

// check conservation law
TEST(CausalLiquifierTest, TEST_GRID_CARTESIAN_CONSERVATION) {
  double ll[3] = {0.05, 0.1, 0.3};

  for (int i = 0; i < 0; i++) {
    double dt = 0.1;
    double l = ll[i];
    double dx = l;
    double dy = l;
    double dz = l;
    int n_cells = 5.0 / l;
    double dV = dx * dx * dz;

    CausalLiquefier lqf(dt, dx, dy, dz);

    // std::o/Users/yasukitachibana/Dropbox
    // (Personal)/Codes/Release2.2CandidateLiquefierUpdate/examples/unittests/causal_liquifier.ccfstream
    // ofs; string filename = "liq_cons_test_dx" + std::to_string(int(l*1000)) +
    // "_tr3000_d2400_delta300.txt"; ofs.open(filename.c_str(),
    // std::ios_base::out);

    for (double t = 0.5; t < 3.0; t += dt) {
      double integrated_value = 0.0;
      for (int ix = 0; ix < n_cells; ix++) {
        double x = (ix - 0.5 * double(n_cells)) * dx;
        for (int iy = 0; iy < n_cells; iy++) {
          double y = (iy - 0.5 * double(n_cells)) * dy;
          for (int iz = 0; iz < n_cells; iz++) {
            double z = (iz - 0.5 * double(n_cells)) * dz;

            double r = sqrt(x * x + y * y + z * z);
            integrated_value += lqf.kernel_rho(t, r) * dV;
          }
        }
      }
      // require conservation within 5%
      // ofs << t << " " << integrated_value <<"\n";
      EXPECT_NEAR(1.0, integrated_value, 0.05);
      // JSINFO <<"t="<<t <<" " <<integrated_value;
    }
    // ofs.close();
  }
}

// check conservation law
TEST(CausalLiquifierTest, TEST_GRID_TAU_ETA_CONSERVATION) {
  std::array<Jetscape::real, 4> x_in = {1.0, 0.0, 0.0, 1.0};  // in tau-eta
  std::array<Jetscape::real, 4> p_in = {1.0, 1.0, 1.0, 1.0};  // in Cartesian
  Droplet a_drop(x_in, p_in);

  double dtau = 0.1;
  double dl = 0.05;
  double dx = dl;
  double dy = dl;
  double deta = 0.05;
  int n_xy = 8.0 / dl;
  int n_eta = 5.0 / deta;

  CausalLiquefier lqf(dtau, dx, dy, deta);

  //    std::ofstream ofs;
  //    string filename =
  //    "liq_cons_tau_eta_test_tdep"+std::to_string(int(x_in[0]))+"_etadep"+std::to_string(int(x_in[3]))+"_dxy"
  //    + std::to_string(int(dl*1000)) +
  //    + "_deta" + std::to_string(int(deta*1000)) + "_tr100_d80_delta100.txt";
  //
  //    ofs.open(filename.c_str(), std::ios_base::out);

  double dt = 0.5;
  for (double t = 0.2; t < 4.0; t += dt) {
    double integrated_value = 0.0;
    double tau_delay = t;
    lqf.set_t_delay(tau_delay);
    // in tau-eta
    std::array<Jetscape::real, 4> total_pmu = {0.0, 0.0, 0.0, 0.0};
    std::array<Jetscape::real, 4> x_hydro = {0.0, 0.0, 0.0, 0.0};
    x_hydro[0] = x_in[0] + tau_delay;
    double dvolume = x_hydro[0] * dx * dy * deta;

    for (int ix = 0; ix < n_xy; ix++) {
      x_hydro[1] = (ix - 0.5 * double(n_xy)) * dx;
      for (int iy = 0; iy < n_xy; iy++) {
        x_hydro[2] = (iy - 0.5 * double(n_xy)) * dy;
        for (int ieta = 0; ieta < n_eta; ieta++) {
          x_hydro[3] = (ieta - 0.5 * double(n_eta)) * deta;

          std::array<Jetscape::real, 4> jmu = {0.0, 0.0, 0.0, 0.0};

          lqf.smearing_kernel(x_hydro[0], x_hydro[1], x_hydro[2], x_hydro[3],
                              a_drop, jmu);

          total_pmu[0] +=
              dtau * (jmu[0] * cosh(x_hydro[3]) + jmu[3] * sinh(x_hydro[3])) *
              dvolume;
          total_pmu[1] += dtau * jmu[1] * dvolume;
          total_pmu[2] += dtau * jmu[2] * dvolume;
          total_pmu[3] +=
              dtau * (jmu[0] * sinh(x_hydro[3]) + jmu[3] * cosh(x_hydro[3])) *
              dvolume;

          integrated_value += dtau * jmu[1] * dvolume;
        }
      }
    }
    //        ofs << t << " " << integrated_value <<"\n";
    //        JSINFO
    //        << total_pmu[0] << " "
    //        << total_pmu[1] << " "
    //        << total_pmu[2] << " "
    //        << total_pmu[3];
    //        ofs.close();
    if (t > 3.) {
      EXPECT_NEAR(1.0, integrated_value, 0.05);
    }
  }
}
