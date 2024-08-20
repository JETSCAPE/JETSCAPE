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
// -----------------------------------------
// This is a causal liquefier with the JETSCAPE framework
// -----------------------------------------

#include "CausalLiquefier.h"

#include <cfloat>

#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {

CausalLiquefier::CausalLiquefier() {
  VERBOSE(8);
  dtau = 0.6;
  dx = 0.3;
  dy = 0.3;
  deta = 0.3;
  tau_delay = 2.0;
  time_relax = 0.1;
  d_diff = 0.08;
  width_delta = 0.1;
  Init();  // Get values of parameters from XML
  c_diff = sqrt(d_diff / time_relax);
  gamma_relax = 0.5 / time_relax;
  if (c_diff > 1.0) {
    JSWARN << "Bad Signal Velocity in CausalLiquefier";
  }
  //    else{
  //        //for debug
  //        JSINFO << "c_diff = " << c_diff;
  //    }
}

CausalLiquefier::CausalLiquefier(double dtau_in, double dx_in, double dy_in,
                                 double deta_in) {
  VERBOSE(8);
  JSINFO << "Initialize CausalLiquefier (for Unit Test) ...";
  dtau = dtau_in;
  dx = dx_in;
  dy = dy_in;
  deta = deta_in;
  tau_delay = 2.0;
  time_relax = 0.1;
  d_diff = 0.08;
  width_delta = 0.1;

  c_diff = sqrt(d_diff / time_relax);
  gamma_relax = 0.5 / time_relax;

  JSINFO << "<CausalLiquefier> Fluid Time Step and Cell Size: dtau=" << dtau
         << " fm, dx=" << dx << " fm, dy=" << dy << " fm, deta=" << deta;
  JSINFO << "<CausalLiquefier> Parameters: tau_delay=" << tau_delay
         << " fm, time_relax=" << time_relax << " fm, d_diff=" << d_diff
         << " fm, width_delta=" << width_delta << " fm";
}

void CausalLiquefier::Init() {
  // Initialize parameter with values in XML
  JSINFO << "Initialize CausalLiquefier ...";

  dtau = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "dtau"});
  dx = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "dx"});
  dy = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "dy"});
  deta = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "deta"});
  tau_delay = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "tau_delay"});  // in [fm]
  time_relax = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "time_relax"});  // in [fm]
  d_diff = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "d_diff"});  // in [fm]
  width_delta = JetScapeXML::Instance()->GetElementDouble(
      {"Liquefier", "CausalLiquefier", "width_delta"});  // in [fm]

  // for debug
  //    JSINFO
  //    << "<CausalLiquefier> Fluid Time Step and Cell Size: dtau="
  //    << dtau << " fm, dx="
  //    << dx << " fm, dy="
  //    << dy << " fm, deta="
  //    << deta;
  //    JSINFO
  //    << "<CausalLiquefier> Parameters: tau_delay="
  //    << tau_delay << " fm, time_relax="
  //    << time_relax << " fm, d_diff="
  //    << d_diff << " fm, width_delta="
  //    << width_delta <<" fm";
}

// CausalLiquefier::~CausalLiquefier(){};

void CausalLiquefier::smearing_kernel(
    Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real eta,
    const Droplet drop_i, std::array<Jetscape::real, 4> &jmu) const {
  jmu = {0., 0, 0, 0};  // source in tau-eta coordinates

  const auto p_drop = drop_i.get_pmu();
  auto x_drop =
      drop_i.get_xmu();  // position of the doloplet in the Cartesian coodinates
  double tau_drop = x_drop[0];
  double eta_drop = x_drop[3];
  x_drop[0] = get_t(tau_drop, eta_drop);
  x_drop[3] = get_z(tau_drop, eta_drop);

  if (tau - 0.5 * dtau <= tau_drop + tau_delay &&
      tau + 0.5 * dtau > tau_drop + tau_delay) {
    double t =
        get_t(tau, eta);  // position in the fluid in the Cartesian coodinates
    double z = get_z(tau, eta);

    double delta_t = t - x_drop[0];
    double delta_r = sqrt((x - x_drop[1]) * (x - x_drop[1]) +
                          (y - x_drop[2]) * (y - x_drop[2]) +
                          (z - x_drop[3]) * (z - x_drop[3]));

    // get solutions of the diffusion equation in the Cartesian coordinates
    double jt = kernel_rho(delta_t, delta_r) / dtau;
    double jz;
    if (delta_r <= DBL_MIN) {
      jz = 0.0;
    } else {
      jz = ((z - x_drop[3]) / delta_r) * kernel_j(delta_t, delta_r) / dtau;
    }
    // get flux for the constant-tau surface
    double jtau = get_ptau(jt, jz, eta);

    // get source in tau-eta coordinates
    jmu[0] = jtau * get_ptau(p_drop[0], p_drop[3], eta);
    jmu[1] = jtau * p_drop[1];
    jmu[2] = jtau * p_drop[2];
    jmu[3] = jtau * get_peta(p_drop[0], p_drop[3], eta);
  }
}

// Charge density rho in causal diffusion
double CausalLiquefier::kernel_rho(double t, double r) const {
  return dumping(t) * (rho_smooth(t, r) + rho_delta(t, r));
}

// Radial component of current j in causal diffusion
double CausalLiquefier::kernel_j(double t, double r) const {
  return dumping(t) * (j_smooth(t, r) + j_delta(t, r));
}

// Dumping factor in solutions of rho and j
double CausalLiquefier::dumping(double t) const {
  return exp(-gamma_relax * t) / (4.0 * M_PI);
}

// Smooth component of rho
double CausalLiquefier::rho_smooth(double t, double r) const {
  if (r < c_diff * t) {
    double u = sqrt(c_diff * c_diff * t * t - r * r);
    double x = gamma_relax * u / c_diff;  // unitless

    double i1 = gsl_sf_bessel_I1(x) / (c_diff * u);
    double i2 = gsl_sf_bessel_In(2, x) * t / u / u;
    double f = gamma_relax * gamma_relax / c_diff;

    return f * (i1 + i2);
  } else {
    return 0.0;
  }
}

// Smooth component of j
double CausalLiquefier::j_smooth(double t, double r) const {
  if (r < c_diff * t) {
    double u = sqrt(c_diff * c_diff * t * t - r * r);
    double x = gamma_relax * u / c_diff;  // unitless

    double i2 = gsl_sf_bessel_In(2, x) * r / u / u;
    double f = gamma_relax * gamma_relax / c_diff;

    return f * i2;

  } else {
    return 0.0;
  }
}

// Wave front component of rho
double CausalLiquefier::rho_delta(double t, double r) const {
  double r_w = width_delta;
  if (c_diff * t <= width_delta) {
    r_w = c_diff * t;
  }

  if (r >= c_diff * t - r_w && r < c_diff * t) {
    double x = gamma_relax * t;
    return (1.0 + x + x * x / 2.0) / r_w / r / r;
  } else {
    return 0.0;
  }
}

// Wave front component of j
double CausalLiquefier::j_delta(double t, double r) const {
  return c_diff * rho_delta(t, r);
}

// Get Cartesian time t from tau and eta
double CausalLiquefier::get_t(double tau, double eta) const {
  return tau * cosh(eta);
}

// Get Cartesian coordinate z from tau and eta
double CausalLiquefier::get_z(double tau, double eta) const {
  return tau * sinh(eta);
}

// Lorentz Transformation to get tau component of four vector
double CausalLiquefier::get_ptau(double p0, double p3, double eta) const {
  return p0 * cosh(eta) - p3 * sinh(eta);
}

// Lorentz Transformation to get eta component of four vector
double CausalLiquefier::get_peta(double p0, double p3, double eta) const {
  return p3 * cosh(eta) - p0 * sinh(eta);
}

// For debug, Change tau_delay
void CausalLiquefier::set_t_delay(double new_tau_delay) {
  tau_delay = new_tau_delay;
}

};  // namespace Jetscape
