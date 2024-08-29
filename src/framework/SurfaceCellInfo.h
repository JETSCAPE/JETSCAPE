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
// This is a general basic class for hydrodynamics

#ifndef SURFACECELLINFO_H
#define SURFACECELLINFO_H

#include "RealType.h"
#include <string>

namespace Jetscape {

class SurfaceCellInfo {
 public:
  // data structure for outputing hyper-surface information
  Jetscape::real tau;
  Jetscape::real x;
  Jetscape::real y;
  Jetscape::real eta;
  Jetscape::real d3sigma_mu[4];    //!< Surface vector.
  Jetscape::real energy_density;   //!< Local energy density [GeV/fm^3].
  Jetscape::real entropy_density;  //!< Local entropy density [1/fm^3].
  Jetscape::real temperature;      //!< Local temperature [GeV].
  Jetscape::real pressure;         //!< Thermal pressure [GeV/fm^3].
  Jetscape::real qgp_fraction;     //!< Fraction of quark gluon plasma assuming
                                   //!< medium is in QGP+HRG phase.
  Jetscape::real mu_B;             //!< Net baryon chemical potential [GeV].
  Jetscape::real mu_Q;             //!< Net charge chemical potential [GeV].
  Jetscape::real mu_S;     //!< Net strangeness chemical potential [GeV].
  Jetscape::real umu[4];   //!< Flow velocity.
  Jetscape::real pi[10];   //!< Shear stress tensor [GeV/fm^3].
  Jetscape::real bulk_Pi;  //!< Bulk viscous pressure [GeV/fm^3].

  /** Default constructor. */
  SurfaceCellInfo() = default;

  /** Destructor. */
  ~SurfaceCellInfo(){};

  /** Function to return member variables in a string */
  std::string sfi_to_string() const {
    std::string str = "tau = " + std::to_string(tau) + ", x = " +
                      std::to_string(x) + ", y = " + std::to_string(y) +
                      ", eta = " + std::to_string(eta) + "\n";
    str += "d3sigma_mu = " + std::to_string(d3sigma_mu[0]) + ", " +
           std::to_string(d3sigma_mu[1]) + ", " + std::to_string(d3sigma_mu[2]) +
           ", " + std::to_string(d3sigma_mu[3]) + "\n";
    str += "energy_density = " + std::to_string(energy_density) +
           ", entropy_density = " + std::to_string(entropy_density) +
           ", temperature = " + std::to_string(temperature) +
           ", pressure = " + std::to_string(pressure) +
           ", qgp_fraction = " + std::to_string(qgp_fraction) +
           ", mu_B = " + std::to_string(mu_B) + ", mu_Q = " + std::to_string(mu_Q) +
           ", mu_S = " + std::to_string(mu_S) + "\n";
    str += "umu = " + std::to_string(umu[0]) + ", " + std::to_string(umu[1]) +
           ", " + std::to_string(umu[2]) + ", " + std::to_string(umu[3]) + "\n";
    str += "pi = " + std::to_string(pi[0]) + ", " + std::to_string(pi[1]) + ", " +
           std::to_string(pi[2]) + ", " + std::to_string(pi[3]) + ", " +
           std::to_string(pi[4]) + ", " + std::to_string(pi[5]) + ", " +
           std::to_string(pi[6]) + ", " + std::to_string(pi[7]) + ", " +
           std::to_string(pi[8]) + ", " + std::to_string(pi[9]) + "\n";
    str += "bulk_Pi = " + std::to_string(bulk_Pi) + "\n";
    return str;
  }

  /** Function to compare two SurfaceCellInfo objects */
  bool operator==(const SurfaceCellInfo &rhs) const {
       return (tau == rhs.tau && x == rhs.x && y == rhs.y && eta == rhs.eta &&
              d3sigma_mu[0] == rhs.d3sigma_mu[0] &&
              d3sigma_mu[1] == rhs.d3sigma_mu[1] &&
              d3sigma_mu[2] == rhs.d3sigma_mu[2] &&
              d3sigma_mu[3] == rhs.d3sigma_mu[3] &&
              energy_density == rhs.energy_density &&
              entropy_density == rhs.entropy_density &&
              temperature == rhs.temperature && pressure == rhs.pressure &&
              qgp_fraction == rhs.qgp_fraction && mu_B == rhs.mu_B &&
              mu_Q == rhs.mu_Q && mu_S == rhs.mu_S && umu[0] == rhs.umu[0] &&
              umu[1] == rhs.umu[1] && umu[2] == rhs.umu[2] &&
              umu[3] == rhs.umu[3] && pi[0] == rhs.pi[0] && pi[1] == rhs.pi[1] &&
              pi[2] == rhs.pi[2] && pi[3] == rhs.pi[3] && pi[4] == rhs.pi[4] &&
              pi[5] == rhs.pi[5] && pi[6] == rhs.pi[6] && pi[7] == rhs.pi[7] &&
              pi[8] == rhs.pi[8] && pi[9] == rhs.pi[9] && bulk_Pi == rhs.bulk_Pi);
  }
};

}  // namespace Jetscape

#endif  // SURFACECELLINFO_H
