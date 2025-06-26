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

#ifndef FLUIDCELLINFO_H
#define FLUIDCELLINFO_H

#include "RealType.h"

namespace Jetscape {

class FluidCellInfo {
 public:
  // data structure for outputing fluid cell information
  Jetscape::real energy_density = 0.;  //!< Local energy density [GeV/fm^3].
  Jetscape::real entropy_density = 0.; //!< Local entropy density [1/fm^3].
  Jetscape::real temperature = 0.;     //!< Local temperature [GeV].
  Jetscape::real pressure = 0.;        //!< Thermal pressure [GeV/fm^3].
  Jetscape::real
      qgp_fraction = 0.; //!< Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
  Jetscape::real mu_B = 0.;       //!< Net baryon chemical potential [GeV].
  Jetscape::real mu_C = 0.;       //!< Net charge chemical potential [GeV]
  Jetscape::real mu_S = 0.;       //!< Net strangeness chemical potential [GeV].
  Jetscape::real vx = 0., vy = 0., vz = 0.; //!< Flow velocity.
  Jetscape::real pi[4][4] = {{0., 0., 0., 0.},
                             {0., 0., 0., 0.},
                             {0., 0., 0., 0.},
                             {0., 0., 0., 0.}};   //!< Shear stress tensor [GeV/fm^3].
  Jetscape::real bulk_Pi = 0.;    //!< Bulk viscous pressure [GeV/fm^3].

  /** Default constructor.*/
  FluidCellInfo();

  /** @param b Multiply the fluid cell by scalar factor b. */
  FluidCellInfo inline operator*=(Jetscape::real b);

  /** Prints fluid cell properties to the screen. */
  // void Print();
};

// overload +-*/ for easier linear interpolation
/// adds \f$ c = a + b \f$
inline FluidCellInfo operator+(FluidCellInfo a, const FluidCellInfo &b) {
  a.energy_density += b.energy_density;
  a.entropy_density += b.entropy_density;
  a.temperature += b.temperature;
  a.pressure += b.pressure;
  a.qgp_fraction += b.qgp_fraction;
  a.mu_B += b.mu_B;
  a.mu_C += b.mu_C;
  a.mu_S += b.mu_S;
  a.vx += b.vx;
  a.vy += b.vy;
  a.vz += b.vz;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      a.pi[i][j] += b.pi[i][j];
    }
  }
  a.bulk_Pi += b.bulk_Pi;
  return a;
}

// Multiply the fluid cell with a scalar factor
FluidCellInfo inline FluidCellInfo::operator*=(Jetscape::real b) {
  this->energy_density *= b;
  this->entropy_density *= b;
  this->temperature *= b;
  this->pressure *= b;
  this->qgp_fraction *= b;
  this->mu_B *= b;
  this->mu_C *= b;
  this->mu_S *= b;
  this->vx *= b;
  this->vy *= b;
  this->vz *= b;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      this->pi[i][j] *= b;
    }
  }
  this->bulk_Pi *= b;
  return *this;
}

/// multiply \f$ c = a * b \f$
inline FluidCellInfo operator*(Jetscape::real a, FluidCellInfo b) {
  b *= a;
  return b;
}

/// multiply \f$ c = a * b \f$
inline FluidCellInfo operator*(FluidCellInfo a, Jetscape::real b) {
  a *= b;
  return a;
}

/// division \f$ c = a / b \f$
inline FluidCellInfo operator/(FluidCellInfo a, Jetscape::real b) {
  a *= 1.0 / b;
  return a;
}

// print the fluid cell information for debuging
// this function has bugs
// std::ostream &operator<<(std::ostream &os, const FluidCellInfo &cell) {
//    os << "energy_density=" << cell.energy_density << std::endl;
//    os << "entropy_density=" << cell.entropy_density << std::endl;
//    os << "temperature=" << cell.temperature << std::endl;
//    os << "pressure=" << cell.pressure << std::endl;
//    os << "qgp_fraction=" << cell.qgp_fraction << std::endl;
//    os << "mu_B=" << cell.mu_B << std::endl;
//    os << "mu_C=" << cell.mu_C << std::endl;
//    os << "mu_S=" << cell.mu_S << std::endl;
//    os << "vx=" << cell.vx << std::endl;
//    os << "vy=" << cell.vy << std::endl;
//    os << "vz=" << cell.vz << std::endl;
//    os << "pi[mu][nu]=" << std::endl;
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            os << cell.pi[i][j] << ' ';
//        }
//        os << std::endl;
//    }
//    os << "bulk_Pi=" << cell.bulk_Pi;
//    return os << std::endl;
//}

}  // end namespace Jetscape

#endif  // FluidCellInfo
