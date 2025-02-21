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

/**
 * @class FluidCellInfo
 * @brief Data structure for outputting fluid cell information.
 *
 * This class contains various properties of a fluid cell, including energy density,
 * entropy density, temperature, pressure, and chemical potentials, among others.
 * It is used to store and manipulate information about the state of a fluid cell
 * in a simulation.
 *
 * @var Jetscape::real FluidCellInfo::energy_density
 * Local energy density [GeV/fm^3].
 *
 * @var Jetscape::real FluidCellInfo::entropy_density
 * Local entropy density [1/fm^3].
 *
 * @var Jetscape::real FluidCellInfo::temperature
 * Local temperature [GeV].
 *
 * @var Jetscape::real FluidCellInfo::pressure
 * Thermal pressure [GeV/fm^3].
 *
 * @var Jetscape::real FluidCellInfo::qgp_fraction
 * Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
 *
 * @var Jetscape::real FluidCellInfo::mu_B
 * Net baryon chemical potential [GeV].
 *
 * @var Jetscape::real FluidCellInfo::mu_C
 * Net charge chemical potential [GeV].
 *
 * @var Jetscape::real FluidCellInfo::mu_S
 * Net strangeness chemical potential [GeV].
 *
 * @var Jetscape::real FluidCellInfo::vx
 * Flow velocity in the x direction.
 *
 * @var Jetscape::real FluidCellInfo::vy
 * Flow velocity in the y direction.
 *
 * @var Jetscape::real FluidCellInfo::vz
 * Flow velocity in the z direction.
 *
 * @var Jetscape::real FluidCellInfo::pi
 * Shear stress tensor [GeV/fm^3].
 *
 * @var Jetscape::real FluidCellInfo::bulk_Pi
 * Bulk viscous pressure [GeV/fm^3].
 *
 * @fn FluidCellInfo::FluidCellInfo()
 * Default constructor.
 *
 * @fn FluidCellInfo FluidCellInfo::operator*=(Jetscape::real b)
 * Multiply the fluid cell by scalar factor b.
 *
 * @fn void FluidCellInfo::Print()
 * Prints fluid cell properties to the screen.
 */
class FluidCellInfo {
 public:
  // data structure for outputing fluid cell information
  Jetscape::real energy_density;   //!< Local energy density [GeV/fm^3].
  Jetscape::real entropy_density;  //!< Local entropy density [1/fm^3].
  Jetscape::real temperature;      //!< Local temperature [GeV].
  Jetscape::real pressure;         //!< Thermal pressure [GeV/fm^3].
  Jetscape::real qgp_fraction;     //!< Fraction of quark gluon plasma assuming
                                   //!< medium is in QGP+HRG phase.
  Jetscape::real mu_B;             //!< Net baryon chemical potential [GeV].
  Jetscape::real mu_C;             //!< Net charge chemical potential [GeV]
  Jetscape::real mu_S;        //!< Net strangeness chemical potential [GeV].
  Jetscape::real vx, vy, vz;  //!< Flow velocity.
  Jetscape::real pi[4][4];    //!< Shear stress tensor [GeV/fm^3].
  Jetscape::real bulk_Pi;     //!< Bulk viscous pressure [GeV/fm^3].

  /** Default constructor.*/
  FluidCellInfo();

  /** @param b Multiply the fluid cell by scalar factor b. */
  FluidCellInfo inline operator*=(Jetscape::real b);

  /** Prints fluid cell properties to the screen. */
  // void Print();
};

/**
 * @brief Overloads the `+` operator for FluidCellInfo to facilitate linear interpolation.
 *
 * This operator adds corresponding attributes of two FluidCellInfo objects.
 * The operation is performed element-wise on all scalar and matrix attributes.
 *
 * @param a The first FluidCellInfo object (copied by value).
 * @param b The second FluidCellInfo object (const reference).
 * @return The result of element-wise addition of `a` and `b`.
 *
 * @note The input `a` is modified during the addition before returning.
 */
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

/**
 * @brief Overloaded multiplication assignment operator for FluidCellInfo.
 *
 * This operator scales all the properties of the FluidCellInfo object by a given scalar value.
 *
 * @param b The scalar value by which to multiply the properties of the FluidCellInfo object.
 * @return A reference to the modified FluidCellInfo object.
 */
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

/**
 * @brief Multiplies a scalar value with a FluidCellInfo object.
 *
 * This operator overload allows for the multiplication of a scalar value
 * (of type Jetscape::real) with a FluidCellInfo object. The multiplication
 * is performed by scaling the FluidCellInfo object by the scalar value.
 *
 * @param a The scalar value to multiply with the FluidCellInfo object.
 * @param b The FluidCellInfo object to be scaled.
 * @return A new FluidCellInfo object that is the result of the multiplication.
 */
inline FluidCellInfo operator*(Jetscape::real a, FluidCellInfo b) {
  b *= a;
  return b;
}

/**
 * @brief Multiplies a FluidCellInfo object by a scalar value.
 *
 * This operator overloads the multiplication operator for the FluidCellInfo class,
 * allowing a FluidCellInfo object to be multiplied by a scalar value of type Jetscape::real.
 *
 * @param a The FluidCellInfo object to be multiplied.
 * @param b The scalar value of type Jetscape::real to multiply with.
 * @return A new FluidCellInfo object that is the result of the multiplication.
 */
inline FluidCellInfo operator*(FluidCellInfo a, Jetscape::real b) {
  a *= b;
  return a;
}

/**
 * @brief Overloads the division operator for FluidCellInfo.
 *
 * This function allows dividing a FluidCellInfo object by a scalar value of type Jetscape::real.
 * It modifies the FluidCellInfo object by multiplying it with the reciprocal of the scalar value.
 *
 * @param a The FluidCellInfo object to be divided.
 * @param b The scalar value of type Jetscape::real by which the FluidCellInfo object is divided.
 * @return A new FluidCellInfo object that is the result of the division.
 */
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
