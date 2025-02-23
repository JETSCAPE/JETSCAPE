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

#ifndef FOURVECTOR_H
#define FOURVECTOR_H

#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "JetScapeConstants.h"

using std::cerr;
using std::cout;
using std::endl;
using std::sqrt;

namespace Jetscape {

/**
 * @class FourVector
 * @brief Represents a four-vector with time and spatial components.
 */
class FourVector {
 public:
  /**
   * @brief Default constructor initializing all components to zero.
   */
  FourVector() : tv(0.0), xv(0.0), yv(0.0), zv(0.0){};

  /**
   * @brief Copy constructor.
   * @param srv FourVector to copy from.
   */
  FourVector(const FourVector &srv)
      : xv(srv.xv), yv(srv.yv), zv(srv.zv), tv(srv.tv){};  // copy constructor

  /**
   * @brief Constructor with an array input.
   * @param a Array of four elements representing (t, x, y, z).
   */
  FourVector(double a[4])  // constructor with array input
  {
    tv = a[0];
    xv = a[1];
    yv = a[2];
    zv = a[3];
  };

  /**
   * @brief Constructor with individual component inputs.
   * @param x_in x component.
   * @param y_in y component.
   * @param z_in z component.
   * @param t_in t component.
   */
  FourVector(double x_in, double y_in, double z_in, double t_in) {
    tv = t_in;
    xv = x_in;
    yv = y_in;
    zv = z_in;
  };

  /**
   * @brief Clears the four-vector, setting all components to zero.
   */
  void clear() { tv = xv = yv = zv = 0.0; }

  // constructors do all sets

  /**
   * @brief Sets the four-vector components.
   * @param x_in x component.
   * @param y_in y component.
   * @param z_in z component.
   * @param t_in t component.
   */
  void Set(double x_in, double y_in, double z_in, double t_in) {
    tv = t_in;
    xv = x_in;
    yv = y_in;
    zv = z_in;
  }

  /**
   * @brief Sets the four-vector components using an array.
   * @param a Array of four elements representing (t, x, y, z).
   */
  void Set(double a[4]) {
    tv = a[0];
    xv = a[1];
    yv = a[2];
    zv = a[3];
  };

  // all gets are done with name calls e.g., vec.x()

  /**
   * @brief Returns the x component.
   * @return x component.
   */
  double x() const { return (xv); };

  /**
   * @brief Returns the y component.
   * @return y component.
   */
  double y() const { return (yv); };

  /**
   * @brief Returns the z component.
   * @return z component.
   */
  double z() const { return (zv); };

  /**
   * @brief Returns the t component.
   * @return t component.
   */
  double t() const { return (tv); };

  /**
   * @brief Retrieves a component value based on the given index.
   *
   * This function returns one of four possible values (tv, xv, yv, zv)
   * based on the provided index. If the index is out of range (not between 0
   * and 3), it prints an error message and returns a large number.
   *
   * @param i The index of the component to retrieve (valid values: 0 to 3).
   * @return The corresponding component value:
   *         - 0 -> tv
   *         - 1 -> xv
   *         - 2 -> yv
   *         - 3 -> zv
   *         - Otherwise, returns a large number.
   */
  const double comp(int i) const {
    switch (i) {
      case 0:
        return (tv);
        break;
      case 1:
        return (xv);
        break;
      case 2:
        return (yv);
        break;
      case 3:
        return (zv);
        break;
      default:
        cout << " component index beyond 0-3! Returning garbage ..." << endl;
        return (a_very_large_number);
        break;
    }
  }

  /**
   * @brief Computes the plus component of the four-vector.
   *
   * This function calculates (zv + tv) / sqrt(2.0),
   * which is commonly used in light-cone coordinates.
   *
   * @return The computed plus component.
   */
  double plus() { return ((zv + tv) / sqrt(2.0)); };

  /**
   * @brief Computes the minus component of the four-vector.
   *
   * This function calculates (tv - zv) / sqrt(2.0),
   * which is commonly used in light-cone coordinates.
   *
   * @return The computed minus component.
   */
  double minus() { return ((tv - zv) / sqrt(2.0)); };

  /**
   * @brief Computes the rapidity of the particle.
   *
   * The rapidity is defined as (1/2) * log(plus/minus) when minus > 0.
   * If the minus component is non-positive, an error is printed.
   *
   * @return The rapidity value if valid; otherwise, returns 0.
   */
  double rapidity() {
    if (this->minus() > 0.0)
      return (std::log(this->plus() / this->minus()) / 2.0);
    cout << endl
         << "ERROR: z component exceeds t component, cannot calculate rapidity"
         << endl;
    return (0);
  };

  /**
   * @brief Computes the pseudorapidity (eta) of the particle.
   *
   * If the particle is moving strictly in the z direction, the function
   * returns a very large number or zero based on the sign of zv.
   * Otherwise, it calculates eta using the formula:
   * eta = (1/2) * log((v + zv) / (v - zv)), where v is the magnitude
   * of the three-momentum.
   *
   * @return The computed pseudorapidity value.
   */
  double eta() {
    if ((xv == 0) && (yv == 0)) {
      cout << " particle strictly in z direction " << endl;
      if (zv > 0)
        return (a_very_large_number);
      if (zv < 0)
        return (-1 * a_very_large_number);
      if (zv == 0)
        return (0);
    }

    double v = sqrt(xv * xv + yv * yv + zv * zv);

    double eta = std::log((v + zv) / (v - zv)) / 2.0;

    return (eta);
  }

  /**
   * @brief Computes the azimuthal angle (φ) in the XY-plane.
   *
   * The function calculates the angle φ in radians using the `atan2` function,
   * ensuring it falls within the range [0, 2π]. If both `x()` and `y()` are
   * approximately zero (within `rounding_error`), the function returns 0.
   *
   * @return The azimuthal angle φ in radians, within the range [0, 2π].
   */
  double phi() {
    if (fabs(x()) < rounding_error && fabs(y()) < rounding_error) {
      return 0;
    }
    double phi = atan2(y(), x());
    while (phi < 0)
      phi += 2.0 * pi;
    return phi;
  };

  /**
   * @brief Computes the Minkowski inner product of two FourVectors.
   *
   * @param c The FourVector to be multiplied.
   * @return The scalar result of the Minkowski inner product.
   */
  double operator*(FourVector &c) {
    return (tv * c.t() - xv * c.x() - yv * c.y() - zv * c.z());
  };

  /**
   * @brief Adds another FourVector to this FourVector.
   *
   * @param c The FourVector to add.
   * @return A reference to the updated FourVector.
   */
  FourVector &operator+=(FourVector &c) {
    tv += c.t();
    xv += c.x();
    yv += c.y();
    zv += c.z();

    return (*this);
  };

  /**
   * @brief Subtracts another FourVector from this FourVector.
   *
   * @param c The FourVector to subtract.
   * @return A reference to the updated FourVector.
   */
  FourVector &operator-=(FourVector &c) {
    tv -= c.t();
    xv -= c.x();
    yv -= c.y();
    zv -= c.z();

    return (*this);
  };

  /**
   * @brief Assigns values from another FourVector to this FourVector.
   *
   * @param c The FourVector to copy from.
   * @return A reference to the updated FourVector.
   */
  FourVector &operator=(FourVector &c) {
    tv = c.t();
    xv = c.x();
    yv = c.y();
    zv = c.z();
    return (*this);
  };

  /**
   * @brief Assigns values from another FourVector to this FourVector (const
   * version).
   *
   * @param c The FourVector to copy from.
   * @return A reference to the updated FourVector.
   */
  FourVector &operator=(const FourVector &c) {
    tv = c.tv;
    xv = c.xv;
    yv = c.yv;
    zv = c.zv;
    return (*this);
  };

  /**
   * @brief Rotates the FourVector around the x-axis by a given angle.
   *
   * @param theta The angle in radians by which to rotate around the x-axis.
   */
  void rotate_around_x(double theta) {
    double new_zv, new_yv;

    new_yv = yv * cos(theta) - zv * sin(theta);
    new_zv = zv * cos(theta) + yv * sin(theta);

    zv = new_zv;
    yv = new_yv;
  };

  /**
   * @brief Rotates the object around the Z-axis by a given angle.
   *
   * This function updates the `xv` and `yv` coordinates based on a
   * counterclockwise rotation about the Z-axis by the angle `theta`.
   *
   * @param theta The rotation angle in radians.
   */
  void rotate_around_z(double theta) {
    double new_xv, new_yv;

    new_xv = xv * cos(theta) - yv * sin(theta);
    new_yv = yv * cos(theta) + xv * sin(theta);

    xv = new_xv;
    yv = new_yv;
  };

  /**
   * @brief Rotates the object around the Y-axis by a given angle.
   *
   * This function updates the `xv` and `zv` coordinates based on a
   * counterclockwise rotation about the Y-axis by the angle `theta`.
   *
   * @param theta The rotation angle in radians.
   */
  void rotate_around_y(double theta) {
    double new_zv, new_xv;

    new_zv = zv * cos(theta) - xv * sin(theta);
    new_xv = xv * cos(theta) + zv * sin(theta);

    xv = new_xv;
    zv = new_zv;
  };

  /**
   * @brief Applies a Lorentz boost in the eta direction.
   *
   * This function modifies the time (tv) and spatial (zv) coordinates
   * using a hyperbolic transformation based on the provided rapidity
   * difference.
   *
   * @param deta The rapidity difference for the boost.
   */
  void eta_boost(double deta) {
    double new_zv, new_tv;

    new_tv = tv * cosh(deta) - zv * sinh(deta);
    new_zv = zv * cosh(deta) - tv * sinh(deta);

    tv = new_tv;
    zv = new_zv;
  };

  /**
   * @brief Applies a Lorentz boost in the rapidity (y) direction.
   *
   * This function updates the time (tv) and spatial (zv) coordinates
   * using a hyperbolic transformation with the given rapidity shift.
   *
   * @param dy The rapidity shift for the boost.
   */
  void y_boost(double dy) {
    double new_zv, new_tv;

    new_tv = tv * cosh(dy) - zv * sinh(dy);
    new_zv = zv * cosh(dy) - tv * sinh(dy);

    tv = new_tv;
    zv = new_zv;
  };

  /**
   * @brief Applies a Lorentz boost with a given velocity vector.
   *
   * This function modifies the four-vector (tv, xv, yv, zv) using a boost
   * transformation along the velocity components (vx, vy, vz). If the velocity
   * magnitude is greater than or equal to 1 (speed of light), a large gamma
   * value is assigned, and a warning message is displayed.
   *
   * @param vx Velocity component in the x-direction.
   * @param vy Velocity component in the y-direction.
   * @param vz Velocity component in the z-direction.
   */
  void boost(double vx, double vy, double vz) {
    double gamma, v;

    v = sqrt(vx * vx + vy * vy + vz * vz);
    if (v < 1) {
      gamma = 1 / sqrt(1 - v * v);
    } else {
      gamma = a_very_large_number;
      std::cout << " light like boost, setting gamma = " << gamma << std::endl;
    };

    double x = xv;
    double y = yv;
    double z = zv;
    double t = tv;

    tv = gamma * t - gamma * (vx * x + vy * y + vz * z);
    xv = -1 * gamma * vx * t + x + (gamma - 1) * vx * vx * x / (v * v) +
         (gamma - 1) * vx * vy * y / (v * v) +
         (gamma - 1) * vx * vz * z / (v * v);
    yv = -1 * gamma * vy * t + y + (gamma - 1) * vy * vy * y / (v * v) +
         (gamma - 1) * vy * vx * x / (v * v) +
         (gamma - 1) * vy * vz * z / (v * v);
    zv = -1 * gamma * vz * t + z + (gamma - 1) * vz * vz * z / (v * v) +
         (gamma - 1) * vz * vx * x / (v * v) +
         (gamma - 1) * vz * vy * y / (v * v);

    double old_inv = t * t - x * x - y * y - z * z;
    double new_inv = tv * tv - xv * xv - yv * yv - zv * zv;

    if (std::abs(old_inv - new_inv) > 1 / a_very_large_number)
      std::cout << " invariants dont match after boost " << std::endl;
  }

 private:
  /**
   * @brief Private member variables representing coordinates and time.
   *
   * The 'v'is for vector.  We call the private variables, xv, tv etc., so that
   * getter functions will be called x, t etc.
   */
  double xv;  ///< X-coordinate component.
  double yv;  ///< Y-coordinate component.
  double zv;  ///< Z-coordinate component.
  double tv;  ///< Time component.
};

};  // namespace Jetscape

#endif  // FOURVECTOR_H
