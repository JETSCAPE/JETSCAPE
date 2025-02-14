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

class FourVector {
  // the class of four vectors

 public:
  FourVector()  // default constructor
  {
    tv = xv = yv = zv = 0.0;
  };

  FourVector(const FourVector &srv)
      : xv(srv.xv), yv(srv.yv), zv(srv.zv), tv(srv.tv){};  // copy constructor

  FourVector(double a[4])  // constructor with array input
  {
    tv = a[0];
    xv = a[1];
    yv = a[2];
    zv = a[3];
  };

  FourVector(double x_in, double y_in, double z_in,
             double t_in)  // constructor with component input
  {
    tv = t_in;
    xv = x_in;
    yv = y_in;
    zv = z_in;
  };

  void clear() { tv = xv = yv = zv = 0.0; }

  // constructors do all sets

  void Set(double x_in, double y_in, double z_in, double t_in) {
    tv = t_in;
    xv = x_in;
    yv = y_in;
    zv = z_in;
  }

  void Set(double a[4]) {
    tv = a[0];
    xv = a[1];
    yv = a[2];
    zv = a[3];
  };

  // all gets are done with name calls e.g., vec.x()
  double x() const { return (xv); };

  double y() const { return (yv); };

  double z() const { return (zv); };

  double t() const { return (tv); };

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

  double plus() { return ((zv + tv) / sqrt(2.0)); };

  double minus() { return ((tv - zv) / sqrt(2.0)); };

  double rapidity() {
    if (this->minus() > 0.0)
      return (std::log(this->plus() / this->minus()) / 2.0);
    cout << endl
         << "ERROR: z component exceeds t component, cannot calculate rapidity"
         << endl;
    return (0);
  };

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

  double phi() {
    if (fabs(x()) < rounding_error && fabs(y()) < rounding_error) {
      return 0;
    }
    double phi = atan2(y(), x());
    while (phi < 0)
      phi += 2.0 * pi;
    return phi;
  };

  double operator*(FourVector &c) {
    return (tv * c.t() - xv * c.x() - yv * c.y() - zv * c.z());
  };

  FourVector &operator+=(FourVector &c) {
    tv += c.t();
    xv += c.x();
    yv += c.y();
    zv += c.z();

    return (*this);
  };

  FourVector &operator-=(FourVector &c) {
    tv -= c.t();
    xv -= c.x();
    yv -= c.y();
    zv -= c.z();

    return (*this);
  };

  FourVector &operator=(FourVector &c) {
    tv = c.t();
    xv = c.x();
    yv = c.y();
    zv = c.z();
    return (*this);
  };

  FourVector &operator=(const FourVector &c) {
    tv = c.tv;
    xv = c.xv;
    yv = c.yv;
    zv = c.zv;
    return (*this);
  };

  void rotate_around_x(double theta) {
    double new_zv, new_yv;

    new_yv = yv * cos(theta) - zv * sin(theta);
    new_zv = zv * cos(theta) + yv * sin(theta);

    zv = new_zv;
    yv = new_yv;
  };

  void rotate_around_z(double theta) {
    double new_xv, new_yv;

    new_xv = xv * cos(theta) - yv * sin(theta);
    new_yv = yv * cos(theta) + xv * sin(theta);

    xv = new_xv;
    yv = new_yv;
  };

  void rotate_around_y(double theta) {
    double new_zv, new_xv;

    new_zv = zv * cos(theta) - xv * sin(theta);
    new_xv = xv * cos(theta) + zv * sin(theta);

    xv = new_xv;
    zv = new_zv;
  };

  void eta_boost(double deta) {
    double new_zv, new_tv;

    new_tv = tv * cosh(deta) - zv * sinh(deta);
    new_zv = zv * cosh(deta) - tv * sinh(deta);

    tv = new_tv;
    zv = new_zv;
  };

  void y_boost(double dy) {
    double new_zv, new_tv;

    new_tv = tv * cosh(dy) - zv * sinh(dy);
    new_zv = zv * cosh(dy) - tv * sinh(dy);

    tv = new_tv;
    zv = new_zv;
  };

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
  // the v is for vector, we call the private variables, xv, tv etc., so that
  // get function calls will be called x, t etc.
  double xv, yv, zv, tv;
};

};  // namespace Jetscape

#endif  // FOURVECTOR_H
