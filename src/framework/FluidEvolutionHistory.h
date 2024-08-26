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

#ifndef EVOLUTIONHISTORY_H
#define EVOLUTIONHISTORY_H

#include <stdexcept>
#include <vector>

#include "FluidCellInfo.h"
#include "RealType.h"

namespace Jetscape {

// The simplest way for 3rd party hydro to provide evolution history is to
// send 2 vectors to the framework. One is a 1D float vector stores all
// the evolution history to data_vector.
// The other is the description of the content of the data to data_info.
// E.g., data_info should store a vector of 3 strings ['energy_density', 'vx,
// 'vy'] if data_vector = [ed0, vx0, vy0, ed1, vx1, vy1, ..., edN, vxN, vyN,
// vzN], where N = ntau * nx * ny * netas is the number of total cells. One can
// pass these 3 data or more data to the framework, by storing them in any
// order, as long as the description matches the data content. The following is
// a list possible valid data names. ValidNames = [ "energy_density",
// "entropy_density", "temperature",
//                "pressure", "qgp_fraction", "mu_b", "mu_c", "mu_s",
//                "vx", "vy", "vz", "pi00", "pi01", "pi02", "pi03",
//                "pi11", "pi12", "pi13", "pi22", "pi23", "pi33", "bulk_pi"];
//

enum EntryName {
  ENTRY_ENERGY_DENSITY,
  ENTRY_ENTROPY_DENSITY,
  ENTRY_TEMPERATURE,
  ENTRY_PRESSURE,
  ENTRY_QGP_FRACTION,
  ENTRY_MU_B,
  ENTRY_MU_C,
  ENTRY_MU_S,
  ENTRY_VX,
  ENTRY_VY,
  ENTRY_VZ,
  ENTRY_PI00,
  ENTRY_PI01,
  ENTRY_PI02,
  ENTRY_PI03,
  ENTRY_PI11,
  ENTRY_PI12,
  ENTRY_PI13,
  ENTRY_PI22,
  ENTRY_PI23,
  ENTRY_PI33,
  ENTRY_BULK_PI,
  ENTRY_INVALID
};

EntryName ResolveEntryName(std::string input);

class InvalidSpaceTimeRange : public std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

class EvolutionHistory {
 public:
  /** @param tau_min Minimum value of tau.*/
  Jetscape::real tau_min, dtau;  //!< @param dtau Step-size for tau.
  /** @param x_min Minimum value of x. */
  Jetscape::real x_min, dx;  //!< @param dx Step-size for x.
  /** @param y_min Minimum value of y. */
  Jetscape::real y_min, dy;  //!< @param dy Step-size for y.
  /** @param eta_min Minimum value of eta. */
  Jetscape::real eta_min, deta;  //!< @param deta Step-size for eta.
  int ntau;  //!< @param ntau Number of grid points in tau-axis.
  int nx;    //!< @param nx Number of grid points in x-axis.
  int ny;    //!< @param ny Number of grid points in y-axis.
  int neta;  //!< @param neta Number of grid points in eta-axis.

  /** Default is set to false. Set flag tau_eta_is_tz to true if hydro dynamics
   * is setup in (t,x,y,z) coordinate. */
  bool tau_eta_is_tz;

  bool boost_invariant;

  /** The bulk information of hydro dynamics, in the form of
   * vector of FluidCellInfo objects. */
  std::vector<FluidCellInfo> data;

  /** The bulk information of hydro dynamics, in the form of 1D float vector.
   * It is easy to pass vector of float from 3rd party hydro module
   * to Jetscape Fluid Evolution History, where bulk info should be stored
   * in orders of ed0, sd0, temp0, ..., ed1, sd1, temp1, ..., edn, sdn, tempn.
   * The content and order of data entres are given by data_info */
  std::vector<float> data_vector;

  /** Store the entry names of one record in the data array*/
  std::vector<std::string> data_info;

  /** Default constructor. */
  EvolutionHistory() = default;

  /** Auxiliary function that helps to read hydro evolution history from
   * external fluid dynamic module. The history from 3rd party hydro should be
   * stored in a std::vector<float> container. The size of the vector should
   * equal ntau * nx * ny * neta * data_info_.size(); For each piece with length
   * data_info_.size(), the vector should store the corresponding float number
   * whose name is written in std::vector<string> data_info_.
   * E.g., data_info_ can be a vector of ['
   * */
  void FromVector(const std::vector<float> &data_,
                  const std::vector<std::string> &data_info_, float tau_min,
                  float dtau, float x_min, float dx, int nx, float y_min,
                  float dy, int ny, float eta_min, float deta, int neta,
                  bool tau_eta_is_tz);

  /** Default destructor. */
  ~EvolutionHistory() {
    data.clear();
    data_vector.clear();
    data_info.clear();
  }

  void clear_up_evolution_data() { data.clear(); }

  int get_data_size() const { return (data.size()); }
  bool is_boost_invariant() const { return (boost_invariant); }
  bool is_Cartesian() const { return (tau_eta_is_tz); }

  Jetscape::real Tau0() const { return (tau_min); }
  Jetscape::real XMin() const { return (x_min); }
  Jetscape::real YMin() const { return (y_min); }
  Jetscape::real EtaMin() const { return (eta_min); }

  /** Maximum value of tau. */
  inline Jetscape::real TauMax() const { return (tau_min + (ntau - 1) * dtau); }

  /** Maximum value of x. */
  inline Jetscape::real XMax() const { return (x_min + (nx - 1) * dx); }

  /** Maximum value of y. */
  inline Jetscape::real YMax() const { return (y_min + (ny - 1) * dy); }

  /** Maximum value of eta. */
  inline Jetscape::real EtaMax() const { return (eta_min + (neta - 1) * deta); }

  /** It checks whether a space-time point (tau, x, y, eta) is inside evolution
     history or outside.
        @param tau Light-cone coordinate.
        @param x  Space coordinate.
        @param y  Space coordinate.
        @param eta Light-cone coordinate.
    */
  int CheckInRange(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                   Jetscape::real eta) const;

  // get the lower bound of the fluid cell along tau
  /** @return Fluid cell number along the tau-grid.
        @param tau Light-cone coordinate.
    */
  inline int GetIdTau(Jetscape::real tau) const {
    return (static_cast<int>((tau - tau_min) / dtau));
  }

  // get the lower bound of the fluid cell along x
  /** @return Fluid cell number along the x-grid.
        @param x Space coordinate.
    */
  inline int GetIdX(Jetscape::real x) const {
    return (static_cast<int>((x - x_min) / dx));
  }

  // get the lower bound of the fluid cell along y
  /** @return Fluid cell number along the y-grid.
        @param y Space coordinate.
    */
  inline int GetIdY(Jetscape::real y) const {
    return (static_cast<int>((y - y_min) / dy));
  }

  // get the lower bound of the fluid cell along eta
  /** @return Fluid cell number along the eta-grid.
        @param eta Light-cone coordinate.
    */
  inline int GetIdEta(Jetscape::real eta) const {
    return (static_cast<int>((eta - eta_min) / deta));
  }

  // get the coordinate of tau, x, y, eta on grid
  /** @param id_tau Fluid cell number along tau-grid.
        @return The tau coordinate for fluid cell number.
    */
  inline Jetscape::real TauCoord(int id_tau) const {
    return (tau_min + id_tau * dtau);
  }

  /** @param id_x Fluid cell number along x-grid.
        @return The x coordinate for fluid cell number.
    */
  inline Jetscape::real XCoord(int id_x) const { return (x_min + id_x * dx); }

  /** @param id_y Fluid cell number along y-grid.
        @return The y coordinate for fluid cell number.
    */
  inline Jetscape::real YCoord(int id_y) const { return (y_min + id_y * dy); }

  /** @param id_eta Fluid cell number along eta-grid.
        @return The eta coordinate for fluid cell number.
    */
  inline Jetscape::real EtaCoord(int id_eta) const {
    return (eta_min + id_eta * deta);
  }

  // get the FluidCellInfo index in data
  /** @return FluidCellInfo index in the data.
        @param id_tau Fluid cell number along tau-grid.
        @param id_x Fluid cell number along x-grid.
        @param id_y Fluid cell number along y-grid.
        @param id_eta Fluid cell number along eta-grid.
    */
  inline int CellIndex(int id_tau, int id_x, int id_y, int id_eta) const {
    id_tau = std::min(ntau - 1, std::max(0, id_tau));
    id_x = std::min(nx - 1, std::max(0, id_x));
    id_y = std::min(ny - 1, std::max(0, id_y));
    id_eta = std::min(neta - 1, std::max(0, id_eta));
    return (id_tau * nx * ny * neta + id_x * ny * neta + id_y * neta + id_eta);
  }

  /* Read fluid cell info for a given lattice cell*/
  FluidCellInfo GetFluidCell(int id_tau, int id_x, int id_y, int id_eta) const;

  // get the FluidCellInfo at space point given time step
  /** @return FluidCellInfo at a point (x,y,eta) and time-step id_tau.
        @param id_tau tau-step number.
        @param x Space coordinate.
        @param y Space coordinate.
        @param eta Light-cone coordinate.
    */
  FluidCellInfo GetAtTimeStep(int id_tau, Jetscape::real x, Jetscape::real y,
                              Jetscape::real etas) const;

  // get the FluidCellInfo at given space time point
  /** @return FluidCellInfo at a point (tau, x, y, eta).
        @param tau Light-cone coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param eta Light-cone coordinate.
    */
  FluidCellInfo get(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                    Jetscape::real etas) const;
  FluidCellInfo get_tz(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                       Jetscape::real z) const;
};

}  // namespace Jetscape

#endif  // EVOLUTIONHISTORY_H
