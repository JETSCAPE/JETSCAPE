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

/**
 * @note
 * The simplest way for 3rd party hydro to provide evolution history is to
 * send 2 vectors to the framework. One is a 1D float vector stores all
 * the evolution history to data_vector.
 * The other is the description of the content of the data to data_info.
 * E.g., data_info should store a vector of 3 strings ['energy_density', 'vx,
 * 'vy'] if data_vector = [ed0, vx0, vy0, ed1, vx1, vy1, ..., edN, vxN, vyN,
 * vzN], where N = ntau * nx * ny * netas is the number of total cells. One can
 * pass these 3 data or more data to the framework, by storing them in any
 * order, as long as the description matches the data content. The following is
 * a list possible valid data names. ValidNames = [ "energy_density",
 * "entropy_density", "temperature",
 *                "pressure", "qgp_fraction", "mu_b", "mu_c", "mu_s",
 *                "vx", "vy", "vz", "pi00", "pi01", "pi02", "pi03",
 *                "pi11", "pi12", "pi13", "pi22", "pi23", "pi33", "bulk_pi"];
 */

/**
 * @enum EntryName
 * @brief Enumeration of various physical quantities in a dataset.
 *
 * This enumeration defines different types of entries that can be present in
 * a dataset related to fluid dynamics or high-energy physics simulations.
 */
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

/**
 * @brief Converts a string representation of an entry name to the corresponding
 * enum type EntryName.
 *
 * This function takes a string as input and attempts to resolve it to a
 * corresponding EntryName enum value. If the input string matches a known entry
 * name, the associated enum value is returned. If no match is found,
 * ENTRY_INVALID is returned.
 *
 * @param input A string representing the entry name.
 * @return The corresponding EntryName enum value, or ENTRY_INVALID if the input
 * is not recognized.
 */
EntryName ResolveEntryName(std::string input);

/**
 * @class InvalidSpaceTimeRange
 * @brief Exception class to represent an invalid space-time range error.
 *
 * This class inherits from `std::invalid_argument` and is used to indicate
 * that a space-time range is invalid.
 *
 * @note This class does not add any new functionality to the base class, but it
 *       provides a more specific exception type that can be caught in exception
 *       handling code.
 *
 * @see std::invalid_argument
 */
class InvalidSpaceTimeRange : public std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

class EvolutionHistory {
 public:
  /** Minimum value of tau. */
  Jetscape::real tau_min;

  /** Step-size for tau. */
  Jetscape::real dtau;

  /** Minimum value of x. */
  Jetscape::real x_min;
  /** Step-size for x. */
  Jetscape::real dx;

  /** Minimum value of y. */
  Jetscape::real y_min;
  /** Step-size for y. */
  Jetscape::real dy;

  /** Minimum value of eta. */
  Jetscape::real eta_min;
  /** Step-size for eta. */
  Jetscape::real deta;

  /** Number of grid points in tau-axis. */
  int ntau;

  /** Number of grid points in x-axis. */
  int nx;

  /** Number of grid points in y-axis. */
  int ny;

  /** Number of grid points in eta-axis. */
  int neta;

  /**
   * Default is set to false.
   * Set flag tau_eta_is_tz to true if hydro dynamics
   * is setup in (t,x,y,z) coordinate.
   */
  bool tau_eta_is_tz;

  /** Boost_invariant flag. */
  bool boost_invariant;

  /**
   * The bulk information of hydro dynamics, in the form of
   * vector of FluidCellInfo objects. */
  std::vector<FluidCellInfo> data;

  /**
   * The bulk information of hydro dynamics, in the form of 1D float vector.
   * It is easy to pass vector of float from 3rd party hydro module
   * to Jetscape Fluid Evolution History, where bulk info should be stored
   * in orders of ed0, sd0, temp0, ..., ed1, sd1, temp1, ..., edn, sdn, tempn.
   * The content and order of data entres are given by data_info */
  std::vector<float> data_vector;

  /** Store the entry names of one record in the data array*/
  std::vector<std::string> data_info;

  /** Default constructor. */
  EvolutionHistory() = default;

  /**
   * @brief Reads hydro evolution history from an external fluid dynamic module.
   *
   * This function processes evolution history data stored in a
   * std::vector<float>. The data must be structured such that its size equals
   * @f$ ntau \times nx \times ny \times neta \times data\_info\_.size() @f$.
   * Each segment of length `data_info_.size()` contains float values
   * corresponding to the names in `data_info_`.
   *
   * @param data_ The input vector containing hydro evolution history data.
   * @param data_info_ A vector of strings indicating the meaning of each data
   * field.
   * @param tau_min Minimum proper time value.
   * @param dtau Proper time step size.
   * @param x_min Minimum x-coordinate.
   * @param dx Step size in x-direction.
   * @param nx Number of grid points in the x-direction.
   * @param y_min Minimum y-coordinate.
   * @param dy Step size in y-direction.
   * @param ny Number of grid points in the y-direction.
   * @param eta_min Minimum space-time rapidity.
   * @param deta Step size in the eta direction.
   * @param neta Number of grid points in the eta direction.
   * @param tau_eta_is_tz Flag indicating whether tau-eta coordinates are used.
   */
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

  /**
   * @brief Clear the evolution history data.
   */
  void clear_up_evolution_data() { data.clear(); }

  /**
   * @brief Get the size of the data vector.
   * @return The size of the data vector.
   */
  int get_data_size() const { return (data.size()); }

  /**
   * @brief Check the boost invariance flag.
   * @return True if the boost invariance flag is set, false otherwise.
   */
  bool is_boost_invariant() const { return (boost_invariant); }

  /**
   * @brief Check if the evolution history is in Cartesian coordinates.
   * @return True if the evolution history is in Cartesian coordinates, false
   * otherwise.
   */
  bool is_Cartesian() const { return (tau_eta_is_tz); }

  /**
   * @brief Gets the minimum value of tau.
   * @return The minimum value of tau.
   */
  Jetscape::real Tau0() const { return (tau_min); }

  /**
   * @brief Gets the minimum value of x.
   * @return The minimum value of x.
   */
  Jetscape::real XMin() const { return (x_min); }

  /**
   * @brief Gets the minimum value of y.
   * @return The minimum value of y.
   */
  Jetscape::real YMin() const { return (y_min); }

  /**
   * @brief Gets the minimum value of eta.
   * @return The minimum value of eta.
   */
  Jetscape::real EtaMin() const { return (eta_min); }

  /**
   * @brief Maximum value of tau.
   * @return The maximum value of tau.
   */
  inline Jetscape::real TauMax() const { return (tau_min + (ntau - 1) * dtau); }

  /**
   * @brief Calculates the maximum value of x.
   * @return The maximum value of x as a Jetscape::real type.
   */
  inline Jetscape::real XMax() const { return (x_min + (nx - 1) * dx); }

  /**
   * @brief Returns the maximum value of y.
   * @return The maximum value of y.
   */
  inline Jetscape::real YMax() const { return (y_min + (ny - 1) * dy); }

  /**
   * @brief Computes the maximum value of eta.
   * @return The maximum value of eta.
   */
  inline Jetscape::real EtaMax() const { return (eta_min + (neta - 1) * deta); }

  /**
   * @brief Checks whether a given space-time point is within the evolution
   * history range.
   *
   * This function determines if the provided light-cone coordinates (tau, eta)
   * and spatial coordinates (x, y) fall within the defined limits of the
   * evolution history. If any of the coordinates are out of bounds, a warning
   * message is generated, and the function returns 0. Otherwise, it returns 1.
   *
   * @param tau Light-cone coordinate representing proper time.
   * @param x   Spatial coordinate in the x-direction.
   * @param y   Spatial coordinate in the y-direction.
   * @param eta Light-cone coordinate representing space-time rapidity.
   * @return int Returns 1 if the point is within range, 0 otherwise.
   */
  int CheckInRange(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                   Jetscape::real eta) const;

  /**
   * @brief Gets the lower bound of the fluid cell along the tau-grid.
   *
   * @param tau The light-cone coordinate.
   * @return The fluid cell number along the tau-grid.
   */
  inline int GetIdTau(Jetscape::real tau) const {
    return (static_cast<int>((tau - tau_min) / dtau));
  }

  /**
   * @brief Gets the lower bound of the fluid cell along the x-axis.
   *
   * This function computes the fluid cell number along the x-grid based on a
   * given space coordinate.
   *
   * @param x The space coordinate along the x-axis.
   * @return The fluid cell number along the x-axis.
   */
  inline int GetIdX(Jetscape::real x) const {
    return (static_cast<int>((x - x_min) / dx));
  }

  /**
   * @brief Get the lower bound of the fluid cell along the y-axis.
   *
   * This function computes the fluid cell number along the y-grid based on the
   * provided space coordinate, using the minimum y-coordinate and the cell
   * size.
   *
   * @param y Space coordinate along the y-axis.
   *
   * @return Fluid cell number along the y-grid.
   */
  inline int GetIdY(Jetscape::real y) const {
    return (static_cast<int>((y - y_min) / dy));
  }

  /**
   * @brief Get the lower bound of the fluid cell along the eta grid.
   *
   * This function computes the fluid cell number along the eta-grid based on
   * the provided light-cone coordinate.
   *
   * @param eta The light-cone coordinate.
   *
   * @return The fluid cell number along the eta grid.
   */
  inline int GetIdEta(Jetscape::real eta) const {
    return (static_cast<int>((eta - eta_min) / deta));
  }

  /**
   * @brief Get the coordinate of tau on the grid.
   *
   * This function calculates the tau coordinate for a given fluid cell number
   * along the tau-grid.
   *
   * @param id_tau Fluid cell number along the tau-grid.
   * @return The tau coordinate for the given fluid cell number.
   */
  inline Jetscape::real TauCoord(int id_tau) const {
    return (tau_min + id_tau * dtau);
  }

  /**
   * @brief Calculates the x-coordinate for a given fluid cell number along the
   * x-grid.
   *
   * This function computes the x-coordinate corresponding to the given fluid
   * cell number along the x-axis.
   *
   * @param id_x The fluid cell number along the x-grid.
   * @return The x-coordinate for the given fluid cell number.
   */
  inline Jetscape::real XCoord(int id_x) const { return (x_min + id_x * dx); }

  /**
   * @brief Calculates the y-coordinate for a given fluid cell number.
   *
   * This function returns the y-coordinate corresponding to the fluid cell
   * number along the y-grid, based on the minimum y-coordinate and the grid
   * spacing.
   *
   * @param id_y Fluid cell number along the y-grid.
   *
   * @return The y-coordinate for the given fluid cell number.
   */
  inline Jetscape::real YCoord(int id_y) const { return (y_min + id_y * dy); }

  /**
   * @brief Computes the eta coordinate for a given fluid cell number.
   *
   * This function calculates the eta coordinate corresponding to the fluid cell
   * number along the eta-grid based on the minimum eta value and the grid
   * spacing.
   *
   * @param id_eta The fluid cell number along the eta-grid.
   * @return The eta coordinate for the specified fluid cell number.
   */
  inline Jetscape::real EtaCoord(int id_eta) const {
    return (eta_min + id_eta * deta);
  }

  /**
   * @brief Calculates the index of the FluidCellInfo in the data array.
   *
   * This function computes the index of a fluid cell based on its coordinates
   * along the tau, x, y, and eta grids, ensuring that the provided indices are
   * within bounds.
   *
   * @param id_tau The fluid cell number along the tau-grid.
   * @param id_x The fluid cell number along the x-grid.
   * @param id_y The fluid cell number along the y-grid.
   * @param id_eta The fluid cell number along the eta-grid.
   *
   * @return The index of the FluidCellInfo in the data array.
   */
  inline int CellIndex(int id_tau, int id_x, int id_y, int id_eta) const {
    id_tau = std::min(ntau - 1, std::max(0, id_tau));
    id_x = std::min(nx - 1, std::max(0, id_x));
    id_y = std::min(ny - 1, std::max(0, id_y));
    id_eta = std::min(neta - 1, std::max(0, id_eta));
    return (id_tau * nx * ny * neta + id_x * ny * neta + id_y * neta + id_eta);
  }

  /**
   * @brief Retrieves fluid cell information for a given lattice cell.
   *
   * This function reads the sparse data stored in `data_` with associated
   * information from `data_info_` to construct a `FluidCellInfo` object.
   * If the data is not stored in a sparse format, it retrieves the
   * corresponding `FluidCellInfo` object directly from `data_`.
   *
   * @param id_tau Temporal index of the cell.
   * @param id_x Spatial index along the x-axis.
   * @param id_y Spatial index along the y-axis.
   * @param id_eta Spatial index along the eta-axis. In 2+1D mode, this
   *               value is set to zero internally.
   * @return FluidCellInfo The reconstructed fluid cell information.
   */
  FluidCellInfo GetFluidCell(int id_tau, int id_x, int id_y, int id_eta) const;

  /**
   * @brief Retrieves the FluidCellInfo at a given spatial point and time step.
   *
   * This function fetches the fluid cell information at a specific
   * (x, y, eta) coordinate for a given time step id_tau. If boost invariance
   * is not assumed, it interpolates in the eta direction as well.
   *
   * @param id_tau The tau-step number (time step index).
   * @param x The x-coordinate of the spatial point.
   * @param y The y-coordinate of the spatial point.
   * @param eta The light-cone coordinate.
   * @return The FluidCellInfo at the specified (x, y, eta) location and time
   * step.
   */
  FluidCellInfo GetAtTimeStep(int id_tau, Jetscape::real x, Jetscape::real y,
                              Jetscape::real etas) const;

  /**
   * @brief Interpolates FluidCellInfo along the time direction.
   *
   * This function performs interpolation to retrieve the FluidCellInfo at a
   * specified space-time point. It interpolates along the time direction (tau)
   * between the two closest time steps. If the requested point is outside the
   * valid range, a default "zero" fluid cell is returned.
   *
   * @param tau The light-cone time coordinate to retrieve the fluid information
   * at.
   * @param x The spatial x-coordinate.
   * @param y The spatial y-coordinate.
   * @param eta The eta coordinate, which is a light-cone coordinate.
   * @return FluidCellInfo The interpolated fluid information at the given
   * space-time point.
   */
  FluidCellInfo get(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                    Jetscape::real etas) const;

  /**
   * @brief Computes the fluid cell information at a given spacetime point.
   *
   * This function calculates the fluid cell information for a given time,
   * spatial coordinates (x, y), and the longitudinal position (z) in spacetime.
   * It first checks if the point lies within the light cone (i.e., if t^2 >
   * z^2) and, if valid, computes the proper time (tau) and the rapidity (eta)
   * based on the provided coordinates. If the point is outside the light cone,
   * a warning is logged.
   *
   * @param t The time coordinate of the spacetime point.
   * @param x The spatial x-coordinate of the point.
   * @param y The spatial y-coordinate of the point.
   * @param z The longitudinal spatial coordinate of the point.
   *
   * @return A FluidCellInfo object containing the computed fluid cell data for
   * the given spacetime point.
   *
   * @note The function assumes the spacetime coordinates are consistent with
   * the relativistic framework and that the point is within the light cone for
   * valid results.
   */
  FluidCellInfo get_tz(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                       Jetscape::real z) const;
};

}  // namespace Jetscape

#endif  // EVOLUTIONHISTORY_H
