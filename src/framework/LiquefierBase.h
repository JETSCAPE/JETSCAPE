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

#ifndef LIQUEFIERBASE_H
#define LIQUEFIERBASE_H

#include <array>
#include <vector>

#include "FluidCellInfo.h"
#include "JetClass.h"
#include "RealType.h"
#include "sigslot.h"

namespace Jetscape {

/**
 * @brief Represents a localized energy-momentum contribution from a parton to the fluid medium.
 *
 * In the JETSCAPE framework, a Droplet is a conceptual object used to bridge the gap
 * between discrete partonic information and continuous hydrodynamic fields.
 * Each droplet represents a localized "chunk" of energy and momentum to be deposited
 * into the hydrodynamic grid. 
 */
class Droplet {
 private:
  std::array<Jetscape::real, 4> xmu; ///< Position 4-vector
  std::array<Jetscape::real, 4> pmu; ///< Momentum 4-vector

 public:
  /**
   * @brief Default constructor.
   *
   * Constructs an uninitialized droplet. The position and momentum vectors will
   * contain undefined values until explicitly set.
   */
  Droplet() = default;

  /**
   * @brief Construct a droplet from given position and momentum vectors.
   *
   * @param x_in The initial position 4-vector.
   * @param p_in The initial momentum 4-vector.
   */
  Droplet(std::array<Jetscape::real, 4> x_in,
          std::array<Jetscape::real, 4> p_in) {
    xmu = x_in;
    pmu = p_in;
  }

  /**
   * @brief Destructor.
   */
  ~Droplet(){};

  /**
   * @brief Get the position 4-vector of the droplet.
   *
   * @return A copy of the position vector.
   */
  std::array<Jetscape::real, 4> get_xmu() const { return (xmu); }

  /**
   * @brief Get the momentum 4-vector of the droplet.
   *
   * @return A copy of the momentum vector.
   */
  std::array<Jetscape::real, 4> get_pmu() const { return (pmu); }
};

/**
 * @brief Base class for converting partonic energy/momentum into hydrodynamic sources ("liquefying").
 */
class LiquefierBase {
 private:
  std::vector<Droplet> dropletlist; ///< List of droplets representing source contributions
  bool GetHydroCellSignalConnected; ///< Flag for whether signal connection to hydro exists
  
  const int drop_stat; ///< Droplet statistics
  const int miss_stat; ///< Missed parton statistics
  const int neg_stat; ///< Negative energy statistics
  const Jetscape::real hydro_source_abs_err; ///< Error tolerance for hydro sources
  bool threshold_energy_switch; ///< Whether to apply energy threshold filtering
  double e_threshold; ///< Energy threshold value

 public:
  /**
   * @brief Constructor.
   */
  LiquefierBase();

  /**
   * @brief Destructor that clears droplet list.
   */
  ~LiquefierBase() { Clear(); }

  /**
   * @brief Add a droplet to the internal list.
   * @param droplet_in Droplet to add
   */
  void add_a_droplet(Droplet droplet_in) { dropletlist.push_back(droplet_in); }

  /**
   * @brief Get number of droplet conversions performed.
   * @return Droplet statistic count
   */
  int get_drop_stat() const { return (drop_stat); }

  /**
   * @brief Get number of partons missed in processing.
   * @return Missed statistic count
   */
  int get_miss_stat() const { return (miss_stat); }

  /**
   * @brief Get number of partons with negative energy.
   * @return Negative statistic count
   */
  int get_neg_stat() const { return (neg_stat); }

  /**
   * @brief Get a specific droplet by index.
   * @param idx Index of the droplet
   * @return Droplet at that index
   */
  Droplet get_a_droplet(const int idx) const { return (dropletlist[idx]); }

  /**
   * @brief Check energy-momentum conservation between partons.
   * @param pIn Input partons
   * @param pOut Output partons
   */
  void check_energy_momentum_conservation(const std::vector<Parton> &pIn,
                                          std::vector<Parton> &pOut);

  /**
   * @brief Apply filtering to remove partons based on criteria.
   * @param pOut Partons to be filtered in-place
   */
  void filter_partons(std::vector<Parton> &pOut);

  /**
   * @brief Add hydrodynamic sources based on input/output partons.
   * @param pIn Input partons
   * @param pOut Output partons
   */
  void add_hydro_sources(std::vector<Parton> &pIn, std::vector<Parton> &pOut);

  /**
   * @brief Signal used to query hydro cell information.
   */
  sigslot::signal5<double, double, double, double,
                   std::unique_ptr<FluidCellInfo> &,
                   sigslot::multi_threaded_local>
      GetHydroCellSignal;

  /**
   * @brief Check if hydro cell signal is connected.
   * @return True if signal is connected
   */
  const bool get_GetHydroCellSignalConnected() {
    return GetHydroCellSignalConnected;
  }

  /**
   * @brief Set the signal connection flag.
   * @param connected Connection status
   */
  void set_GetHydroCellSignalConnected(bool m_GetHydroCellSignalConnected) {
    GetHydroCellSignalConnected = m_GetHydroCellSignalConnected;
  }

  /**
   * @brief Get number of droplets in the internal list.
   * @return Number of droplets
   */
  int get_dropletlist_size() const { return (dropletlist.size()); }

  /**
   * @brief Compute total energy from all droplets.
   * @return Total energy
   */
  Jetscape::real get_dropletlist_total_energy() const;

  /**
   * @brief Apply smearing kernel to a single droplet.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @param drop_i Droplet to smear
   * @param jmu Output source current 4-vector
   */
  virtual void smearing_kernel(Jetscape::real tau, Jetscape::real x,
                               Jetscape::real y, Jetscape::real eta,
                               const Droplet drop_i,
                               std::array<Jetscape::real, 4> &jmu) const {
    jmu = {0, 0, 0, 0};
  }

  /**
   * @brief Accumulate source term at a given space-time point.
   * @param tau Proper time
   * @param x Transverse x
   * @param y Transverse y
   * @param eta Space-time rapidity
   * @param jmu Output source current 4-vector
   */
  void get_source(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                  Jetscape::real eta, std::array<Jetscape::real, 4> &jmu) const;

  /**
   * @brief Clear all droplets from internal list.
   */
  virtual void Clear();
};

};  // namespace Jetscape

#endif  // LIQUEFIERBASE_H
