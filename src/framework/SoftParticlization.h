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
// JETSCAPE module for soft particlization
// This module will generate Monte-Carlo samples for soft hadrons
// -----------------------------------------

#ifndef SOFTPARTICLIZATION_H_
#define SOFTPARTICLIZATION_H_

#include <vector>

#include "JetClass.h"
#include "JetScapeModuleBase.h"
#include "JetScapeWriter.h"
#include "SurfaceCellInfo.h"

namespace Jetscape {

/**
 * @brief JETSCAPE module for soft particlization
 * 
 * This module will generate Monte-Carlo samples for soft hadrons from the
 * hydrodynamic output.
 */
class SoftParticlization : public JetScapeModuleBase {
 private:
  /// Flag for the connection status of the GetHydroHyperSurface signal
  bool HydroHyperSurfaceConnected_;

  /// Flag for the connection status of the ClearHydroHyperSurface signal
  bool ClearHydroHyperSurfaceConnected_;

 public:
  /**
   * @brief Construct a new SoftParticlization object
   */
  SoftParticlization();

  /**
   * @brief Destroy the SoftParticlization object
   */
  ~SoftParticlization();

  /**
   * @brief Initialize the SoftParticlization module
   */
  virtual void Init();

  /**
   * @brief Execute the SoftParticlization module
   */
  virtual void Exec();

  /**
   * @brief Clear the SoftParticlization module
   */
  virtual void Clear();

  /**
   * @brief Signal for getting the hydrodynamic hypersurface
   */
  sigslot::signal1<std::vector<SurfaceCellInfo> &, multi_threaded_local>
      GetHydroHyperSurface;
  
  /**
   * @brief Signal for clearing the hydrodynamic hypersurface
   */
  sigslot::signal0<multi_threaded_local> ClearHydroHyperSurface;

  /**
   * @brief Set the GetHydroHyperSurfaceConnected flag
   * 
   * @param m_GetHydroHyperSurfaceConnected Boolean flag
   */
  void SetGetHydroHyperSurfaceConnected(bool m_GetHydroHyperSurfaceConnected) {
    HydroHyperSurfaceConnected_ = m_GetHydroHyperSurfaceConnected;
  }

  /**
   * @brief Set the ClearHydroHyperSurfaceConnected flag
   * 
   * @param m_ClearHydroHyperSurfaceConnected Boolean flag
   */
  void SetClearHydroHyperSurfaceConnected(
      bool m_ClearHydroHyperSurfaceConnected) {
    ClearHydroHyperSurfaceConnected_ = m_ClearHydroHyperSurfaceConnected;
  }

  /**
   * @brief Get the GetHydroHyperSurfaceConnected flag
   * 
   * @return Boolean
   */
  bool GetGetHydroHyperSurfaceConnected() const {
    return HydroHyperSurfaceConnected_;
  }

  /**
   * @brief Get the ClearHydroHyperSurfaceConnected flag
   * 
   * @return Boolean
   */
  bool GetClearHydroHyperSurfaceConnected() const {
    return ClearHydroHyperSurfaceConnected_;
  }

  /// List of hadrons
  std::vector<std::vector<shared_ptr<Hadron>>> Hadron_list_;

  /// Flag for boost invariance
  bool boost_invariance;

  /**
   * @brief Check the boost invariance
   * 
   * @return Boolean
   */
  bool check_boost_invariance();
};

}  // end namespace Jetscape

#endif  // SOFTPARTICLIZATION_H_
