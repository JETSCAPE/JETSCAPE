/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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

#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "SurfaceCellInfo.h"

namespace Jetscape {

class SoftParticlization : public JetScapeModuleBase {
private:
  bool HydroHyperSurfaceConnected_;
  bool ClearHydroHyperSurfaceConnected_;

public:
  SoftParticlization();
  ~SoftParticlization();

  virtual void Init();
  virtual void Exec();
  virtual void Clear();

  sigslot::signal1<std::vector<SurfaceCellInfo> &,
                   multi_threaded_local> GetHydroHyperSurface;
  sigslot::signal0<multi_threaded_local> ClearHydroHyperSurface;

  void SetGetHydroHyperSurfaceConnected(bool m_GetHydroHyperSurfaceConnected) {
    HydroHyperSurfaceConnected_ = m_GetHydroHyperSurfaceConnected;
  }

  void SetClearHydroHyperSurfaceConnected(bool m_ClearHydroHyperSurfaceConnected) {
    ClearHydroHyperSurfaceConnected_ = m_ClearHydroHyperSurfaceConnected;
  }

  bool GetGetHydroHyperSurfaceConnected() const {
    return HydroHyperSurfaceConnected_;
  }

  bool GetClearHydroHyperSurfaceConnected() const {
    return ClearHydroHyperSurfaceConnected_;
  }

  std::vector<std::vector<shared_ptr<Hadron>>> Hadron_list_;

  bool boost_invariance;
  bool check_boost_invariance();
};

} // end namespace Jetscape

#endif // SOFTPARTICLIZATION_H_
