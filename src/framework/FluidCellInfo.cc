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

#include "FluidCellInfo.h"

namespace Jetscape {

FluidCellInfo::FluidCellInfo() {
  energy_density = 0.0;
  entropy_density = 0.0;
  temperature = 0.0;
  pressure = 0.0;
  qgp_fraction = 0.0;
  mu_B = 0.0;
  mu_S = 0.0;
  mu_C = 0.0;

  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      pi[i][j] = 0.0;
    }
  }
  bulk_Pi = 0.0;
}

}  // namespace Jetscape
