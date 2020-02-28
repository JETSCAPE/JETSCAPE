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

#ifndef GUBSERHYDRO_H
#define GUBSERHYDRO_H

#include "JetScapeLogger.h"

#include "FluidDynamics.h"

using namespace Jetscape;

class GubserHydro : public FluidDynamics {
  // this is wrapper class for a simple brick
  // so that it can be used within the JETSCAPE framework
private:
  double q;
  double e_0;
  double temperature(double e_local);

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<GubserHydro> reg;

public:
  GubserHydro();
  ~GubserHydro();

  void InitializeHydro(Parameter parameter_list);

  void EvolveHydro();

  void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                    Jetscape::real z,
                    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);

  void GetHyperSurface(Jetscape::real T_cut,
                       SurfaceCellInfo *surface_list_ptr){};
};

#endif // GUBSERHYDRO_H
