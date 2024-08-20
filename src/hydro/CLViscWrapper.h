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

#ifndef CLVISC_WRAPPER_H
#define CLVISC_WRAPPER_H

#include "FluidDynamics.h"
#include "clvisc.h"

using namespace Jetscape;

// CLViscWrapper for jetscape
class CLVisc : public FluidDynamics {
 private:
  std::unique_ptr<clvisc::CLVisc> hydro_;
  int doCooperFrye;
  // scale the initial density because Trento only provides
  // the density distribution profile without scaling factor
  double initial_condition_scale_factor;

  // Allows the registration of the module
  static RegisterJetScapeModule<CLVisc> reg;

 public:
  CLVisc();
  ~CLVisc();

  void InitializeHydro(Parameter parameter_list);
  void EvolveHydro();
  void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y,
                    Jetscape::real z,
                    std::unique_ptr<FluidCellInfo> &fluid_cell_info_ptr);
  void GetHyperSurface(Jetscape::real T_cut,
                       SurfaceCellInfo *surface_list_ptr){};
};

#endif  // CLViscWRAPPER_H
