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

#ifndef MUSICWRAPPER_H
#define MUSICWRAPPER_H

#include <memory>

#include "FluidDynamics.h"
#include "music.h"
#include "hydro_source_base.h"
#include "LiquefierBase.h"
#include "data_struct.h"

using namespace Jetscape;


class HydroSourceJETSCAPE : public HydroSourceBase {
 private:
    std::weak_ptr<LiquefierBase> liquefier_ptr;

 public:
    HydroSourceJETSCAPE() = default;
    ~HydroSourceJETSCAPE() {}

    void add_a_liqueifier(std::shared_ptr<LiquefierBase> new_liqueifier) {
        liquefier_ptr = new_liqueifier;
    }

    int get_number_of_sources() const {
        return(liquefier_ptr.lock()->get_dropletlist_size());
    }
    
    //! this function returns the energy source term J^\mu at a given point
    //! (tau, x, y, eta_s)
    void get_hydro_energy_source(
        const double tau, const double x, const double y, const double eta_s,
        const FlowVec &u_mu, EnergyFlowVec &j_mu) const {
        std::array<Jetscape::real, 4> jmu_tmp = {0.0};
        liquefier_ptr.lock()->get_source(tau, x, y, eta_s, jmu_tmp);
        for (int i = 0; i < 4; i++)
            j_mu[i] = jmu_tmp[i];
    }

};

//! this is wrapper class for MUSIC so that it can be used as a external
//! library for the JETSCAPE integrated framework
class MpiMusic: public FluidDynamics {
 private:
    // int mode;            //!< records running mode
    MUSIC *music_hydro_ptr;
    int doCooperFrye;    //!< flag to run Cooper-Frye freeze-out
                         //!< for soft particles
    int flag_output_evo_to_file;
    std::shared_ptr<HydroSourceJETSCAPE> hydro_source_terms_ptr;

 public:
     MpiMusic();
     ~MpiMusic();

     void InitializeHydro(Parameter parameter_list);

     void EvolveHydro();
     void GetHydroInfo(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
		std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);
     
     void GetHydroInfo_JETSCAPE(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
		std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);
     void GetHydroInfo_MUSIC(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
		std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr);

     void SetHydroGridInfo();
     void PassHydroEvolutionHistoryToFramework();
    
     void add_a_liqueifier(std::shared_ptr<LiquefierBase> new_liqueifier) {
        liquefier_ptr = new_liqueifier;
        hydro_source_terms_ptr->add_a_liqueifier(liquefier_ptr.lock());
    }

     void GetHyperSurface(Jetscape::real T_cut,
                          SurfaceCellInfo* surface_list_ptr) {};
     void collect_freeze_out_surface();

};



#endif // MUSICWRAPPER_H
