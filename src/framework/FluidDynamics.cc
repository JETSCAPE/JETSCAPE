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
// This is a general basic class for hydrodynamics

#include <iostream>
#include <array>
#include "FluidDynamics.h"
#include "LinearInterpolation.h"
#include "JetScapeSignalManager.h"

#include "MakeUniqueHelper.h"
#include "SurfaceFinder.h"


#define MAGENTA "\033[35m"

using namespace std;

namespace Jetscape {

FluidDynamics::FluidDynamics(){
    VERBOSE(8);
    eta = -99.99;
    SetId("FluidDynamics");
}

FluidDynamics::FluidDynamics(string m_name) : JetScapeModuleBase (m_name) {
    VERBOSE(8);
    eta = -99.99;
    SetId("FluidDynamics");
}


FluidDynamics::~FluidDynamics() {
    VERBOSE(8);
    disconnect_all();
}


void FluidDynamics::Init() {
    JetScapeModuleBase::Init();
    JSINFO << "Intialize FluidDynamics : "<< GetId() << " ...";
    fd = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro");
  
    if (!fd) {
        JSWARN << "Not a valid JetScape XML Hydro section file or "
               << "no XML file loaded!";
        exit(-1);
    }
  
    VERBOSE(8);  
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
        JSWARN << "No initialization module, "
               << "try: auto trento = make_shared<TrentoInitial>(); "
               << "jetscape->Add(trento);";
    }

    pre_eq_ptr = JetScapeSignalManager::Instance()->GetPreEquilibriumPointer().lock();
    if (!pre_eq_ptr) {
        JSWARN << "No Pre-equilibrium module";
    }
  
    InitializeHydro(parameter_list);
    InitTask();

    JetScapeTask::InitTasks();
}

void FluidDynamics::Exec() {
    JSINFO << "Run Hydro : " << GetId() << " ...";
    VERBOSE(8) << "Current Event #" << GetCurrentEvent();

    if (ini) {
        VERBOSE(3) << "length of entropy density vector="
                   << ini->GetEntropyDensityDistribution().size();
    }

    EvolveHydro();  
    JetScapeTask::ExecuteTasks();
}

void FluidDynamics::Clear() {
    clear_up_evolution_data();
    if (!weak_ptr_is_uninitialized(liquefier_ptr)) {
        liquefier_ptr.lock()->Clear();
    }
}


void FluidDynamics::CollectHeader(weak_ptr<JetScapeWriter> w) {
    auto f = w.lock();
    if ( f ) {
        auto& header = f->GetHeader();
        header.SetEventPlaneAngle( GetEventPlaneAngle() );
    }
}

  
// this function returns the energy density [GeV] at a space time point
// (time, x, y, z)
Jetscape::real FluidDynamics::GetEnergyDensity(
                    Jetscape::real time, Jetscape::real x,
                    Jetscape::real y, Jetscape::real z) {
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real energy_density = fluid_cell_ptr->energy_density;
    return(energy_density);
}

// this function returns the entropy density [GeV] at a space time point
// (time, x, y, z)
Jetscape::real FluidDynamics::GetEntropyDensity(
                    Jetscape::real time, Jetscape::real x,
                    Jetscape::real y, Jetscape::real z) {
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real entropy_density = fluid_cell_ptr->entropy_density;
    return(entropy_density);
}


// this function returns the temperature [GeV] at a space time point
// (time, x, y, z)
Jetscape::real FluidDynamics::GetTemperature(
                    Jetscape::real time, Jetscape::real x,
                    Jetscape::real y, Jetscape::real z) {
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real temperature = fluid_cell_ptr->temperature;
    return(temperature);
}


// this function returns the QGP fraction at a space time point
// (time, x, y, z)
Jetscape::real FluidDynamics::GetQgpFraction(
                    Jetscape::real time, Jetscape::real x,
                    Jetscape::real y, Jetscape::real z) {
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real qgp_fraction = fluid_cell_ptr->qgp_fraction;
    return(qgp_fraction);
}
    

void FluidDynamics::get_source_term(Jetscape::real tau, Jetscape::real x,
                                    Jetscape::real y, Jetscape::real eta,
                                    std::array<Jetscape::real, 4> jmu) const {
    liquefier_ptr.lock()->get_source(tau, x, y, eta, jmu);
}

  
void FluidDynamics::PrintFluidCellInformation(
						   FluidCellInfo* fluid_cell_info_ptr) {
    // this function print out the information of the fluid cell to the screen
    JSINFO << "=======================================================";
    JSINFO << "print out cell information:";
    JSINFO << "=======================================================";
    JSINFO << "energy density = " << fluid_cell_info_ptr->energy_density
           << " GeV/fm^3.";
    JSINFO << "entropy density = " << fluid_cell_info_ptr->entropy_density
           << " 1/fm^3.";
    JSINFO << "temperature = " << fluid_cell_info_ptr->temperature << " GeV.";
    JSINFO << "pressure = " << fluid_cell_info_ptr->pressure << " GeV/fm^3.";
    JSINFO << "QGP_fraction = " << fluid_cell_info_ptr->qgp_fraction;
    JSINFO << "mu_B = " << fluid_cell_info_ptr->mu_B << " GeV.";
    JSINFO << "mu_S = " << fluid_cell_info_ptr->mu_S << " GeV.";
    JSINFO << "mu_C = " << fluid_cell_info_ptr->mu_C << " GeV.";
    JSINFO << "vx = " << fluid_cell_info_ptr->vx;
    JSINFO << "vy = " << fluid_cell_info_ptr->vy;
    JSINFO << "vz = " << fluid_cell_info_ptr->vz;
    JSINFO << "shear viscous pi^{munu} (GeV/fm^3): ";
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            JSINFO << fluid_cell_info_ptr->pi[i][j];
        }
    }
    JSINFO << "bulk_Pi = " << fluid_cell_info_ptr->bulk_Pi << " GeV/fm^3";
    JSINFO << "=======================================================";
}


void FluidDynamics::UpdateEnergyDeposit(int t, double edop) {
    //sigslot::lock_block<multi_threaded_local> lock(this);
    JSDEBUG << MAGENTA << "Jet Signal received : "<< t << " " << edop;
}


} // end namespace Jetscape
