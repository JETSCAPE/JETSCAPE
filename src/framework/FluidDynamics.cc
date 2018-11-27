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
#include "./FluidDynamics.h"
#include "./LinearInterpolation.h"
#include "JetScapeSignalManager.h"

#define MAGENTA "\033[35m"

using namespace std;

namespace Jetscape {
  /** For one given time step id_tau,
   * get FluidCellInfo at spatial point (x, y, eta)*/
  FluidCellInfo EvolutionHistory::GetAtTimeStep(int id_tau,
						   real x, real y, real eta) {
    int id_x = GetIdX(x);
    int id_y = GetIdY(y);
    int id_eta = GetIdEta(eta);
    // cijk for idx=i, idy=j and id_eta=k
    int c000 = CellIndex(id_tau, id_x, id_y, id_eta);
    int c001 = CellIndex(id_tau, id_x, id_y, id_eta+1);
    int c010 = CellIndex(id_tau, id_x, id_y+1, id_eta);
    int c011 = CellIndex(id_tau, id_x, id_y+1, id_eta+1);
    int c100 = CellIndex(id_tau, id_x+1, id_y, id_eta);
    int c101 = CellIndex(id_tau, id_x+1, id_y, id_eta+1);
    int c110 = CellIndex(id_tau, id_x+1, id_y+1, id_eta);
    int c111 = CellIndex(id_tau, id_x+1, id_y+1, id_eta+1);
    real x0 = XCoord(id_x);
    real x1 = XCoord(id_x + 1);
    real y0 = YCoord(id_y);
    real y1 = YCoord(id_y + 1);
    real eta0 = EtaCoord(id_eta);
    real eta1 = EtaCoord(id_eta + 1);

    return TrilinearInt(x0, x1, y0, y1, eta0, eta1,
			 data.at(c000), data.at(c001), data.at(c010), data.at(c011),
			 data.at(c100), data.at(c101), data.at(c110), data.at(c111),
			 x, y, eta);
  }

  // do interpolation along time direction; we may also need high order
  // interpolation functions 
  FluidCellInfo EvolutionHistory::get(real tau, real x, real y, real eta){
    CheckInRange(tau, x, y, eta);
    int id_tau = GetIdTau(tau);
    real tau0 = TauCoord(id_tau);
    real tau1 = TauCoord(id_tau + 1);
    FluidCellInfo bulk0 = GetAtTimeStep(id_tau, x, y, eta);
    FluidCellInfo bulk1 = GetAtTimeStep(id_tau+1, x, y, eta);
    return LinearInt(tau0, tau1, bulk0, bulk1, tau);
  }

  FluidDynamics::FluidDynamics(){
    VERBOSE(8);
    eta=-99.99;
    SetId("FluidDynamics");
  }

  FluidDynamics::FluidDynamics(string m_name) : JetScapeModuleBase (m_name) {
    VERBOSE(8);
    eta=-99.99;
    SetId("FluidDynamics");
  }


  FluidDynamics::~FluidDynamics()
  {
    VERBOSE(8);
    disconnect_all();
  }

  void FluidDynamics::Init()
  {
    JetScapeModuleBase::Init();
    JSINFO<<"Intialize FluidDynamics : "<<GetId()<< " ...";
    fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" );
  
    if (!fd) {
      JSWARN << "Not a valid JetScape XML Hydro section file or no XML file loaded!";
      exit(-1);
    }
  
    VERBOSE(8);  
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
      JSWARN << "No initialization module, try: auto trento = make_shared<TrentoInitial>(); jetscape->Add(trento);";
    }

    pre_eq_ptr = JetScapeSignalManager::Instance()->GetPreEquilibriumPointer().lock();
    if (!pre_eq_ptr) {
      JSWARN << "No Pre-equilibrium module";
    }
  
    InitializeHydro(parameter_list);
    InitTask();

    JetScapeTask::InitTasks();
  }

  void FluidDynamics::Exec()
  {
    JSINFO <<"Run Hydro : "<<GetId()<< " ...";
    VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

    if (ini) {
      VERBOSE(3) << "length of entropy density vector=" << ini->GetEntropyDensityDistribution().size();
    }

    EvolveHydro();  
    JetScapeTask::ExecuteTasks();
  }

  void FluidDynamics::UpdateEnergyDeposit(int t, double edop) {
    //sigslot::lock_block<multi_threaded_local> lock(this);
    JSDEBUG<<MAGENTA<<"Jet Signal received : "<<t<<" "<<edop;
  }
  
  void FluidDynamics::GetEnergyDensity(int t,double &edensity) {
    //sigslot::lock_block<multi_threaded_local> lock(this);
    edensity=0.5;
    JSDEBUG<<"Edensity to Jet = "<<edensity<<" at t="<<t;
  }

  void FluidDynamics::CollectHeader( weak_ptr<JetScapeWriter> w ){
    auto f = w.lock();
    if ( f ){
      auto& header = f->GetHeader();
      header.SetEventPlaneAngle( GetEventPlaneAngle() );
    }
  }

  
  real FluidDynamics::GetEnergyDensity(real time, real x, real y, real z) {
    // this function returns the energy density [GeV] at a space time point
    // (time, x, y, z)
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real energy_density = fluid_cell_ptr->energy_density;
    // delete fluid_cell_ptr;
    return(energy_density);
  }

  real FluidDynamics::GetEntropyDensity(real time, real x, real y, real z) {
    // this function returns the entropy density [GeV] at a space time point
    // (time, x, y, z)
    // FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real entropy_density = fluid_cell_ptr->entropy_density;
    //delete fluid_cell_ptr;
    return(entropy_density);
  }

  real FluidDynamics::GetTemperature(real time, real x, real y, real z) {
    // this function returns the temperature [GeV] at a space time point
    // (time, x, y, z)
    // FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real temperature = fluid_cell_ptr->temperature;
    // delete fluid_cell_ptr;
    return(temperature);
  }

  real FluidDynamics::GetQgpFraction(real time, real x, real y, real z) {
    // this function returns the QGP fraction at a space time point
    // (time, x, y, z)
    // FluidCellInfo *fluid_cell_ptr = new FluidCellInfo;
    std::unique_ptr<FluidCellInfo> fluid_cell_ptr;
    GetHydroInfo(time, x, y, z, fluid_cell_ptr);
    real qgp_fraction = fluid_cell_ptr->qgp_fraction;
    // delete fluid_cell_ptr;
    return(qgp_fraction);
  }

  void FluidDynamics::PrintFluidCellInformation(
						   FluidCellInfo* fluid_cell_info_ptr) {
    // this function print out the information of the fluid cell to the screen
    cout << "=======================================================" << endl;
    cout << "print out cell information:" << endl;
    cout << "=======================================================" << endl;
    cout << "energy density = " << fluid_cell_info_ptr->energy_density
         << " GeV/fm^3." << endl;
    cout << "entropy density = " << fluid_cell_info_ptr->entropy_density
         << " 1/fm^3." << endl;
    cout << "temperature = " << fluid_cell_info_ptr->temperature
         << " GeV." << endl;
    cout << "pressure = " << fluid_cell_info_ptr->pressure
         << " GeV/fm^3." << endl;
    cout << "QGP_fraction = " << fluid_cell_info_ptr->qgp_fraction
         << endl;
    cout << "mu_B = " << fluid_cell_info_ptr->mu_B << " GeV." << endl;
    cout << "mu_S = " << fluid_cell_info_ptr->mu_S << " GeV." << endl;
    cout << "mu_C = " << fluid_cell_info_ptr->mu_C << " GeV." << endl;
    cout << "vx = " << fluid_cell_info_ptr->vx << endl;
    cout << "vy = " << fluid_cell_info_ptr->vy << endl;
    cout << "vz = " << fluid_cell_info_ptr->vz << endl;
    cout << "shear viscous pi^{munu} (GeV/fm^3): " << endl;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	cout << fluid_cell_info_ptr->pi[i][j] << "   ";
      }
      cout << endl;
    }
    cout << "bulk_Pi = " << fluid_cell_info_ptr->bulk_Pi << " GeV/fm^3"
         << endl;
    cout << "=======================================================" << endl;
  }

  void FluidCellInfo::Print()
  {
    // this function print out the information of the fluid cell to the screen
    cout << "=======================================================" << endl;
    cout << "print out cell information:" << endl;
    cout << "=======================================================" << endl;
    cout << "energy density = " << energy_density
         << " GeV/fm^3." << endl;
    cout << "entropy density = " << entropy_density
         << " 1/fm^3." << endl;
    cout << "temperature = " << temperature
         << " GeV." << endl;
    cout << "pressure = " << pressure
         << " GeV/fm^3." << endl;
    cout << "QGP_fraction = " << qgp_fraction
         << endl;
    cout << "mu_B = " << mu_B << " GeV." << endl;
    cout << "mu_S = " << mu_S << " GeV." << endl;
    cout << "mu_C = " << mu_C << " GeV." << endl;
    cout << "vx = " << vx << endl;
    cout << "vy = " << vy << endl;
    cout << "vz = " << vz << endl;
    cout << "shear viscous pi^{munu} (GeV/fm^3): " << endl;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	cout << pi[i][j] << "   ";
      }
      cout << endl;
    }
    cout << "bulk_Pi = " << bulk_Pi << " GeV/fm^3"
         << endl;
    cout << "=======================================================" << endl;
  }

} // end namespace Jetscape
