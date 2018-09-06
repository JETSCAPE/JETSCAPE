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
  // convert the string type entry name to enum type EntryNames
  EntryName ResolveEntryName(std::string input) {
      static const std::map<std::string, EntryName> optionStrings = {
            { "energy_density", ENTRY_ENERGY_DENSITY},
            { "entropy_density", ENTRY_ENTROPY_DENSITY},
            { "temperature", ENTRY_TEMPERATURE},
            { "pressure", ENTRY_PRESSURE},
            { "qgp_fraction", ENTRY_QGP_FRACTION},
            { "mu_b", ENTRY_MU_B},
            { "mu_c", ENTRY_MU_C},
            { "mu_s", ENTRY_MU_S},
            { "vx", ENTRY_VX},
            { "vy", ENTRY_VY},
            { "vz", ENTRY_VZ},
            { "pi00", ENTRY_PI00},
            { "pi01", ENTRY_PI01},
            { "pi02", ENTRY_PI02},
            { "pi03", ENTRY_PI03},
            { "pi11", ENTRY_PI11},
            { "pi12", ENTRY_PI12},
            { "pi13", ENTRY_PI13},
            { "pi22", ENTRY_PI22},
            { "pi23", ENTRY_PI23},
            { "pi33", ENTRY_PI33},
            { "bulk_pi", ENTRY_BULK_PI},
        };
      auto itr = optionStrings.find(input);
      if (itr != optionStrings.end()) {
          return itr->second;
      } else {
          return ENTRY_INVALID;
      }
 }

  /* This function will read the sparse data stored in data_ with associated 
   * information data_info_ into to FluidCellInfo object */
  FluidCellInfo EvolutionHistory::GetFluidCell(int id_tau,
          int id_x, int id_y, int id_eta) {
      auto fluid_cell_ptr = std::make_unique<FluidCellInfo>();
      int entries_per_record = data_info.size();
      int id_eta_corrected = id_eta;
      // set id_eta=0 if hydro is in 2+1D mode
      if (neta == 0 || neta == 1) { id_eta_corrected = 0; }
      int record_starting_id = CellIndex(id_tau, id_x, id_y, id_eta_corrected)
                               *entries_per_record;
      for (int i=0; i<entries_per_record; i++) {
          auto entry_name = ResolveEntryName(data_info.at(i));
          auto entry_data = data.at(record_starting_id + i);
          switch (entry_name) {
              case ENTRY_ENERGY_DENSITY:
                  fluid_cell_ptr->energy_density = entry_data;
                  break;
              case ENTRY_ENTROPY_DENSITY:
                  fluid_cell_ptr->entropy_density = entry_data;
                  break;
              case ENTRY_TEMPERATURE:
                  fluid_cell_ptr->temperature = entry_data;
                  break;
              case ENTRY_PRESSURE:
                  fluid_cell_ptr->pressure = entry_data;
                  break;
              case ENTRY_QGP_FRACTION:
                  fluid_cell_ptr->qgp_fraction = entry_data;
                  break;
              case ENTRY_MU_B:
                  fluid_cell_ptr->mu_B = entry_data;
                  break;
              case ENTRY_MU_C:
                  fluid_cell_ptr->mu_C = entry_data;
                  break;
              case ENTRY_MU_S:
                  fluid_cell_ptr->mu_S = entry_data;
                  break;
              case ENTRY_VX:
                  fluid_cell_ptr->vx = entry_data;
                  break;
              case ENTRY_VY:
                  fluid_cell_ptr->vy = entry_data;
                  break;
              case ENTRY_VZ:
                  fluid_cell_ptr->vz = entry_data;
                  break;
              case ENTRY_PI00:
                  fluid_cell_ptr->pi[0][0] = entry_data;
                  break;
              case ENTRY_PI01:
                  fluid_cell_ptr->pi[0][1] = entry_data;
                  fluid_cell_ptr->pi[1][0] = entry_data;
                  break;
              case ENTRY_PI02:
                  fluid_cell_ptr->pi[0][2] = entry_data;
                  fluid_cell_ptr->pi[2][0] = entry_data;
                  break;
              case ENTRY_PI03:
                  fluid_cell_ptr->pi[0][3] = entry_data;
                  fluid_cell_ptr->pi[3][0] = entry_data;
                  break;
              case ENTRY_PI11:
                  fluid_cell_ptr->pi[1][1] = entry_data;
                  break;
              case ENTRY_PI12:
                  fluid_cell_ptr->pi[1][2] = entry_data;
                  fluid_cell_ptr->pi[2][1] = entry_data;
                  break;
              case ENTRY_PI13:
                  fluid_cell_ptr->pi[1][3] = entry_data;
                  fluid_cell_ptr->pi[3][1] = entry_data;
                  break;
              case ENTRY_PI22:
                  fluid_cell_ptr->pi[2][2] = entry_data;
                  break;
              case ENTRY_PI23:
                  fluid_cell_ptr->pi[2][3] = entry_data;
                  fluid_cell_ptr->pi[3][2] = entry_data;
                  break;
              case ENTRY_PI33:
                  fluid_cell_ptr->pi[3][3] = entry_data;
                  break;
              case ENTRY_BULK_PI:
                  fluid_cell_ptr->bulk_Pi = entry_data;
                  break;
              default:
                  WARN << "The entry name in data_info_ must be one of the \
                          energy_density, entropy_density, temperature, pressure, qgp_fraction, \
                          mu_b, mu_c, mu_s, vx, vy, vz, pi00, pi01, pi02, pi03, pi11, pi12, \
                          pi13, pi22, pi23, pi33, bulk_pi";
                  break;
          }
      }

      return *fluid_cell_ptr;
  }


  /** Construct evolution history given the bulk_data and the data_info */
  void EvolutionHistory::Construct(const std::vector<float> & data_,
                   const std::vector<std::string> & data_info_,
                   float tau_min_, float dtau_, 
                   float x_min_, float dx_, int nx_,
                   float y_min_, float dy_, int ny_,
                   float eta_min_, float deta_, int neta_,
                   bool tau_eta_is_tz_) {
      data = data_;
      data_info = data_info_;
      tau_min = tau_min_;
      x_min = x_min_;
      y_min = y_min_;
      eta_min = eta_min_;
      dtau = dtau_;
      dx = dx_;
      dy = dy_;
      deta = deta_;
      nx = nx_;
      ny = ny_;
      neta = neta_;
      tau_eta_is_tz = tau_eta_is_tz_;
      ntau = data_.size() / (data_info_.size() * nx * ny * neta);
  }


  /** For one given time step id_tau,
   * get FluidCellInfo at spatial point (x, y, eta)*/
  FluidCellInfo EvolutionHistory::GetAtTimeStep(int id_tau,
						   real x, real y, real eta) {
    int id_x = GetIdX(x);
    int id_y = GetIdY(y);
    int id_eta = GetIdEta(eta);
    auto c000 = GetFluidCell(id_tau, id_x, id_y, id_eta);
    auto c001 = GetFluidCell(id_tau, id_x, id_y, id_eta+1);
    auto c010 = GetFluidCell(id_tau, id_x, id_y+1, id_eta);
    auto c011 = GetFluidCell(id_tau, id_x, id_y+1, id_eta+1);
    auto c100 = GetFluidCell(id_tau, id_x+1, id_y, id_eta);
    auto c101 = GetFluidCell(id_tau, id_x+1, id_y, id_eta+1);
    auto c110 = GetFluidCell(id_tau, id_x+1, id_y+1, id_eta);
    auto c111 = GetFluidCell(id_tau, id_x+1, id_y+1, id_eta+1);
    real x0 = XCoord(id_x);
    real x1 = XCoord(id_x + 1);
    real y0 = YCoord(id_y);
    real y1 = YCoord(id_y + 1);
    real eta0 = EtaCoord(id_eta);
    real eta1 = EtaCoord(id_eta + 1);

    return TrilinearInt(x0, x1, y0, y1, eta0, eta1,
			 c000, c001, c010, c011, c100, c101, c110, c111,
			 x, y, eta);
  }

  // do interpolation along time direction; we may also need high order
  // interpolation functions 
  FluidCellInfo EvolutionHistory::Get(real tau, real x, real y, real eta){
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
    INFO<<"Intialize FluidDynamics : "<<GetId()<< " ...";
    fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Hydro" );
  
    if (!fd) {
      WARN << "Not a valid JetScape XML Hydro section file or no XML file loaded!";
      exit(-1);
    }
  
    VERBOSE(8);  
    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
      WARN << "No initialization module, try: auto trento = make_shared<TrentoInitial>(); jetscape->Add(trento);";
    }

    pre_eq_ptr = JetScapeSignalManager::Instance()->GetPreEquilibriumPointer().lock();
    if (!pre_eq_ptr) {
      WARN << "No Pre-equilibrium module";
    }
  
    InitializeHydro(parameter_list);
    InitTask();

    JetScapeTask::InitTasks();
  }

  void FluidDynamics::Exec()
  {
    INFO <<"Run Hydro : "<<GetId()<< " ...";
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
                const std::unique_ptr<FluidCellInfo> & fluid_cell_info_ptr) {
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
