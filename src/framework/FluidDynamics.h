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

#ifndef FLUIDDYNAMICS_H
#define FLUIDDYNAMICS_H

#include <memory>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <map>

#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "PreequilibriumDynamics.h"
#include "RealType.h"

namespace Jetscape {

  /// Flags for hydrodynamics status.  
  enum HydroStatus {NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR};

  /// placeholder for the future
  class JetSource {
  public:
    JetSource() {}
    // JetSource():j0(0.), j1(0.), j2(0.), j3(0.) {}
  private:
    // Jetscape::real j0, j1, j2, j3;
  };

  enum EntryName {ENTRY_ENERGY_DENSITY, ENTRY_ENTROPY_DENSITY, ENTRY_TEMPERATURE,
      ENTRY_PRESSURE, ENTRY_QGP_FRACTION, ENTRY_MU_B, ENTRY_MU_C, ENTRY_MU_S,
      ENTRY_VX, ENTRY_VY, ENTRY_VZ,
      ENTRY_PI00, ENTRY_PI01, ENTRY_PI02, ENTRY_PI03,
      ENTRY_PI11, ENTRY_PI12, ENTRY_PI13, 
      ENTRY_PI22, ENTRY_PI23, ENTRY_PI33, ENTRY_BULK_PI, ENTRY_INVALID
  };

  EntryName ResolveEntryName(std::string input);

  //overload +-*/ for easier linear interpolation
  class FluidCellInfo {
  public:
    // data structure for outputing fluid cell information
    Jetscape::real energy_density;    //!< Local energy density [GeV/fm^3].
    Jetscape::real entropy_density;   //!< Local entropy density [1/fm^3].
    Jetscape::real temperature;       //!< Local temperature [GeV].
    Jetscape::real pressure;          //!< Thermal pressure [GeV/fm^3].
    Jetscape::real qgp_fraction;      //!< Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
    Jetscape::real mu_B;              //!< Net baryon chemical potential [GeV].
    Jetscape::real mu_C;              //!< Net charge chemical potential [GeV]
    Jetscape::real mu_S;              //!< Net strangeness chemical potential [GeV].
    Jetscape::real vx, vy, vz;        //!< Flow velocity.
    Jetscape::real pi[4][4];          //!< Shear stress tensor [GeV/fm^3].
    Jetscape::real bulk_Pi;           //!< Bulk viscous pressure [GeV/fm^3].
    /** Default constructor.*/
    FluidCellInfo() = default;    

    /** @param b Multiply the fluid cell by scalar factor b. */ 
    FluidCellInfo inline operator*=(Jetscape::real b);

    /** Prints fluid cell properties to the screen. */
    void Print();
  };

  /// adds \f$ c = a + b \f$
  inline FluidCellInfo operator+(FluidCellInfo a, const FluidCellInfo & b) {
    a.energy_density += b.energy_density;
    a.entropy_density += b.entropy_density;
    a.temperature += b.temperature;
    a.pressure += b.pressure;
    a.qgp_fraction += b.qgp_fraction;
    a.mu_B += b.mu_B;
    a.mu_C += b.mu_C;
    a.mu_S += b.mu_S;
    a.vx += b.vx;
    a.vy += b.vy;
    a.vz += b.vz;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	a.pi[i][j] += b.pi[i][j];
      }
    }
    a.bulk_Pi += b.bulk_Pi;
    return a;
  }

  // Multiply the fluid cell with a scalar factor
  FluidCellInfo inline FluidCellInfo::operator*=(Jetscape::real b){
    this->energy_density *= b;
    this->entropy_density *= b;
    this->temperature *= b;
    this->pressure *= b;
    this->qgp_fraction *= b;
    this->mu_B *= b;
    this->mu_C *= b;
    this->mu_S *= b;
    this->vx *= b;
    this->vy *= b;
    this->vz *= b;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	this->pi[i][j] *= b;
      }
    }
    this->bulk_Pi *= b;
    return *this;
  }

  /// multiply \f$ c = a * b \f$
  inline FluidCellInfo operator*(Jetscape::real a, FluidCellInfo b){
    b *= a;
    return b;
  }

  /// multiply \f$ c = a * b \f$
  inline FluidCellInfo operator*(FluidCellInfo a, Jetscape::real b){
    a *= b;
    return a;
  }

  /// division \f$ c = a / b \f$
  inline FluidCellInfo operator/(FluidCellInfo a, Jetscape::real b){
    a *= 1.0/b;
    return a;
  }

  // print the fluid cell information for debuging
  // this function has bugs
  //std::ostream &operator<<(std::ostream &os, const FluidCellInfo &cell) {
  //    os << "energy_density=" << cell.energy_density << std::endl; 
  //    os << "entropy_density=" << cell.entropy_density << std::endl; 
  //    os << "temperature=" << cell.temperature << std::endl; 
  //    os << "pressure=" << cell.pressure << std::endl;
  //    os << "qgp_fraction=" << cell.qgp_fraction << std::endl;
  //    os << "mu_B=" << cell.mu_B << std::endl;
  //    os << "mu_C=" << cell.mu_C << std::endl;
  //    os << "mu_S=" << cell.mu_S << std::endl;
  //    os << "vx=" << cell.vx << std::endl;
  //    os << "vy=" << cell.vy << std::endl;
  //    os << "vz=" << cell.vz << std::endl;
  //    os << "pi[mu][nu]=" << std::endl;
  //    for (int i = 0; i < 4; i++) {
  //        for (int j = 0; j < 4; j++) {
  //            os << cell.pi[i][j] << ' ';
  //        }
  //        os << std::endl;
  //    }
  //    os << "bulk_Pi=" << cell.bulk_Pi;
  //    return os << std::endl;
  //}

  class SurfaceCellInfo {
  public:
    // data structure for outputing hyper-surface information
    Jetscape::real d3sigma_mu[4];     //!< Surface vector.
    Jetscape::real energy_density;    //!< Local energy density [GeV/fm^3].
    Jetscape::real entropy_density;   //!< Local entropy density [1/fm^3].
    Jetscape::real temperature;       //!< Local temperature [GeV].
    Jetscape::real pressure;          //!< Thermal pressure [GeV/fm^3].
    Jetscape::real qgp_fraction;      //!< Fraction of quark gluon plasma assuming medium is in QGP+HRG phase.
    Jetscape::real mu_B;              //!< Net baryon chemical potential [GeV].
    Jetscape::real mu_C;              //!< Net charge chemical potential [GeV].
    Jetscape::real mu_S;              //!< Net strangeness chemical potential [GeV].
    Jetscape::real vx, vy, vz;        //!< Flow velocity.
    Jetscape::real pi[4][4];          //!< Shear stress tensor [GeV/fm^3].
    Jetscape::real bulk_Pi;           //!< Bulk viscous pressure [GeV/fm^3].

    /** Default constructor. */
    SurfaceCellInfo() {};

    /** Destructor. */
    ~SurfaceCellInfo() {};
  };


  class Parameter{
  public:
    /** Hydro dynamics parameters file name.*/
    char* hydro_input_filename;
  };

  class InvalidSpaceTimeRange : public std::invalid_argument {
    using std::invalid_argument::invalid_argument;
  };

  class EvolutionHistory{
  public:
    /** @param tau_min Minimum value of tau.*/
    Jetscape::real tau_min, dtau; //!< @param dtau Step-size for tau.
    /** @param x_min Minimum value of x. */
    Jetscape::real x_min, dx;     //!< @param dx Step-size for x. 
    /** @param y_min Minimum value of y. */
    Jetscape::real y_min, dy;     //!< @param dy Step-size for y.
    /** @param eta_min Minimum value of eta. */
    Jetscape::real eta_min, deta; //!< @param deta Step-size for eta. 
    int ntau;   //!< @param ntau Number of grid points in tau-axis.
    int nx;    //!< @param nx Number of grid points in x-axis. 
    int ny;   //!< @param ny Number of grid points in y-axis.
    int neta;  //!< @param neta Number of grid points in eta-axis.   

    /** Default is set to false. Set flag tau_eta_is_tz to true if hydro dynamics is setup in (t,x,y,z) coordinate. */
    bool tau_eta_is_tz;   

    /** The bulk information of hydro dynamics.
     * It is easiest to pass vector of float from 3rd party hydro module
     * to Jetscape Evolution History, where bulk info should be stored
     * in orders of ed0, sd0, temp0, ..., ed1, sd1, temp1, ..., edn, sdn, tempn.
     * The order and data entry is given by data_info */
    std::vector<float> data;

    /** Store the entry names of one record in the data array*/
    std::vector<std::string> data_info;

    /** Default constructor. */
    EvolutionHistory() {};

    void Construct(const std::vector<float> & data_,
                   const std::vector<std::string> & data_info_,
                   float tau_min, float dtau,
                   float x_min, float dx, int nx,
                   float y_min, float dy, int ny,
                   float eta_min, float deta, int neta,
                   bool tau_eta_is_tz);

    /** Default destructor. */
    ~EvolutionHistory() {}

    /** Maximum value of tau. */
    inline bool DataSizeMatch() {return data.size() == (ntau * nx * ny * neta * data_info.size());}

    /** Maximum value of tau. */
    inline Jetscape::real TauMax() {return tau_min + (ntau - 1) * dtau;}

    /** Maximum value of x. */
    inline Jetscape::real XMax() {return x_min + (nx - 1) * dx;}

    /** Maximum value of y. */
    inline Jetscape::real YMax() {return y_min + (ny - 1) * dy;}

    /** Maximum value of eta. */
    inline Jetscape::real EtaMax() {return eta_min + (neta - 1) * deta;}

    /** It checks whether a space-time point (tau, x, y, eta) is inside evolution history or outside.
	@param tau Light-cone coordinate.
	@param x  Space coordinate.   
	@param y  Space coordinate.
	@param eta Light-cone coordinate.
    */
    void CheckInRange(Jetscape::real tau, Jetscape::real x, Jetscape::real y,
                      Jetscape::real eta)
    {
      if (tau < tau_min || tau > TauMax()) {
        throw InvalidSpaceTimeRange("tau=" + std::to_string(tau)
				    + " is not in range [" + std::to_string(tau_min) + "," 
				    + std::to_string(TauMax()) + "]");
      }
      if (x < x_min || x > XMax()) {
        throw InvalidSpaceTimeRange("x=" + std::to_string(x)
				    + " is not in range [" + std::to_string(x_min) + "," 
				    + std::to_string(XMax()) + "]");
      }
      if (y < y_min || y > YMax()) {
        throw InvalidSpaceTimeRange("y=" + std::to_string(y)
				    + " is not in range [" + std::to_string(y_min) + "," 
				    + std::to_string(YMax()) + "]");
      }
      if (eta < eta_min || eta > EtaMax()) {
        if (neta == 0 || neta == 1) return;
        throw InvalidSpaceTimeRange("eta=" + std::to_string(eta)
				    + " is not in range [" + std::to_string(eta_min) + "," 
				    + std::to_string(EtaMax()) + "]");
      }
    }

    // get the lower bound of the fluid cell along tau
    /** @return Fluid cell number along the tau-grid.
	@param tau Light-cone coordinate.
    */
    inline int GetIdTau(Jetscape::real tau){
      return(static_cast<int>((tau - tau_min)/dtau));
    }

    // get the lower bound of the fluid cell along x
    /** @return Fluid cell number along the x-grid.           
        @param x Space coordinate.                          
    */
    inline int GetIdX(Jetscape::real x) { 
      return(static_cast<int>((x - x_min)/dx));
    }

    // get the lower bound of the fluid cell along y
    /** @return Fluid cell number along the y-grid.                
        @param y Space coordinate. 
    */
    inline int GetIdY(Jetscape::real y) {
      return(static_cast<int>((y - y_min)/dy));
    }

    // get the lower bound of the fluid cell along eta
    /** @return Fluid cell number along the eta-grid.
        @param eta Light-cone coordinate.
    */
    inline int GetIdEta(Jetscape::real eta) {
      return(static_cast<int>((eta - eta_min)/deta));
    }

    // get the coordinate of tau, x, y, eta on grid
    /** @param id_tau Fluid cell number along tau-grid.
	@return The tau coordinate for fluid cell number.
    */
    inline Jetscape::real TauCoord(int id_tau) { return tau_min + id_tau * dtau; }

    /** @param id_x Fluid cell number along x-grid.
        @return The x coordinate for fluid cell number.
    */
    inline Jetscape::real XCoord(int id_x) { return x_min + id_x * dx; }

    /** @param id_y Fluid cell number along y-grid.
        @return The y coordinate for fluid cell number.
    */
    inline Jetscape::real YCoord(int id_y) { return y_min + id_y * dy; }

    /** @param id_eta Fluid cell number along eta-grid.
        @return The eta coordinate for fluid cell number.
    */
    inline Jetscape::real EtaCoord(int id_eta) { return eta_min + id_eta * deta; }

    // get the FluidCellInfo index in data
    /** @return FluidCellInfo index in the data.
	@param id_tau Fluid cell number along tau-grid.
	@param id_x Fluid cell number along x-grid.
	@param id_y Fluid cell number along y-grid.
	@param id_eta Fluid cell number along eta-grid.
    */
    inline int CellIndex(int id_tau, int id_x, int id_y, int id_eta) {
      return  id_tau * nx * ny * neta + id_x * ny * neta
	+ id_y * neta + id_eta;
    }

    // get the FluidCellInfo at space point given time step
    /** @return FluidCellInfo at a point (id_tau, id_x, id_y, id_etas)
	@param id_tau tau-step number.
	@param x Space coordinate.
	@param y Space coordinate.
	@param eta Light-cone coordinate.
    */
    FluidCellInfo GetFluidCell(int id_tau, int id_x, int id_y, int id_etas);


    // get the FluidCellInfo at space point given time step
    /** @return FluidCellInfo at a point (x,y,eta) and time-step id_tau.
	@param id_tau tau-step number.
	@param x Space coordinate.
	@param y Space coordinate.
	@param eta Light-cone coordinate.
    */
    FluidCellInfo GetAtTimeStep(int id_tau, Jetscape::real x, Jetscape::real y, Jetscape::real etas);

    // get the FluidCellInfo at given space time point
    /** @return FluidCellInfo at a point (tau, x, y, eta).
        @param tau Light-cone coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param eta Light-cone coordinate. 
    */
    FluidCellInfo Get(Jetscape::real tau, Jetscape::real x, Jetscape::real y, Jetscape::real etas);
  };


  class FluidDynamics : public JetScapeModuleBase {


  public:
    /** Default constructor. task ID as "FluidDynamics",  
        eta is initialized to -99.99.
    */  
    FluidDynamics();
    
    /** Standard constructor to create a Fluid Dynamics task. Sets the task ID as "FluidDynamics", and the value of rapidity (eta) to -99.99. 
	@param m_name is a name of the control XML file which contains the input parameters under the tag <Hydro>.
    */
    FluidDynamics(string m_name);
    
    /** Default destructor. */
    virtual ~FluidDynamics();
    
    /** Reads the input parameters from the XML file under the tag <Hydro>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls InitializeHydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.   
	@sa Read about @a polymorphism in C++.
    */
    virtual void Init();
    
    /** Calls EvolveHydro(); This explicit call can be used for actual execution of hydrodynamic evolution defined in the modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
	@sa Read about @a polymorphism in C++.
    */
    virtual void Exec();
    
    /// slots for "jet" signals (future)
    virtual void UpdateEnergyDeposit(int t, double edop);
    /// slots for "jet" signals (future)
    virtual void GetEnergyDensity(int t,double& edensity);

    /** @return parameter_list A pointer to the class Parameter which contains a file name for the fluid dynamics task.
	@sa Implementation of the class Parameter.
    */
    Parameter& GetParameterList() {return parameter_list;}
  
    /** Collect header information for writer modules
	@param w is a pointer of type JetScapeWrite class.
    */
    virtual void CollectHeader( weak_ptr<JetScapeWriter> w );
    
    /** Generated event plane angle
	To be overwritten by implementations that have such information.
    */
    virtual double GetEventPlaneAngle(){ return -999; };
    
    /// Jet signals (future)
    virtual void AddJetSource(double t, double x, double y, double z, JetSource jS) {}; // to be implemented ...

    /** It stores the temperature of the fluid cell at location (t or tau,x,y,z or eta) into an input variable "mT". It can be overridden by modules attached to the FluidDynamics class.
	@param t  tau or t coordinate.
	@param x  space x coordinate. 
	@param y  space y coordinate.
	@param z  rapidity eta or space z coordinate.
	@param mT temperature.
    */
    virtual void GetTemperature(double t, double x, double y, double z, double &mT) {mT=GetTemperature(t,x,y,z);}

    /** It calls GetHydroInfo(t,x,y,z,fCell) to retrieve the properties of the fluid cell at location (t or tau,x,y,z or eta). It can be overridden by modules attached to the FluidDynamics class. 
	@param t  tau or t coordinate.
	@param x  space x coordinate.      
	@param y  space y coordinate.       
	@param z  rapidity eta or space z coordinate.
	@param fCell A pointer of type FluidCellInfo class.  
    */
    virtual void GetHydroCell(double t, double x, double y, double z, std::unique_ptr<FluidCellInfo>& fCell) {GetHydroInfo(t,x,y,z,fCell);} 

    /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in an XML file under the tag <Hydro>.
     */
    tinyxml2::XMLElement* GetHydroXML() {return fd;}

    // currently we have no standard for passing configurations 
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    /** Default function to initialize the hydrodynamics. It can be overridden by different modules. 
	@param parameter_list An object of the class Parameter. */
    virtual void InitializeHydro(Parameter parameter_list) {};
    
    /** Default function to evolve the hydrodynamics. It can be overridden by different modules. */
    virtual void EvolveHydro() {};
    
    /** Default function to evolve the hydrodynamics by one-time (or tau) step. It can be overridden by different modules.
	@param jmu An object to a JetSource class.
    */
    virtual void EvolveHydroOneStep(JetSource jmu) {};

    /** @return Status of the hydrodynamics (NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR). */
    int GetHydroStatus() {return(hydro_status);}

    /** @return Start time (or tau) for hydrodynamic evolution.
     */
    Jetscape::real GetHydroStartTime() {return(hydro_tau_0);}

    /** @return End time (or tau) for hydrodynamic evolution.
     */
    Jetscape::real GetHydroEndTime() {return(hydro_tau_max);}
    /** @return Freeze-out temperature.
     */
    Jetscape::real GetHydroFreezeOutTemperature() {
      return(hydro_freeze_out_temperature);
    }

    /** Retrieves the hydro information at a given space-time point.
     * It throws a InvalidSpaceTimeRange message when
     * (t or tau, x, y, z or eta) is out of the evolution history range. 
     @param time Time or tau coordinate. 
     @param x Space coordinate.
     @param y Space coordinate.
     @param z Space or eta coordinate.
     @param fluid_cell_info_ptr A pointer to the FluidCellInfo class.
    */
    virtual void GetHydroInfo(Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
                                std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr){
      if (hydro_status != FINISHED || !bulk_info.DataSizeMatch()) {
        throw std::runtime_error("Hydro evolution is not finished "
				 "or EvolutionHistory data size does no match description");
      }
      // judge whether to use 2D interpolation or 3D interpolation
      if (!bulk_info.tau_eta_is_tz) {
        Jetscape::real tau = std::sqrt(t * t - z * z);
        Jetscape::real eta = 0.5 * (std::log(t + z) - std::log(t - z));
        bulk_info.CheckInRange(tau, x, y, eta);
        *fluid_cell_info_ptr = bulk_info.Get(tau, x, y, eta);
      } else {
        bulk_info.CheckInRange(t, x, y, z);
        *fluid_cell_info_ptr = bulk_info.Get(t, x, y, z);
      }
    }

    // this function print out the information of the fluid cell to the screen
    /** It prints out the information of the fluid cell.
	@param fluid_cell_info_ptr A pointer to FluidCellInfor class.
    */
    void PrintFluidCellInformation(const std::unique_ptr<FluidCellInfo> & fluid_cell_info_ptr);

    // this function returns hypersurface for Cooper-Frye or recombination
    // the detailed implementation is left to the hydro developper
    /** @return Default function to get the hypersurface for Cooper-Frye or recombination model. It can overridden by different modules.
     */
    virtual void GetHyperSurface(Jetscape::real T_cut,
                                  SurfaceCellInfo* surface_list_ptr) {};

    // all the following functions will call function GetHydroInfo()
    // to get thermaldynamic and dynamical information at a space-time point
    // (time, x, y, z)

    /** @return Energy density at point (t or tau, x, y, z or eta)
        @param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate. 
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetEnergyDensity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    /** @return Entropy density at point (t or tau, x, y, z or eta)
        @param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate. 
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetEntropyDensity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    /** @return Temperature at point (t or tau, x, y, z or eta)
	@param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetTemperature(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    /** @return Fraction of quark gluon plasma assuming medium is in QGP+HRG phase at point (t or tau, x, y, z or eta).
	@param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetQgpFraction(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    // These have no default implementation
    // /** @return 3-component (vx,vy,vz) fluid velocity at point (t or tau, x, y, z or eta).
    //     @param time Time or tau coordinate.
    //     @param x Space coordinate.
    //     @param y Space coordinate.
    //     @param z Space or eta coordinate.
    // */
    // virtual real3 Get3FluidVelocity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z)=0;

    // /** @return 4-component fluid velocity at point (t or tau, x, y, zor eta).
    //     @param time Time or tau coordinate.
    //     @param x Space coordinate.
    //     @param y Space coordinate.
    //     @param z Space or eta coordinate. 
    // */
    // virtual real4 Get4FluidVelocity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    // /** @return Net baryon density at point (t or tau, x, y, z or eta).
    //     @param time Time or tau coordinate.
    //     @param x Space coordinate.
    //     @param y Space coordinate.
    //     @param z Space or eta coordinate. 
    // */
    // virtual Jetscape::real GetNetBaryonDensity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);

    // /** @return Net charge density at point (t or tau, x, y, z or eta).
    // 	@param time Time or tau coordinate.
    // 	@param x Space coordinate.
    // 	@param y Space coordinate.
    // 	@param z Space or eta coordinate.
    // */
    // virtual Jetscape::real GetNetChargeDensity(Jetscape::real time, Jetscape::real x, Jetscape::real y, Jetscape::real z);
    
    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory,
    // for large dataset, std::deque is better than std::vector.
    /** Stores the evolution history. */
    EvolutionHistory bulk_info;

  protected:
    // record hydro start and end proper time [fm/c]
    Jetscape::real hydro_tau_0, hydro_tau_max;
    // record hydro freeze out temperature [GeV]
    Jetscape::real hydro_freeze_out_temperature;
    // record hydro running status
    HydroStatus hydro_status;

    // add initial state shared pointer
    /** A pointer of type InitialState class.
     */
    std::shared_ptr<InitialState> ini;
    std::shared_ptr<PreequilibriumDynamics> pre_eq_ptr;

    double eta;
    Parameter parameter_list;
    tinyxml2::XMLElement *fd;
    
  }; // end class FluidDynamics

} // end namespace Jetscape

#endif  // FLUIDDYNAMICS_H
