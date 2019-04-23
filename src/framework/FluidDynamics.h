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
#include <array>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <iostream>

#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "PreequilibriumDynamics.h"
#include "RealType.h"
#include "FluidCellInfo.h"
#include "SurfaceCellInfo.h"
#include "FluidEvolutionHistory.h"
#include "LiquefierBase.h"

namespace Jetscape {

/// Flags for hydrodynamics status.  
enum HydroStatus {NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR};

/** A helper class for hydro parameters file name.*/
class Parameter {
 public:
    char* hydro_input_filename;
};


class FluidDynamics : public JetScapeModuleBase {
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
    
    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory,
    // for large dataset, std::deque is better than std::vector.
    /** Stores the evolution history. */
    EvolutionHistory bulk_info;

    std::weak_ptr<LiquefierBase> liquefier_ptr;

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

    virtual void Clear();
    
    /** @return parameter_list A pointer to the class Parameter which contains a file name for the fluid dynamics task.
	@sa Implementation of the class Parameter.
    */
    Parameter& GetParameterList() {return(parameter_list);}
  
    /** Collect header information for writer modules
	@param w is a pointer of type JetScapeWrite class.
    */
    virtual void CollectHeader( weak_ptr<JetScapeWriter> w );
    
    /** Generated event plane angle
	To be overwritten by implementations that have such information.
    */
    virtual double GetEventPlaneAngle() {return(-999);}
    
    /** It stores the temperature of the fluid cell at location (t or tau,x,y,z or eta) into an input variable "mT". It can be overridden by modules attached to the FluidDynamics class.
	@param t  tau or t coordinate.
	@param x  space x coordinate. 
	@param y  space y coordinate.
	@param z  rapidity eta or space z coordinate.
	@param mT temperature.
    */
    virtual void GetTemperature(double t, double x, double y, double z,
                                double &mT) {
        mT = GetTemperature(t,x,y,z);
    }

    /** It calls GetHydroInfo(t,x,y,z,fCell) to retrieve the properties of the fluid cell at location (t or tau,x,y,z or eta). It can be overridden by modules attached to the FluidDynamics class. 
	@param t  tau or t coordinate.
	@param x  space x coordinate.      
	@param y  space y coordinate.       
	@param z  rapidity eta or space z coordinate.
	@param fCell A pointer of type FluidCellInfo class.  
    */
    virtual void GetHydroCell(double t, double x, double y, double z,
                              std::unique_ptr<FluidCellInfo>& fCell) {
        GetHydroInfo(t, x, y, z, fCell);
    } 

    /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in an XML file under the tag <Hydro>.
     */
    tinyxml2::XMLElement* GetHydroXML() {return(fd);}

    // currently we have no standard for passing configurations 
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    /** Default function to initialize the hydrodynamics. It can be overridden by different modules. 
	@param parameter_list An object of the class Parameter. */
    virtual void InitializeHydro(Parameter parameter_list) {};
    
    /** Default function to evolve the hydrodynamics. It can be overridden by different modules. */
    virtual void EvolveHydro() {};
    

    /** @return Status of the hydrodynamics (NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR). */
    int GetHydroStatus() const {return(hydro_status);}

    void StoreHydroEvolutionHistory(
                        std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr) {
        bulk_info.data.push_back(*fluid_cell_info_ptr);
    }

    void clear_up_evolution_data() {bulk_info.clear_up_evolution_data();}

    /** @return Start time (or tau) for hydrodynamic evolution.
     */
    Jetscape::real GetHydroStartTime() const {return(hydro_tau_0);}

    /** @return End time (or tau) for hydrodynamic evolution.
     */
    Jetscape::real GetHydroEndTime() const {return(hydro_tau_max);}
    /** @return Freeze-out temperature.
     */
    Jetscape::real GetHydroFreezeOutTemperature() const {
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
    virtual void GetHydroInfo(
        Jetscape::real t, Jetscape::real x, Jetscape::real y, Jetscape::real z,
        std::unique_ptr<FluidCellInfo>& fluid_cell_info_ptr) {
        if (hydro_status != FINISHED || bulk_info.data.size() == 0) {
	        throw std::runtime_error("Hydro evolution is not finished "
				                     "or EvolutionHistory is empty");
        }
        // judge whether to use 2D interpolation or 3D interpolation
        if (!bulk_info.tau_eta_is_tz) {
	        Jetscape::real tau = std::sqrt(t * t - z * z);
	        Jetscape::real eta = 0.5 * (std::log(t + z) - std::log(t - z));
	        bulk_info.CheckInRange(tau, x, y, eta);
	        //return bulk_info.get(tau, x, y, eta);
        } else {
	        bulk_info.CheckInRange(t, x, y, z);
	        //return bulk_info.get(t, x, y, z);
        }
    }

    // this function print out the information of the fluid cell to the screen
    /** It prints out the information of the fluid cell.
	@param fluid_cell_info_ptr A pointer to FluidCellInfor class.
    */
    void PrintFluidCellInformation(FluidCellInfo* fluid_cell_info_ptr);

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
    virtual Jetscape::real GetEnergyDensity(
                Jetscape::real time, Jetscape::real x,
                Jetscape::real y, Jetscape::real z);

    /** @return Entropy density at point (t or tau, x, y, z or eta)
        @param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate. 
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetEntropyDensity(
                Jetscape::real time, Jetscape::real x,
                Jetscape::real y, Jetscape::real z);

    /** @return Temperature at point (t or tau, x, y, z or eta)
	@param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetTemperature(
                Jetscape::real time, Jetscape::real x,
                Jetscape::real y, Jetscape::real z);

    /** @return Fraction of quark gluon plasma assuming medium is in QGP+HRG phase at point (t or tau, x, y, z or eta).
	@param time Time or tau coordinate.
        @param x Space coordinate.
        @param y Space coordinate.
        @param z Space or eta coordinate.
    */
    virtual Jetscape::real GetQgpFraction(
                Jetscape::real time, Jetscape::real x,
                Jetscape::real y, Jetscape::real z);

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
    
    virtual void add_a_liqueifier(
                            std::shared_ptr<LiquefierBase> new_liqueifier) {
        liquefier_ptr = new_liqueifier;
    }

    void get_source_term(Jetscape::real tau, Jetscape::real x,
                         Jetscape::real y, Jetscape::real eta,
                         std::array<Jetscape::real, 4> jmu) const;

    /// slots for "jet" signals (future)
    virtual void UpdateEnergyDeposit(int t, double edop);
    /// slots for "jet" signals (future)
    virtual void GetEnergyDensity(int t, double& edensity) {edensity = 0.0;}

}; // end class FluidDynamics
  
} // end namespace Jetscape

#endif  // FLUIDDYNAMICS_H
