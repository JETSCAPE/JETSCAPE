// Copyright JETSCAPE Collaboration @ 2016
// This is a general basic class for preequilibrium dynamics
// This is written by Derek and ...

#ifndef SRC_PREEQUIL_DYNAMICS_H_
#define SRC_PREEQUIL_DYNAMICS_H_

#include <vector>
#include <cstring>
#include <stdexcept>
#include <cmath>
#include <iostream>

#include "realtype.h"

namespace Jetscape {
  ///Flags for preequilibrium dynamics status.
enum PreequilibriumStatus {NOT_STARTED, INIT, DONE, ERR};
/*
class Parameter{
 public:
  //preequilibrium dynamics parameters file name. 
    char* preequilibrium_input_filename;
};
*/
class PreequilibriumDynamicsBase{
 protected:
    // record preequilibrium start and end proper time [fm/c]
    real preequilibrium_tau_0, preequilibrium_tau_max;
    // record preequilibrium running status
    PreequilibriumStatus preequilibrium_status;

 public:
    // Default constructor.
    PreequilibriumDynamicsBase() {};

    // Default destructor.
    ~PreequilibriumDynamicsBase() {};

    // How to store this data? In memory or hard disk?
    // 3D hydro may eat out the memory,
    // for large dataset, std::deque is better than std::vector.
    //Stores the evolution history.
    //EvolutionHistory bulk_info;

    //Keep this interface open in the beginning.
    // FreezeOutHyperSf hyper_sf;

    // currently we have no standard for passing configurations
    // pure virtual function; to be implemented by users
    // should make it easy to save evolution history to bulk_info
    // Default function to initialize the hydrodynamics. It can be overridden by different modules.
    // @param parameter_list An object of the class Parameter.
    virtual void initialize_preequilibrium(Parameter parameter_list) {};

    /** Default function to evolve the hydrodynamics. It can be overridden by different modules. */
    virtual void evolve_preequilibrium() {};

    //Default function to evolve the hydrodynamics by one-time (or tau) step. It can be overridden by different modules.
	//@param jmu An object to a JetSource class.

    //virtual void evolve_hydro_one_step(JetSource jmu) {};

    // the following functions should be implemented in Jetscape
    // @return Status of the hydrodynamics (NOT_START, INITIALIZED, EVOLVING, FINISHED, ERROR).
    int get_preequilibrium_status() {return(preequilibrium_status);}

    // @return Start time (or tau) for hydrodynamic evolution
    real get_preequilibrium_start_time() {return(preequilibrium_tau_0);}

    // @return End time (or tau) for hydrodynamic evolution.
    real get_preequilibrium_end_time() {return(preequilibrium_tau_max);}

}; // end class PreequilibriumDynamicsBase

} // end namespace Jetscape

#endif  // SRC_PREEQUIL_DYNAMICS_H_
