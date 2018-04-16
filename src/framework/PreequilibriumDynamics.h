/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef PREEQUILDYNAMICS_H
#define PREEQUILDYNAMICS_H

#include<vector>
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "RealType.h"

namespace Jetscape {
    // Flags for preequilibrium dynamics status.
    enum PreequilibriumStatus {NOT_STARTED, INIT, DONE, ERR};

class PreEquilibriumParameterFile {
 public:
    // preequilibrium dynamics parameters file name.
    char* preequilibrium_input_filename;
};

// Interface for the Preequilibrium Dynamics of the medium
class PreequilibriumDynamics : public JetScapeModuleBase {
 private:
    tinyxml2::XMLElement *fd;
    PreEquilibriumParameterFile parameter_list_;
    // record preequilibrium start and end proper time [fm/c]
    real preequilibrium_tau_0_, preequilibrium_TauMax_;

 public:
    PreequilibriumDynamics();
    PreequilibriumDynamics(string m_name) : JetScapeModuleBase(m_name)
    {SetId("PreequilibriumDynamics");}

    virtual ~PreequilibriumDynamics();

    /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls InitializeHydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    void Init();

    /** Calls EvolvePreequilibrium(); This explicit call can be used for actual execution of Preequilibrium evolution defined in the modules such as @a Brick, @a MpiMusic, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    void Exec();

    virtual void InitializePreequilibrium(PreEquilibriumParameterFile parameter_list) {}
    virtual void EvolvePreequilibrium() {}

    /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in an XML file under the tag <Hydro>.
    */
    tinyxml2::XMLElement* GetPreequilibriumXML() {return fd;}

    // add initial state shared pointer
    /** A pointer of type InitialState class.
    */
    std::shared_ptr<InitialState> ini;

    PreEquilibriumParameterFile& GetParameterList() {return parameter_list_;}

    int GetPreequilibriumStatus() {return(preequilibrium_status_);}

    // @return Start time (or tau) for hydrodynamic evolution
    real GetPreequilibriumStartTime() {return(preequilibrium_tau_0_);}

    // @return End time (or tau) for hydrodynamic evolution.
    real GetPreequilibriumEndTime() {return(preequilibrium_TauMax_);}

    // record preequilibrium running status
    PreequilibriumStatus preequilibrium_status_;

    std::vector<double> e_;
    std::vector<double> P_;
    std::vector<double> utau_;
    std::vector<double> ux_;
    std::vector<double> uy_;
    std::vector<double> ueta_;
    std::vector<double> pi00_;
    std::vector<double> pi01_;
    std::vector<double> pi02_;
    std::vector<double> pi03_;
    std::vector<double> pi11_;
    std::vector<double> pi12_;
    std::vector<double> pi13_;
    std::vector<double> pi22_;
    std::vector<double> pi23_;
    std::vector<double> pi33_;
    std::vector<double> bulk_Pi_;
};

}  // end namespace Jetscape

#endif  
