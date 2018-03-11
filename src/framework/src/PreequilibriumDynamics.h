// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef PREEQUILDYNAMICS_H
#define PREEQUILDYNAMICS_H

#include<vector>
#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "realtype.h"

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
    real preequilibrium_tau_0_, preequilibrium_tau_max_;

 public:
    PreequilibriumDynamics();
    PreequilibriumDynamics(string m_name) : JetScapeModuleBase(m_name)
    {SetId("PreequilibriumDynamics");}

    virtual ~PreequilibriumDynamics();

    /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls initialize_hydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    void Init();

    /** Calls evolve_preequilibrium(); This explicit call can be used for actual execution of Preequilibrium evolution defined in the modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    void Exec();

    virtual void initialize_preequilibrium(PreEquilibriumParameterFile parameter_list) {}
    virtual void evolve_preequilibrium() {}

    /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in an XML file under the tag <Hydro>.
    */
    tinyxml2::XMLElement* GetPreequilibriumXML() {return fd;}

    // add initial state shared pointer
    /** A pointer of type InitialState class.
    */
    std::shared_ptr<InitialState> ini;

    PreEquilibriumParameterFile& GetParameterList() {return parameter_list_;}

    int get_preequilibrium_status() {return(preequilibrium_status_);}

    // @return Start time (or tau) for hydrodynamic evolution
    real get_preequilibrium_start_time() {return(preequilibrium_tau_0_);}

    // @return End time (or tau) for hydrodynamic evolution.
    real get_preequilibrium_end_time() {return(preequilibrium_tau_max_);}

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
    std::vector<double> pi23_;
    std::vector<double> pi33_;
    std::vector<double> bulk_Pi_;
};

}  // end namespace Jetscape

#endif  
