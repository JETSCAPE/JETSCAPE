// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef PREEQUILDYNAMICS_H
#define PREEQUILDYNAMICS_H

#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "preequilibrium_dynamics.h"

//JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING WORKS !
/*
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
*/
//JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING FILES WORKS !

namespace Jetscape {
  /**
  @class
  Interface for the Preequilibrium Dynamics of the medium
  */
  class PreequilibriumDynamics : public JetScapeModuleBase , public PreequilibriumDynamicsBase
  {

  public:

    /** Default constructor.
    */
    PreequilibriumDynamics();

    /** Standard constructor to create a Preequilibrium Dynamics task. Sets the task ID as "PreequilibriumDynamics".
    @param m_name is a name of the control XML file which contains the input parameters under the tag <Preequilibrium>.
    */
    PreequilibriumDynamics(string m_name) : JetScapeModuleBase (m_name), PreequilibriumDynamicsBase()
    {SetId("PreequilibriumDynamics");}

    /** Destructor for the Preequilibrium Dynamics task.
    */
    virtual ~PreequilibriumDynamics();

    /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls initialize_hydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    virtual void Init();

    /** Calls evolve_preequilibrium(); This explicit call can be used for actual execution of Preequilibrium evolution defined in the modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
    @sa Read about @a polymorphism in C++.
    */
    virtual void Exec();

    /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in an XML file under the tag <Hydro>.
    */
    tinyxml2::XMLElement* GetPreequilibriumXML() {return fd;}

    // add initial state shared pointer
    /** A pointer of type InitialState class.
    */
    std::shared_ptr<InitialState> ini;

    Parameter& GetParameterList() {return parameter_list;}

  private:

    //double eta;
    tinyxml2::XMLElement *fd;
    Parameter parameter_list;

  };

  //JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING FILES WORKS !
  /*
  PreequilibriumDynamics::PreequilibriumDynamics()
  {
    VERBOSE(8);
    SetId("PreequilibriumDynamics");
  }

  PreequilibriumDynamics::~PreequilibriumDynamics()
  {
    VERBOSE(8);
    disconnect_all();
  }

  void PreequilibriumDynamics::Init()
  {
    JetScapeModuleBase::Init();

    INFO<<"Intialize PreequilibriumDynamics : "<<GetId()<< " ...";

    fd= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Preequilibrium" );

    if (!fd) {
      WARN << "Not a valid JetScape XML Preequilibrium Dynamics section file or no XML file loaded!";
      exit(-1);
    }

    VERBOSE(8);

    ini = JetScapeSignalManager::Instance()->GetInitialStatePointer().lock();
    if (!ini) {
      WARN << "No initialization module, try: auto trento = make_shared<TrentoInitial>(); jetscape->Add(trento);";
    }

    initialize_preequilibrium(parameter_list);

    InitTask();

    JetScapeTask::InitTasks();
  }

  void PreequilibriumDynamics::Exec()
  {
    INFO <<"Run Preequilibrium : "<<GetId()<< " ...";
    VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

    if (ini) {
      INFO << "length of energy density vector=" << ini->entropy_density_distribution_.size();
    }

    evolve_preequilibrium();

    JetScapeTask::ExecuteTasks();
  }
  */
  //JUST FOR TESTING! REMOVE THIS BLOCK ONCE CMAKE LINKING FILES WORKS !

} // end namespace Jetscape

#endif
