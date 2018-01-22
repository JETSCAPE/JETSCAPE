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

namespace Jetscape {
// in principle one class level too much ... think about and compress ...
   /**
       @class
       Interface for the Fluid Dynamics of the medium
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
    {SetId("FluidDynamics");}

  /** Destructor for the Preequilibrium Dynamics task.
  */
  virtual ~PreequilibriumDynamics();

  /** Reads the input parameters from the XML file under the tag <Preequilibrium>. Uses JetScapeSingnalManager Instance to retrive the Initial State Physics information. Calls initialize_hydro(parameter_list) and InitTask(); This explicit call can be used for actual initialization of modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
  */
  virtual void Init();

  /** Calls evolve_hydro(); This explicit call can be used for actual execution of hydrodynamic evolution defined in the modules such as @a Brick, @a MPI_MUSIC, or @a OSU-HYDRO if attached as a @a polymorphic class. It also execute the tasks within the current module.
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

  // slots for "jet" signals (will be obsolete ...)
  //void UpdateEnergyDeposit(int t, double edop);
  //void GetEnergyDensity(int t,double& edensity);
  /** @return parameter_list A pointer to the class Parameter which contains a file name for the fluid dynamics task.
      @sa Implementation of the class Parameter.
   */
  Parameter& GetParameterList() {return parameter_list;}

  // real slots based on PreequilibriumDynamics(Test) class
  //virtual void AddJetSource(double t, double x, double y, double z, JetSource jS) {}; // to be implemented ...
  /** It stores the temperature of the fluid cell at location (t or tau,x,y,z or eta) into an input variable "mT". It can be overridden by modules attached to the FluidDynamics class.
      @param t  tau or t coordinate.
      @param x  space x coordinate.
      @param y  space y coordinate.
      @param z  rapidity eta or space z coordinate.
      @param mT temperature.
   */
  //virtual void GetTemperature(double t, double x, double y, double z, double &mT) {mT=get_temperature(t,x,y,z);}
  /** It calls get_hydro_info(t,x,y,z,fCell) to retrieve the properties of the fluid cell at location (t or tau,x,y,z or eta). It can be overridden by modules attached to the FluidDynamics class.
      @param t  tau or t coordinate.
      @param x  space x coordinate.
      @param y  space y coordinate.
      @param z  rapidity eta or space z coordinate.
      @param fCell A pointer of type FluidCellInfo class.
   */
  //virtual void GetHydroCell(double t, double x, double y, double z, FluidCellInfo* fCell) {get_hydro_info(t,x,y,z,fCell);}

 private:

  //double eta;
  tinyxml2::XMLElement *fd;
  Parameter parameter_list;

};

} // end namespace Jetscape

#endif
