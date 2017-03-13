// Framework test (dummy) FluidDynamics class (to be changed with real implemenation)

#ifndef FLUIDDYNAMICS_H
#define FLUIDDYNAMICS_H

#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "fluid_dynamics.h"
//#include "brick_jetscape.h"

// in principle one class level too much ... think about and compress ...

class FluidDynamics : public JetScapeModuleBase , public FluidDynamicsBase //, std::enable_shared_from_this<FluidDynamics>
{
  
 public:
  
  FluidDynamics();
  FluidDynamics(string m_name) : JetScapeModuleBase (m_name), FluidDynamicsBase()
    {eta=-99.99; SetId("FluidDynamics");}
  virtual ~FluidDynamics();

  virtual void Init();
  virtual void Exec();

  tinyxml2::XMLElement* GetHydroXML() {return fd;}
 
  // slots for "jet" signals (will be obsolete ...)
  void UpdateEnergyDeposit(int t, double edop);
  void GetEnergyDensity(int t,double& edensity);
  Parameter& GetParameterList() {return parameter_list;}
  
  // real slots based on FluidDymanics(Test) class
  virtual void AddJetSource(double t, double x, double y, double z, JetSource jS) {}; // check ...
  virtual void GetTemperature(double t, double x, double y, double z, double &mT) {mT=get_temperature(t,x,y,z);}
  virtual void GetHydroCell(double t, double x, double y, double z, FluidCellInfo* fCell) {get_hydro_info(t,x,y,z,fCell);} 
  
 private:

  double eta;
  tinyxml2::XMLElement *fd;
  Parameter parameter_list;
  
};

#endif
