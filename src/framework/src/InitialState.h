// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include "JetScapeModuleBase.h"
#include "JetClass.hpp"

namespace Jetscape {

class InitialState : public JetScapeModuleBase
{
  
 public:

  InitialState();
  InitialState(double collEnergy, Vertex vtx);
  InitialState(string m_name): JetScapeModuleBase(m_name){ SetId("InitialState"); }
  ~InitialState();

  virtual void Init();
  virtual void Exec();
  virtual void Clear();

  virtual void Write(weak_ptr<JetScapeWriter> w);

  inline void setTemperature(double temp) { temperature = temp; }
  inline double getTemperature() {return temperature; }

 private:
    double collEnergy;
    double temperature;
    Vertex initialVtx;

};

} // end namespace Jetscape

#endif
