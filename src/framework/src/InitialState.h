// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include <memory>
#include "JetScapeModuleBase.h"
#include "JetClass.hpp"
#include "tinyxml2.h"
#include "JetScapeXML.h"

namespace Jetscape {

class InitialState : public JetScapeModuleBase {
 public:

  InitialState();
  InitialState(string m_name): JetScapeModuleBase(m_name){ SetId("InitialState"); }
  ~InitialState();

  virtual void Init();
  virtual void Exec();
  virtual void Clear();

  virtual void Write(weak_ptr<JetScapeWriter> w);

  tinyxml2::XMLElement * GetIniStateXML() { return xml_; }

  inline double get_x_max(){ return grid_max_x_; }
  inline double get_x_step(){ return grid_step_x_; }
  inline double get_y_max(){ return grid_max_y_; }
  inline double get_y_step(){ return grid_step_y_; }
  inline double get_z_max(){ return grid_max_z_; }
  inline double get_z_step(){ return grid_step_z_; }

  // initial state entropy density distribution for the given grids
  std::vector<double> entropy_density_distribution_;

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  std::vector<double> num_of_binary_collisions_;

  // xml_ reads the xml configurations for initial states
  tinyxml2::XMLElement * xml_;
 private:
  // default assumption: x range = [-grid_max_x_, grid_max_x_]
  // default assumption: y range = [-grid_max_y_, grid_max_y_]
  // default assumption: z range = [-grid_max_z_, grid_max_z_]
  double grid_max_x_;
  double grid_max_y_;
  double grid_max_z_;     // z can represent spatial rapidity

  // one problem: different hydro codes require different dx_i / dtau
  // matches between (dx_i/dtau, hydro_code) should be checked
  double grid_step_x_;
  double grid_step_y_;
  double grid_step_z_;

};

} // end namespace Jetscape

#endif
