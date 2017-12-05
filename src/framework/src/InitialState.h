// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include <tuple>
#include <memory>
#include "JetScapeModuleBase.h"
#include "JetClass.hpp"
#include "tinyxml2.h"
#include "JetScapeXML.h"

namespace Jetscape {
  /** @class
      Interface for Initial State Physics.
   */
class InitialState : public JetScapeModuleBase {
 public:

  /** Default constructor.      
   */
  InitialState();
  /** Standard constructor to create a Initial State Physics task. Sets the task ID as "InitialState". 
      @param m_name is a name of the control XML file which contains the input parameters under the tag <IS>. 
   */
  InitialState(string m_name): JetScapeModuleBase(m_name){ SetId("InitialState"); }
  
  /**  Destructor for the Initial State Physics task.
   */
  ~InitialState();

  /** Reads the input parameters from the XML file under the tag  <IS>. Calls InitTask(); This explicit call of InitTask() can be used for actual initialization of modules such as @a Trento if attached as a @a polymorphic class. It also initializes the tasks within the current module.
      @sa Read about @a polymorphism in C++.
   */
  virtual void Init();

  /** Default Exec() function. It can be overridden by other tasks.
   */
  virtual void Exec();
  
  /** Default Clear() function. It can be overridden by other tasks.
   */
  virtual void Clear();

  /** Default Write() function. It can be overridden by other tasks.
      @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w);

  /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in the XML file under the tag <IS>.
   */
  tinyxml2::XMLElement * GetIniStateXML() { return xml_; }

  // one can set range by hand if not read from xml file 
  /** Sets the range of the coordinates (xmax, ymax, zmax). 
      @param xmax Maximum value of the coordinate x in the nuclear density profile.
      @param ymax Maximum value of the coordinate y in the nuclear density profile.
      @param zmax Maxium value of the spatial rapidity ( if (tau,x,y,eta) system), or maximum value of the coordinate z (if in (t,x,y,z) system) in the nuclear density profile.  
   */
  inline void set_ranges(double xmax, double ymax, double zmax) {
      grid_max_x_ = xmax;
      grid_max_y_ = ymax;
      grid_max_z_ = zmax;
  }

  // one can set grid steps by hand if not read from xml file 
  /** It sets the step-size (dx, dy, dz).
	@param dx Step-size for x.
	@param dy Step-size for y.
	@param dz Step-size for z or eta.
   */
  inline void set_steps(double dx, double dy, double dz) {
      grid_step_x_ = dx;
      grid_step_y_ = dy;
      grid_step_z_ = dz;
  }
  /**  @return The maximum value of coordinate "x" in the nuclear profile of a nucleus.
   */
  inline double get_x_max(){ return grid_max_x_; }
  /** @return The step-size "dx" used to discretize the nuclear profile of a nucleus in x-direction.
   */
  inline double get_x_step(){ return grid_step_x_; }
  /** @return The maximum value of coordinate "y" in the nuclear profile of a nucleus.                       
  */
  inline double get_y_max(){ return grid_max_y_; }
  /** @return The step-size "dy" used to discretize the nuclear profile of a nucleus in y-direction. 
  */
  inline double get_y_step(){ return grid_step_y_; }
  /** @return The maximum value of coordinate "z or eta" in the nuclear profile of a nucleus.                          
  */
  inline double get_z_max(){ return grid_max_z_; }
  /** @return The step-size "dz" used to discretize the nuclear profile of a nucleus in z or eta direction.              
  */
  inline double get_z_step(){ return grid_step_z_; }

  // get number of grids along x, follow trento convention
  /** @return The number of grid points in x-direction in the nuclear profile of a nucleus.
   */
  inline int get_x_size() {
      return int(std::ceil(2 * grid_max_x_ / grid_step_x_));
  }

  // get number of grids along y
  /** @return The number of grid points in y-direction in the nuclear profile of a nucleus.
   */
  inline int get_y_size() {
      return int(std::ceil(2 * grid_max_y_ / grid_step_y_));
  }

  // get number of grids along z
  /** @return The number of grid points in z or eta direction in the nuclear profile of a nucleus.
   */
  inline int get_z_size() {
      int nz = 1;
      if (grid_step_z_ != 0) {
          int nz = int(std::ceil(2 * grid_max_z_ / grid_step_z_));
          // for 2d case, user may set grid_max_z_ = 0,
          if (nz == 0) nz = 1;
      }
      return nz;
  }

  // initial state entropy density distribution for the given grids
  // stored order: for z { for y {for x } }
  /** It stores the initial state entropy density distribution. The index of the vector is associated to a cell with coordinate (x, y, z or eta).
      @sa Function coord_from_idx(int idx) for the mapping.
   */
  std::vector<double> entropy_density_distribution_;

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  // stored order: for z { for y {for x } }
  /** It represents the un-normalized probability density of binary collisions. The index of the vector is associated to a cell with coordinate (x, y, z or eta).
      @sa Function coord_from_idx(int idx) for the mapping.
   */
  std::vector<double> num_of_binary_collisions_;

  // the above should be private. Only Adding getters for now to not break other people's code

  /**  @return The initial state entropy density distribution.
       @sa Function coord_from_idx(int idx) for mapping of the index of the vector entropy_density_distribution_ to the fluid cell at location (x, y, z or eta).
   */
  inline std::vector<double> get_entropy_density_distribution() {return entropy_density_distribution_;};

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  /** @return The un-normalized probability density of binary collisions.
      @sa Function coord_from_idx(int idx) for mapping of the index of the vector num_of_binary_collisions_ to the fluid cell at location (x, y, z or eta).
   */
  inline std::vector<double> get_num_of_binary_collisions() {return num_of_binary_collisions_;};
    
  // compute 3d coordinates (x, y, z) given the 1D index in vector
  /** @return Grid point (x,y,z or eta). 
      @param idx is an integer which maps to an unique unit cell in the coordinate space (x,y,z or eta). 
   */
  std::tuple<double, double, double> coord_from_idx(int idx);

  // xml_ reads the xml configurations for initial states
  /** It is used to access the input parameters from the XML file relevant to the initial state physics.
   */
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
