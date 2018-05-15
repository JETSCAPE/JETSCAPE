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

#ifndef INITIALSTATE_H
#define INITIALSTATE_H

#include <tuple>
#include <memory>
#include "JetScapeModuleBase.h"
#include "JetClass.h"
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
      After this is run, GetNumBinaryCollisions and GetEntropyDensityDistribution should return sensible values.
   */
  virtual void Exec();
  
  /** Default Clear() function. It can be overridden by other tasks.
   */
  virtual void Clear();

  /** Default Write() function. It can be overridden by other tasks.
      @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w);

  /** Collect header information for writer modules
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void CollectHeader( weak_ptr<JetScapeWriter> w );

  /** Generated number of collision participants.
      To be overwritten by implementations that have such information.
  */
  virtual double GetNpart(){ return -1; };

  /** Generated number of binary collisions.
      To be overwritten by implementations that have such information.
  */
  virtual double GetNcoll(){ return -1; };

  /** Generated total entropy
      To be overwritten by implementations that have such information.
  */
  virtual double GetTotalEntropy(){ return -1; };

  /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in the XML file under the tag <IS>.
   */
  tinyxml2::XMLElement * GetIniStateXML() { return xml_; }

  // one can set range by hand if not read from xml file 
  /** Sets the range of the coordinates (xmax, ymax, zmax). 
      @param xmax Maximum value of the coordinate x in the nuclear density profile.
      @param ymax Maximum value of the coordinate y in the nuclear density profile.
      @param zmax Maxium value of the spatial rapidity ( if (tau,x,y,eta) system), or maximum value of the coordinate z (if in (t,x,y,z) system) in the nuclear density profile.  
   */
  inline void SetRanges(double xmax, double ymax, double zmax) {
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
  inline void SetSteps(double dx, double dy, double dz) {
      grid_step_x_ = dx;
      grid_step_y_ = dy;
      grid_step_z_ = dz;
  }

  /**  @return The initial state entropy density distribution.
       @sa Function CoordFromIdx(int idx) for mapping of the index of the vector entropy_density_distribution_ to the fluid cell at location (x, y, z or eta).
  */
  inline std::vector<double> GetEntropyDensityDistribution() {return entropy_density_distribution_;};
  
  /** one can sample jet production position from Ta * Tb
      where Ta * Tb is the distribution of num_of_binary_collisions
      @return The un-normalized probability density of binary collisions.
      @sa Function CoordFromIdx(int idx) for mapping of the index of the vector num_of_binary_collisions_ to the fluid cell at location (x, y, z or eta).
   */
  inline std::vector<double> GetNumOfBinaryCollisions() {return num_of_binary_collisions_;};

  //! @return the event id 
  int GetEventId() const {return(event_id_);}

  //! set the event id 
  void SetEventId(int event_id_in) {event_id_ = event_id_in;}
 
  /** compute 3d coordinates (x, y, z) given the 1D index in vector
      @return Grid point (x,y,z or eta). 
      @param idx is an integer which maps to an unique unit cell in the coordinate space (x,y,z or eta). 
   */
  std::tuple<double, double, double> CoordFromIdx(int idx);

  /**  @return The maximum value of coordinate "x" in the nuclear profile of a nucleus.
   */
  inline double GetXMax(){ return grid_max_x_; }
  /** @return The step-size "dx" used to discretize the nuclear profile of a nucleus in x-direction.
   */
  inline double GetXStep(){ return grid_step_x_; }
  /** @return The maximum value of coordinate "y" in the nuclear profile of a nucleus.                       
  */
  inline double GetYMax(){ return grid_max_y_; }
  /** @return The step-size "dy" used to discretize the nuclear profile of a nucleus in y-direction. 
  */
  inline double GetYStep(){ return grid_step_y_; }
  /** @return The maximum value of coordinate "z or eta" in the nuclear profile of a nucleus.                          
  */
  inline double GetZMax(){ return grid_max_z_; }
  /** @return The step-size "dz" used to discretize the nuclear profile of a nucleus in z or eta direction.              
  */
  inline double GetZStep(){ return grid_step_z_; }

  // get number of grids along x, follow trento convention
  /** @return The number of grid points in x-direction in the nuclear profile of a nucleus.
   */
  inline int GetXSize() {
      return int(std::ceil(2 * grid_max_x_ / grid_step_x_));
  }

  // get number of grids along y
  /** @return The number of grid points in y-direction in the nuclear profile of a nucleus.
   */
  inline int GetYSize() {
      return int(std::ceil(2 * grid_max_y_ / grid_step_y_));
  }

  // get number of grids along z
  /** @return The number of grid points in z or eta direction in the nuclear profile of a nucleus.
   */
  inline int GetZSize() {
      int nz = 1;
      if (grid_step_z_ != 0) {
          int nz = int(std::ceil(2 * grid_max_z_ / grid_step_z_));
          // for 2d case, user may set grid_max_z_ = 0,
          if (nz == 0) nz = 1;
      }
      return nz;
  }

 protected:

  // initial state entropy density distribution for the given grids
  // stored order: for z { for y {for x } }
  /** It stores the initial state entropy density distribution. The index of the vector is associated to a cell with coordinate (x, y, z or eta).
      @sa Function CoordFromIdx(int idx) for the mapping.
   */
  std::vector<double> entropy_density_distribution_;

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  // stored order: for z { for y {for x } }
  /** It represents the un-normalized probability density of binary collisions. The index of the vector is associated to a cell with coordinate (x, y, z or eta).
      @sa Function CoordFromIdx(int idx) for the mapping.
   */
  std::vector<double> num_of_binary_collisions_;
  // the above should be private. Only Adding getters for now to not break other people's code

  /**  @return The initial state entropy density distribution.
       @sa Function CoordFromIdx(int idx) for mapping of the index of the vector entropy_density_distribution_ to the fluid cell at location (x, y, z or eta).
   */
  //inline std::vector<double> GetEntropyDensityDistribution() {return entropy_density_distribution_;};

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  /** @return The un-normalized probability density of binary collisions.
      @sa Function CoordFromIdx(int idx) for mapping of the index of the vector num_of_binary_collisions_ to the fluid cell at location (x, y, z or eta).
   */
  //inline std::vector<double> GetNumOfBinaryCollisions() {return num_of_binary_collisions_;};
    
  // compute 3d coordinates (x, y, z) given the 1D index in vector
  /** @return Grid point (x,y,z or eta). 
      @param idx is an integer which maps to an unique unit cell in the coordinate space (x,y,z or eta). 
   */
  //std::tuple<double, double, double> CoordFromIdx(int idx);

  // xml_ reads the xml configurations for initial states
  /** It is used to access the input parameters from the XML file relevant to the initial state physics.
   */
  tinyxml2::XMLElement * xml_;
    
  int event_id_;
  //int GetEventId() const {return(event_id_);}
  //void SetEventId(int event_id_in) {event_id_ = event_id_in;}

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
