/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include <memory>
#include <tuple>

#include "JetClass.h"
#include "JetScapeModuleBase.h"

namespace Jetscape {
/** @class InitialState
   * @brief Interface for Initial State Physics.
   * 
   * Inherits from JetScapeModuleBase.
   */
class InitialState : public JetScapeModuleBase {
 public:
  /** @brief Default constructor to create a Initial State Physics task. 
   * 
   * Sets the task ID as "InitialState".
   */
  InitialState() { SetId("InitialState"); }

  /** 
   * @brief Destructor for the Initial State Physics task.
   */
  ~InitialState();

  /** 
   * @brief Initializes the Initial State Physics task.
   * 
   * Reads the input parameters from the XML file under the tag  <IS>. Calls
     InitTask(); This explicit call of InitTask() can be used for actual
     initialization of modules such as @a Trento if attached as a @a polymorphic
     class. It also initializes the tasks within the current module.
      @sa Read about @a polymorphism in C++.
   */
  virtual void Init();

  /** 
   * @brief Default Exec() function. 
   * Can be overridden by other tasks.
      After this is run, GetNumBinaryCollisions and
     GetEntropyDensityDistribution should return sensible values.
   */
  virtual void Exec();

  /** @brief Default Clear() function. 
   * 
   * Can be overridden by other tasks.
   */
  virtual void Clear();

  /** @brief Default Write() function. 
   * Can be overridden by other tasks.
   * 
   * @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w);

  /** 
   * @brief Collect header information for writer modules
   * 
   * @param w is a pointer of type JetScapeWrite class.
  */
  virtual void CollectHeader(weak_ptr<JetScapeWriter> w);

  /** 
   * @brief Gets generated number of collision participants.
   * 
   * To be overwritten by implementations that have such information.
  */
  virtual double GetNpart() { return -1; };

  /** @brief Gets generated number of binary collisions.
   * 
   * To be overwritten by implementations that have such information.
  */
  virtual double GetNcoll() { return -1; };

  /** @brief Gets generated total entropy
   * 
   * To be overwritten by implementations that have such information.
  */
  virtual double GetTotalEntropy() { return -1; };

  /** @brief Gets generated event centrality
   * 
   * To be overwritten by implementations that have such information.
  */
  virtual double GetEventCentrality() { return -1; };

  // one can set range by hand if not read from xml file
  /** 
   * @brief Sets the range of the coordinates (xmax, ymax, zmax).
   * 
   * @param xmax Maximum value of the coordinate x in the nuclear density
     profile.
     @param ymax Maximum value of the coordinate y in the nuclear density
     profile.
     @param zmax Maxium value of the spatial rapidity ( if (tau,x,y,eta)
     system), or maximum value of the c ordinate z (if in (t,x,y,z) system) in
     the nuclear density profile.
   */
  inline void SetRanges(double xmax, double ymax, double zmax) {
    grid_max_x_ = xmax;
    grid_max_y_ = ymax;
    grid_max_z_ = zmax;
  }

  // one can set grid steps by hand if not read from xml file
  /** 
   * @brief Sets the step-size (dx, dy, dz).
   * 
   * @param dx Step-size for x.
   * @param dy Step-size for y.
   * @param dz Step-size for z or eta.
   */
  inline void SetSteps(double dx, double dy, double dz) {
    grid_step_x_ = dx;
    grid_step_y_ = dy;
    grid_step_z_ = dz;
  }

  /** 
   * @brief Gets the initial state entropy density distribution.
   * @return The initial state entropy density distribution.
       @sa Function CoordFromIdx(int idx) for mapping of the index of the vector
     entropy_density_distribution_ to the fluid cell at location (x, y, z or
     eta).
  */
  inline std::vector<double> GetEntropyDensityDistribution() {
    return entropy_density_distribution_;
  };

  /** 
   * @brief Gets the un-normalized probability density of binary collisions.
   * 
   * One can sample jet production position from Ta * Tb where Ta * Tb is 
     the distribution of num_of_binary_collisions
   * 
   * @return The un-normalized probability density of binary collisions.
   * @sa Function CoordFromIdx(int idx) for mapping of the index of the vector
     num_of_binary_collisions_ to the fluid cell at location (x, y, z or eta).
   */
  inline std::vector<double> GetNumOfBinaryCollisions() {
    return num_of_binary_collisions_;
  };

  //! @return the event id
  int GetEventId() const { return (event_id_); }

  /** 
   * @brief Set the event id
   * @param event_id_in The event id.
   */ 
  void SetEventId(int event_id_in) { event_id_ = event_id_in; }

  /** 
   * @brief compute 3d coordinates (x, y, z) given the 1D index in vector
   * @return Grid point (x,y,z or eta).
   * @param idx is an integer which maps to an unique unit cell in the
     coordinate space (x,y,z or eta).
   */
  std::tuple<double, double, double> CoordFromIdx(int idx);
  virtual void SampleABinaryCollisionPoint(double &x, double &y);

  /** 
   * @return The maximum value of coordinate "x" in the nuclear profile of a nucleus.
   */
  inline double GetXMax() { return grid_max_x_; }

  /** 
   * @return The step-size "dx" used to discretize the nuclear profile of a
     nucleus in x-direction.
   */
  inline double GetXStep() { return grid_step_x_; }

  /** 
   * @return The maximum value of coordinate "y" in the nuclear profile of a
     nucleus.
   */
  inline double GetYMax() { return grid_max_y_; }

  /** 
   * @return The step-size "dy" used to discretize the nuclear profile of a
     nucleus in y-direction.
   */
  inline double GetYStep() { return grid_step_y_; }
  
  /** @return The maximum value of coordinate "z or eta" in the nuclear profile
   * of a nucleus.
   */
  inline double GetZMax() { return grid_max_z_; }

  /** @return The step-size "dz" used to discretize the nuclear profile of a
   * nucleus in z or eta direction.
   */
  inline double GetZStep() { return grid_step_z_; }

  /** 
   * @brief Gets numbef or grids along x-direction following trento convention.
   * @return The number of grid points in x-direction in the nuclear profile of
     a nucleus.
   */
  inline int GetXSize() {
    return int(std::ceil(2 * grid_max_x_ / grid_step_x_));
  }

  
  /** 
   * @brief Gets numbef or grids along y-direction following trento convention.
   * @return The number of grid points in y-direction in the nuclear profile of
     a nucleus.
   */
  inline int GetYSize() {
    return int(std::ceil(2 * grid_max_y_ / grid_step_y_));
  }

  
  /** 
   * @brief Gets numbef or grids along z-direction following trento convention.
   * @return The number of grid points in z or eta direction in the nuclear
     profile of a nucleus.
   */
  inline int GetZSize() {
    int nz = 1;
    if (grid_step_z_ != 0) {
      nz = int(std::ceil(2 * grid_max_z_ / grid_step_z_));
      // for 2d case, user may set grid_max_z_ = 0,
      if (nz == 0)
        nz = 1;
    }
    return nz;
  }

 protected:
  // initial state entropy density distribution for the given grids
  // stored order: for z { for y {for x } }
  /** It stores the initial state entropy density distribution. The index of the
     vector is associated to a cell with coordinate (x, y, z or eta).
      @sa Function CoordFromIdx(int idx) for the mapping.
   */
  std::vector<double> entropy_density_distribution_;

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  // stored order: for z { for y {for x } }
  /** It represents the un-normalized probability density of binary collisions.
     The index of the vector is associated to a cell with coordinate (x, y, z or
     eta).
      @sa Function CoordFromIdx(int idx) for the mapping.
   */
  std::vector<double> num_of_binary_collisions_;
  // the above should be private. Only Adding getters for now to not break other
  // people's code

  /**  @return The initial state entropy density distribution.
       @sa Function CoordFromIdx(int idx) for mapping of the index of the vector
     entropy_density_distribution_ to the fluid cell at location (x, y, z or
     eta).
   */
  // inline std::vector<double> GetEntropyDensityDistribution() {return
  // entropy_density_distribution_;};

  // one can sample jet production position from Ta * Tb
  // where Ta * Tb is the distribution of num_of_binary_collisions
  /** @return The un-normalized probability density of binary collisions.
      @sa Function CoordFromIdx(int idx) for mapping of the index of the vector
     num_of_binary_collisions_ to the fluid cell at location (x, y, z or eta).
   */
  // inline std::vector<double> GetNumOfBinaryCollisions() {return
  // num_of_binary_collisions_;};

  // compute 3d coordinates (x, y, z) given the 1D index in vector
  /** @return Grid point (x,y,z or eta).
      @param idx is an integer which maps to an unique unit cell in the
     coordinate space (x,y,z or eta).
   */
  // std::tuple<double, double, double> CoordFromIdx(int idx);

  int event_id_;
  // int GetEventId() const {return(event_id_);}
  // void SetEventId(int event_id_in) {event_id_ = event_id_in;}

  // default assumption: x range = [-grid_max_x_, grid_max_x_]
  // default assumption: y range = [-grid_max_y_, grid_max_y_]
  // default assumption: z range = [-grid_max_z_, grid_max_z_]
  double grid_max_x_;
  double grid_max_y_;
  double grid_max_z_;  // z can represent spatial rapidity

  // one problem: different hydro codes require different dx_i / dtau
  // matches between (dx_i/dtau, hydro_code) should be checked
  double grid_step_x_;
  double grid_step_y_;
  double grid_step_z_;
};

}  // end namespace Jetscape

#endif
