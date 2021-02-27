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

#ifndef NCOLLLISTFROMFILE_H
#define NCOLLLISTFROMFILE_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include "JetScapeModuleBase.h"
#include "InitialState.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

class NcollListFromFile : public Jetscape::InitialState {
  // this is wrapper class to read external files that
  // stores initial number of binary collisions and corresponding
  // configurations
public:
  NcollListFromFile();
  ~NcollListFromFile();

  //void Init();

  /** Default Exec() function. It can be overridden by other tasks.
   */
  void Exec();

  /** Default Clear() function. It can be overridden by other tasks.
   */
  void Clear();

  /** Generated number of binary collisions. */
  double GetNcoll() { return(ncoll_); };

  //! Load saved number of binary collisions
  void ReadNbcList();

  void SampleABinaryCollisionPoint(double &x, double &y) {

private:
  std::vector<double> binary_collision_x_;
  std::vector<double> binary_collision_y_;
  std::shared_ptr<std::uniform_int_distribution<int>> rand_int_ptr_;

  double ncoll_ = -1;

  // Allows the registration of the module so that it is available to be used
  // by the Jetscape framework.
  static RegisterJetScapeModule<InitialFromFile> reg;
};

#endif  // NCOLLLISTFROMFILE_H
