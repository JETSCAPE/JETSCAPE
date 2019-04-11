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

#ifndef INITIALFROMFILE_H
#define INITIALFROMFILE_H

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "hdf5.h"
#include "Hydroinfo_h5.h"
#include "JetScapeModuleBase.h"
#include "InitialState.h"
#include "JetScapeLogger.h"

using namespace Jetscape;

class InitialFromFile: public Jetscape::InitialState {
  // this is wrapper class to read external files that 
  // stores initial number of binary collisions and corresponding
  // configurations
  public:
  InitialFromFile();
  ~InitialFromFile();

  /** Reads the input parameters from the XML file under the tag  <IS>. Calls InitTask(); This explicit call of InitTask() can be used for actual initialization of modules such as @a Trento if attached as a @a polymorphic class. It also initializes the tasks within the current module.
      @sa Read about @a polymorphism in C++.
   */
  //void Init();

  /** Default Exec() function. It can be overridden by other tasks.
   */
  void Exec();
  
  /** Default Clear() function. It can be overridden by other tasks.
   */
  void Clear();

  void InitTask();

  /** Default Write() function. It can be overridden by other tasks.
      @param w A pointer to the JetScapeWriter class.
   */
  virtual void Write(weak_ptr<JetScapeWriter> w);

  /** Generated number of collision participants.
  */
  double GetNpart(){ return npart; };

  /** Generated number of binary collisions.
  */
  double GetNcoll(){ return ncoll; };

  /** Generated total entropy
  */
  double GetTotalEntropy(){ return totalentropy; };
  

   private:

   // the hdf5 file pointer, e.g. *.hdf5
   hid_t H5file_ptr_;

   // the hdf5/group pointer, e.g. /event0
   hid_t H5group_ptr_;

    //! Load saved configurations for each event
   void ReadConfigs();
 
   //! Load saved number of binary collisions
   void ReadNbcDist();

   //! Load saved initial entropy density distribution
   void ReadEntropyDist();

   //! want to use auxiliary hdf5 file readers 
    HydroinfoH5 * h5_helper_;

    int dim_x_, dim_y_;

    double npart=-1;
    double ncoll=-1;
    double totalentropy=-1;
  
    // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<InitialFromFile> reg;
};

#endif  // INITIALFROMFILE_H
