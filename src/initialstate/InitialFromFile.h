/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
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
#include "tinyxml2.h"
#include "InitialState.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

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

  /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in the XML file under the tag <IS>.
   */
  tinyxml2::XMLElement * GetIniStateXML() { return xml_; }

  // one can set range by hand if not read from xml file 
  /** Sets the range of the coordinates (xmax, ymax, zmax). 
      @param xmax Maximum value of the coordinate x in the nuclear density profile.
      @param ymax Maximum value of the coordinate y in the nuclear density profile.
      @param zmax Maxium value of the spatial rapidity ( if (tau,x,y,eta) system), or maximum value of the coordinate z (if in (t,x,y,z) system) in the nuclear density profile.  
   */

   private:
   // the hdf5 file pointer, e.g. *.hdf5
   hid_t H5file_ptr_;

   // the hdf5/group pointer, e.g. /event0
   hid_t H5group_ptr_;

    //! Load saved configurations for each event
   void read_configs_();
 
   //! Load saved number of binary collisions
   void read_nbc_dist_();

   //! Load saved initial entropy density distribution
   void read_entropy_dist_();

   //! want to use auxiliary hdf5 file readers 
    HydroinfoH5 * h5_helper_;

    int dim_x_, dim_y_;
};

#endif  // INITIALFROMFILE_H
