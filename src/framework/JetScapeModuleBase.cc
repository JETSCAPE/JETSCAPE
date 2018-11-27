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

#include "JetScapeModuleBase.h"
#include "JetScapeXML.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"

#include<iostream>

namespace Jetscape {

  int JetScapeModuleBase::current_event=0;

  // ---------------------------------------------------------------------------
  /** Default constructor to create a JetScapeModuleBase. It sets the XML file name to a default string value.                                 
   */
  JetScapeModuleBase::JetScapeModuleBase()
  {
    xml_file_name = "";
  }

  // ---------------------------------------------------------------------------
  /** This is a constructor to create a JetScapeModuleBase. It sets the XML file name to "m_name" to be used to read input parameters.
   */
  JetScapeModuleBase::JetScapeModuleBase(string m_name)
  {
    xml_file_name = m_name;
  }

  // ---------------------------------------------------------------------------
  /** This is a destructor for the JetScapeModuleBase.                       
   */
  JetScapeModuleBase::~JetScapeModuleBase()
  {
    disconnect_all();
  }

  // ---------------------------------------------------------------------------
  /** A virtual function for a default initialization of a JetScapeModuleBase. It also checks whether a XML file is loaded or not.
   */
  void JetScapeModuleBase::Init()
  {
    if (!JetScapeXML::Instance()->GetXMLRoot())
      {
	JSWARN << "Not a valide JetScape XML file or no XML file loaded!";
	exit(-1);
      }
  }

  // ---------------------------------------------------------------------------
  /** This function returns a random number based on Mersenne-Twister algorithm.
   */
  shared_ptr<std::mt19937> JetScapeModuleBase::GetMt19937Generator(){
    // Instantiate if it isn't there yet
    if ( !mt19937_generator_ ){
      mt19937_generator_ = JetScapeTaskSupport::Instance()->GetMt19937Generator( GetMyTaskNumber() );
    }
    return mt19937_generator_;
  }
  

} // end namespace Jetscape
