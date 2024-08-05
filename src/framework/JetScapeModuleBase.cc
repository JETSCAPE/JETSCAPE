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

#include "JetScapeModuleBase.h"

#include <iostream>

#include "JetScapeLogger.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeXML.h"

namespace Jetscape {

// Create an instance of the static map to register modules
JetScapeModuleFactory::map_type *JetScapeModuleFactory::moduleMap =
    new JetScapeModuleFactory::map_type;

int JetScapeModuleBase::current_event = 0;

// ---------------------------------------------------------------------------
/** Default constructor to create a JetScapeModuleBase. It sets the XML file
 * name to a default string value.
 */
JetScapeModuleBase::JetScapeModuleBase()
    : JetScapeTask(),
      xml_main_file_name(""),
      xml_user_file_name(""),
      mt19937_generator_(nullptr) {}

// ---------------------------------------------------------------------------
/** This is a destructor for the JetScapeModuleBase.
 */
JetScapeModuleBase::~JetScapeModuleBase() { disconnect_all(); }

// ---------------------------------------------------------------------------
/** A virtual function for a default initialization of a JetScapeModuleBase. It
 * also checks whether a XML file is loaded or not.
 */
void JetScapeModuleBase::Init() {
  if (!JetScapeXML::Instance()->GetXMLRootMain()) {
    JSWARN << "Not a valid JetScape Main XML file or no XML file loaded!";
    exit(-1);
  }
  if (!JetScapeXML::Instance()->GetXMLRootUser()) {
    JSWARN << "Not a valid JetScape XML file or no XML file loaded!";
    exit(-1);
  }
}

// ---------------------------------------------------------------------------
/** This function returns a random number based on Mersenne-Twister algorithm.
 */
shared_ptr<std::mt19937> JetScapeModuleBase::GetMt19937Generator() {
  // Instantiate if it isn't there yet
  if (!mt19937_generator_) {
    mt19937_generator_ =
        JetScapeTaskSupport::Instance()->GetMt19937Generator(GetMyTaskNumber());
  }
  return mt19937_generator_;
}

}  // end namespace Jetscape
