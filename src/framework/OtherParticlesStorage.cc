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


#include "OtherParticlesStorage.h"
#include "JetScapeLogger.h"
#include <iostream>
#include <vector>


using namespace std;

namespace Jetscape {

  OtherParticlesStorage::OtherParticlesStorage()
  {
    SetId("OtherParticlesStorage");
    GetPhotonListConnected = false;
    VERBOSE(8); 
  }

  OtherParticlesStorage::~OtherParticlesStorage()
  {
    Clear();
  }

  void OtherParticlesStorage::Clear()
  {
    JSDEBUG << "Clearing other particles lists...";

    photons.clear();

    VERBOSE(8)<<photons.size();
  }

  void OtherParticlesStorage::Init()
  {
    JSINFO<<"Intialize Other Particles Storage...";
  }

  void OtherParticlesStorage::WriteTask(weak_ptr<JetScapeWriter> w)
  {
    VERBOSE(8);
  }

} // end of jetscape namespace 
