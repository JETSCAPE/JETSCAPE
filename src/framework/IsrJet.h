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

#ifndef ISRJET_H
#define ISRJET_H

//#include "JetScapeTask.h"
#include "JetEnergyLoss.h"
#include "sigslot.h"

namespace Jetscape {
/** @class Jet energy loss manager.
   */
class IsrJet
    : public JetEnergyLoss
      //public std::enable_shared_from_this<IsrManager>
      {

public:
  /** Default constructor to create a jet energy loss manager. Sets task ID as "JLossManager". Flag GetHardPartonListConnected is set to false.
   */
  IsrJet();

  /** Destructor for the jet energy loss manager.
   */
  virtual ~IsrJet();

  virtual void Init();


private:

};

} // end namespace Jetscape

#endif
