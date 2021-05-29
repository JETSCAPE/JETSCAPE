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

#ifndef ISRMANAGER_H
#define ISRMANAGER_H

//#include "JetScapeTask.h"
#include "JetEnergyLossManager.h"
#include "JetClass.h"
#include "sigslot.h"

#include <vector>

namespace Jetscape {
/** @class Jet energy loss manager.
   */
class IsrManager
    : public JetEnergyLossManager
      //public std::enable_shared_from_this<IsrManager>
      {

public:
  /** Default constructor to create a jet energy loss manager. Sets task ID as "JLossManager". Flag GetHardPartonListConnected is set to false.
   */
  IsrManager();

  /** Destructor for the jet energy loss manager.
   */
  virtual ~IsrManager();


  virtual void Init();


  virtual void Exec();

  /*
  virtual void Clear();

  virtual void CalculateTime();

  virtual void ExecTime();

  virtual void InitPerEvent();

  virtual void FinishPerEvent();

  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  */

private:

};

} // end namespace Jetscape

#endif
